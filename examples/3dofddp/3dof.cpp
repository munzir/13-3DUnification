#include <dart/dart.hpp>
#include <dart/gui/gui.hpp>
#include <dart/utils/urdf/urdf.hpp>
#include <iostream>
#include <fstream>
#include <boost/circular_buffer.hpp>
#include <chrono>
#include <thread>
#include <ddp/costs.hpp>
#include <ddp/ddp.hpp>
#include <ddp/mpc.hpp>
#include <ddp/util.hpp>
#include "krang3d.h"

using namespace std;
using namespace dart::common;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::math;


const double default_speed_increment = 0.5;

const int default_ik_iterations = 4500;

const double default_force =  50.0; // N
const int default_countdown = 100;  // Number of timesteps for applying force


class Controller
{
public:
  /// Constructor
  Controller(const SkeletonPtr& skel)
    : mNDOF(skel),
      mSpeed(0.0)
  {
    int nDofs = mNDOF->getNumDofs();
    
    mForces = Eigen::VectorXd::Zero(nDofs);
    
    mKp = Eigen::MatrixXd::Identity(nDofs, nDofs);
    mKd = Eigen::MatrixXd::Identity(nDofs, nDofs);
  
    for(std::size_t i = 0; i < 6; ++i)
    {
      mKp(i, i) = 0.0;
      mKd(i, i) = 0.0;
    }

    for(std::size_t i = 6; i < mNDOF->getNumDofs(); ++i)
    {
      mKp(i, i) = 1000;
      mKd(i, i) = 50;
    }
    
    setTargetPositions(mNDOF->getPositions());
  }
  
  /// Reset the desired dof position to the current position
  void setTargetPositions(const Eigen::VectorXd& pose)
  {
    mTargetPositions = pose;
  }

  /// Clear commanding forces
  void clearForces()
  {
    mForces.setZero();
  }
  
  /// Add commanding forces from PD controllers (Lesson 2 Answer)
  void addPDForces()
  {
    Eigen::VectorXd q = mNDOF->getPositions();
    Eigen::VectorXd dq = mNDOF->getVelocities();
    
    Eigen::VectorXd p = -mKp * (q - mTargetPositions);
    Eigen::VectorXd d = -mKd * dq;
    
    mForces += p + d;
    mNDOF->setForces(mForces);
  }

  /// Add commanind forces from Stable-PD controllers (Lesson 3 Answer)
  void addSPDForces()
  {
    Eigen::VectorXd q = mNDOF->getPositions();
    Eigen::VectorXd dq = mNDOF->getVelocities();

    Eigen::MatrixXd invM = (mNDOF->getMassMatrix()
                            + mKd * mNDOF->getTimeStep()).inverse();
    Eigen::VectorXd p =
        -mKp * (q + dq * mNDOF->getTimeStep() - mTargetPositions);
    Eigen::VectorXd d = -mKd * dq;
    Eigen::VectorXd qddot =
        invM * (-mNDOF->getCoriolisAndGravityForces()
            + p + d + mNDOF->getConstraintForces());
    
    mForces += p + d - mKd * qddot * mNDOF->getTimeStep();
    mNDOF->setForces(mForces);
  }
  
  /// add commanding forces from ankle strategy (Lesson 4 Answer)
  void addAnkleStrategyForces()
  {
    Eigen::Vector3d COM = mNDOF->getCOM();
    // Approximated center of pressure in sagittal axis
    Eigen::Vector3d offset(0.05, 0, 0);
    Eigen::Vector3d COP = mNDOF->getBodyNode("h_heel_left")->
        getTransform() * offset;
    double diff = COM[0] - COP[0];

    Eigen::Vector3d dCOM = mNDOF->getCOMLinearVelocity();
    Eigen::Vector3d dCOP =  mNDOF->getBodyNode("h_heel_left")->
        getLinearVelocity(offset);
    double dDiff = dCOM[0] - dCOP[0];

    int lHeelIndex = mNDOF->getDof("j_heel_left_1")->getIndexInSkeleton();
    int rHeelIndex = mNDOF->getDof("j_heel_right_1")->getIndexInSkeleton();
    int lToeIndex = mNDOF->getDof("j_toe_left")->getIndexInSkeleton();
    int rToeIndex = mNDOF->getDof("j_toe_right")->getIndexInSkeleton();
    if(diff < 0.1 && diff >= 0.0) {
      // Feedback rule for recovering forward push
      double k1 = 200.0;
      double k2 = 100.0;
      double kd = 10;
      mForces[lHeelIndex] += -k1 * diff - kd * dDiff;
      mForces[lToeIndex] += -k2 * diff - kd * dDiff;
      mForces[rHeelIndex] += -k1 * diff - kd * dDiff;
      mForces[rToeIndex] += -k2 * diff - kd * dDiff;
    }else if(diff > -0.2 && diff < -0.05) {
      // Feedback rule for recovering backward push
      double k1 = 2000.0;
      double k2 = 100.0;
      double kd = 100;
      mForces[lHeelIndex] += -k1 * diff - kd * dDiff;
      mForces[lToeIndex] += -k2 * diff - kd * dDiff;
      mForces[rHeelIndex] += -k1 * diff - kd * dDiff;
      mForces[rToeIndex] += -k2 * diff - kd * dDiff;
    }  
    mNDOF->setForces(mForces);
  }

  void initialize_observer_gains()
  {

  }



  // Send velocity commands on wheel actuators (Lesson 6 Answer)
  void setWheelCommands()
  {
    // Extract Position and Angle
    // Eigen::VectorXd q = m3DOF->getPositions();

    // cout << q << endl;

    // Set Wheel Speed from LQR Controller

    int wheelFirstIndex = mNDOF->getDof("JLWheel")->getIndexInSkeleton();
    for (std::size_t i = wheelFirstIndex; i < mNDOF->getNumDofs(); ++i)
    {
      mKp(i, i) = 0.0;
      mKd(i, i) = 0.0;
    }
    
    int index1 = mNDOF->getDof("JLWheel")->getIndexInSkeleton();
    int index2 = mNDOF->getDof("JRWheel")->getIndexInSkeleton();

    mNDOF->setCommand(index1, mSpeed);
    mNDOF->setCommand(index2, mSpeed);

  }
  
  void changeWheelSpeed(double increment)
  {
    mSpeed += increment;
    std::cout << "wheel speed = " << mSpeed << std::endl;
  }
  
protected:
  /// The biped Skeleton that we will be controlling
  SkeletonPtr mNDOF;
  
  /// Joint forces for the biped (output of the Controller)
  Eigen::VectorXd mForces;
  
  /// Control gains for the proportional error terms in the PD controller
  Eigen::MatrixXd mKp;

  /// Control gains for the derivative error terms in the PD controller
  Eigen::MatrixXd mKd;

  /// Target positions for the PD controllers
  Eigen::VectorXd mTargetPositions;
    
  /// For velocity actuator: Current speed of the skateboard
  double mSpeed;

  /// Pendulum Control Parameters
};



class filter {
  public:
    filter(const int dim, const int n)
    {
      samples.set_capacity(n);
      total = Eigen::VectorXd::Zero(dim,1);
    }
    void AddSample(Eigen::VectorXd v)
    {
      if(samples.full()) 
      {
        total -= samples.front();
      }
      samples.push_back(v);
      total += v;
      average = total/samples.size();
    }
  
    boost::circular_buffer<Eigen::VectorXd> samples;
    Eigen::VectorXd total;
    Eigen::VectorXd average;
    
};

class MyWindow : public dart::gui::SimWindow
{
    using Dynamics = Krang3D<double>;
    using Scalar = double;
    using DDP_Opt = optimizer::DDP<Dynamics>;
    using Cost = Krang3DCost<Scalar>;
    using TerminalCost = Krang3DTerminalCost<Scalar>;
    using StateTrajectory = typename Dynamics::StateTrajectory ;
    using ControlTrajectory= typename Dynamics::ControlTrajectory ;
    using State = typename Dynamics::State;
    using Control = typename Dynamics::Control;

  public:
    MyWindow(const WorldPtr& world)
    {

      setWorld(world);
      m3DOF = world->getSkeleton("m3DOF");
      qInit = m3DOF->getPositions();
      
      psi = 0; // Heading Angle

      R = 0.25;
      L = 0.68;//*6;


      steps = 0;

      outFile.open("constraints.csv");
      dqFilt = new filter(8, 100);
      cFilt = new filter(5, 100);

      mController = dart::common::make_unique<Controller>(m3DOF);

      util::DefaultLogger logger;
      Scalar tf = 20;
      Scalar ddp_dt = 0.01;
      auto time_steps = util::time_steps(tf, ddp_dt);
      int max_iterations = 15;
      bool verbose = true;

      param p;
      p.R = 2.500000e-01; p.mw = 0.5; p.Iw = 5.100000e-03; p.L = 0.8; p.g=9.800000e+00;
      p.m_1 = 115.0;
      p.MX_1 = 0; p.MY_1 = 0.8; p.MZ_1 = 0;
      p.XX_1 = 26; p.YY_1 = 0.0832724; p.ZZ_1 = 0.086493;
      p.XY_1 = 2.45462e-05; p.YZ_1 = -0.00131733; p.XZ_1 = 0.00713022;
      p.fric_1 = 15;
      p.XXw = 0.005; p.YYw=0.0025; p.ZZw=0.0025;

      Dynamics krangDynamics(p);
      Dynamics::State x0 = Dynamics::State::Zero();
      Dynamics::State xf; xf << 0, 0, 0, 0, 0, 0, 5, 0;
      Dynamics::ControlTrajectory u = Dynamics::ControlTrajectory::Zero(2, time_steps);
      Cost::StateHessian Q;
      Q.setZero();
      Q.diagonal() << 0,0.1,0.1,0.1,0.1,0.1,0.1,0.1;

      Cost::ControlHessian R;
      R.setZero();
      R.diagonal() << 0.01, 0.01;

      TerminalCost::Hessian Qf;
      Qf.setZero();
      Qf.diagonal() << 0,1e4,1e4,1e4,1e4,1e4,1e4,1e4;

      Cost cp_cost(xf, Q, R);
      TerminalCost cp_terminal_cost(xf, Qf);

      // initialize DDP for trajectory planning
      DDP_Opt trej_ddp (ddp_dt, time_steps, max_iterations, &logger, verbose);

      // Get initial trajectory from DDP
      OptimizerResult<Dynamics> traj_results = trej_ddp.run(x0, u, krangDynamics, cp_cost, cp_terminal_cost);
      StateTrajectory xs = traj_results.state_trajectory;
      ctl_traj = traj_results.control_trajectory;
      state_traj = traj_results.state_trajectory;

      std::cout << "Obtained initial state trajectory";
      for (int m = 0; m < xs.cols(); ++m) {
        std::cout << "\n";
        for (int n = 0; n < xs.rows(); ++n) {
          std::cout << xs(n, m);
        }
      }

      std::cout << "Obtained initial control trajectory";
      for (int m = 0; m < ctl_traj.cols(); ++m) {
        std::cout << "\n";
        for (int n = 0; n < ctl_traj.rows(); ++n) {
          std::cout << ctl_traj(n, m);
        }
      }
    }

    void timeStepping() override
    {
      // Read Positions, Speeds, Transform speeds to world coordinates and filter the speeds
//      Eigen::Matrix<double, 4, 4> Tf = m3DOF->getBodyNode(0)->getTransform().matrix();
//      psi =  atan2(Tf(0,0),-Tf(1,0));
//      qBody1 = atan2(Tf(0,1)*cos(psi) + Tf(1,1)*sin(psi), Tf(2,1));
      Eigen::VectorXd q = m3DOF->getPositions();
//      Eigen::VectorXd xPlane(3);
//      xPlane << q(0),q(1),0;
//
//      // Distance from Origin
//      dist = xPlane.norm();
//
//      Eigen::VectorXd dq_orig = m3DOF->getVelocities();
//      Eigen::Matrix<double, 8, 1> dq;
//      dq << (Tf.block<3,3>(0,0) * dq_orig.head(3)) , (Tf.block<3,3>(0,0) * dq_orig.segment(3,3)), dq_orig(6), dq_orig(7);
//      dqFilt->AddSample(dq);
//
//      // Wheel Rotation (World Frame)
//      thL = q(6) - std::abs(qBody1);
//      thR = q(7) - std::abs(qBody1);
//
//      // Corresponding Velocities (Filtered)
//      dpsi = dq(2);
//      dpsiFilt = dqFilt->average(2);
//      dqBody1 = -dq_orig(0);
//      dqBody1Filt = (-dqFilt->average(0)*sin(psi) + dqFilt->average(1)*cos(psi));
//      dthL = dq(6) + dqBody1;
//      dthLFilt = dqFilt->average(6) + dqBody1Filt;
//      dthR = dq(7) + dqBody1;
//      dthRFilt = dqFilt->average(7) + dqBody1Filt;
//
//      // thR =
//      // Constraints
//      // 1. dZ0 = 0                                               => dq_orig(4)*cos(qBody1) + dq_orig(5)*sin(qBody1) = 0
//      // 2. da3 + R/L*(dthL - dthR) = 0                           => dq_orig(1)*cos(qBody1) + dq_orig(2)*sin(qBody1) + R/L*(dq_orig(6) - dq_orig(7)) = 0
//      // 3. da1*cos(psii) + da2*sin(psii) = 0                     => dq_orig(1)*sin(qBody1) - dq_orig(2)*cos(qBody1) = 0
//      // 4. dX0*sin(psii) - dY0*cos(psii) = 0                     => dq_orig(3) = 0
//      // 5. dX0*cos(psii) + dY0*sin(psii) - R/2*(dthL + dthR) = 0 => dq_orig(4)*sin(qBody1) - dq_orig(5)*cos(qBody1) - R/2*(dq_orig(6) + dq_orig(7) - 2*dq_orig(0)) = 0
//      Eigen::Matrix<double, 5, 1> c;
//      c << (dq_orig(4)*cos(qBody1) + dq_orig(5)*sin(qBody1)), (dq_orig(1)*cos(qBody1) + dq_orig(2)*sin(qBody1) + R/L*(dq_orig(6) - dq_orig(7))), (dq_orig(1)*sin(qBody1) - dq_orig(2)*cos(qBody1)), dq_orig(3), (dq_orig(4)*sin(qBody1) - dq_orig(5)*cos(qBody1) - R/2*(dq_orig(6) + dq_orig(7) - 2*dq_orig(0)));
//      cFilt->AddSample(c);
//      //if(Tf(2,1) > 0)
//      {
//      for(int i=0; i<8; i++) outFile << dq(i) << ", ";
//      for(int i=0; i<8; i++) outFile << dqFilt->average(i) << ", ";
//      outFile << psi << ", " << dpsi << ", " << qBody1 << ", " << dqBody1 << ", " <<  dthL << ", " << dthR << ", ";
//      outFile << psiFilt << ", " << dpsiFilt << ", " << qBody1Filt << ", " << dqBody1Filt << ", " <<  dthLFilt << ", " << dthRFilt << ", ";
//      for(int i=0; i<5; i++) outFile << c(i) << ", ";
//      for(int i=0; i<5; i++) outFile << cFilt->average(i) << ", ";
//      for(int i=0; i<8; i++) outFile << q(i) << ", ";
//      outFile << std::endl;
//      }
//
//
//      // Wheeled Inverted Pendulum Parameters
//      double I_ra = 0;
//      double gamma = 1.0;
//      double g = 9.81;
//      double c_w = 0.1;
//
//      double r_w, m_w, I_wa; // Wheel Parameters
//      double M_g, l_g; // COM Parameters
//
//      double Iw_xx, Iw_yy, Iw_zz, Iw_xy, Iw_xz, Iw_yz; // Wheel Moment of Inertia
//      double I_xx, I_yy, I_zz, I_xy, I_xz, I_yz; // Base Moment of Inertia
//
//      double delta, c1, c2; // Intermediate Parameters
//
//      double u_theta, u_psi, tau_w; // Control Signals
//
//      Eigen::VectorXd x_dot(10);
//      Eigen::VectorXd x_dot_sys(4); // Observer Dynamics
//      Eigen::VectorXd x_dot_w(3);
//      Eigen::VectorXd x_dot_p(3);
//
//      std::vector< BodyNode * > nodes = m3DOF->getBodyNodes();
//      for(auto const& n: nodes) {
//        std::string node_name = n->getName();
//
//        // Extract Wheel Parameters
//        if(node_name.find("Wheel") != string::npos){
//          m_w = n->getMass();
//          n->getMomentOfInertia(Iw_xx,Iw_yy,Iw_zz,Iw_xy,Iw_xz,Iw_yz);
//          I_wa = Iw_xx;
//          r_w = 0.25;
//        }
//
//        // Extract Body Parameters
//        if(node_name.find("Base") != string::npos){
//          M_g = n->getMass();
//          n->getMomentOfInertia(I_xx,I_yy,I_zz,I_xy,I_xz,I_yz);
//          Eigen::Vector3d com = n->getLocalCOM();
//          l_g = com[1];
//        }
//      }
//
//      delta = (M_g*l_g+I_yy+pow(gamma,2)*I_ra)*(M_g+m_w)*pow(r_w,2)+I_wa+I_ra*pow(gamma,2)-pow(M_g*r_w*l_g-I_ra*pow(gamma,2),2);
//      c1 = (M_g+m_w)*pow(r_w,2)+I_wa+I_ra*pow(gamma,2)+M_g*r_w*l_g+I_ra*pow(gamma,2);
//      c2 = M_g*r_w*l_g+M_g*pow(l_g,2)+I_yy;
//
//      // Hardcode Dynamics + Feedback for Now
//      Eigen::MatrixXd A(4,4);
//      Eigen::MatrixXd B(4,1);
//
//      Eigen::MatrixXd A_(3,3);
//      Eigen::MatrixXd B_1(3,1);
//      Eigen::MatrixXd B_2(3,1);
//      Eigen::MatrixXd L_1(3,1);
//      Eigen::MatrixXd L_2(3,1);
//      // Eigen::MatrixXi C_(4,4);
//
//      A << 0, 0, 1, 0,
//           0, 0, 0, 1,
//           17.7828531704201,  0, -0.00858417221300891,  0.00858417221300891,
//           47.8688365622367,  0, 0.0288387521185202,  -0.0288387521185202;
//
//      B << 0,
//           0,
//           -0.0858417221300890,
//           0.288387521185202;
//
//      A_ << 0, 1, 0,
//            0, 0, 1,
//            0, 0, 0;
//
//      B_1 << 0,
//             0.288387521185202,
//             0;
//
//      B_2 << 0,
//             -0.0858417221300890,
//             0;
//
//      L_1 << 1159.99999999673,173438.396407957,1343839.4084839;
//      L_2 = L_1;
//
//
//      Eigen::VectorXd xL(4);
//      xL << qBody1, thR, dqBody1Filt, dthRFilt;
//
//      Eigen::VectorXd xR(4);
//      xR << qBody1, thR, dqBody1Filt, dthRFilt;
//
//      if(steps==0){
//        thR_hat = thR;
//        dthR_hat = dthRFilt;
//        f_thR = 0.0;
//
//        qBody1_hat = qBody1;
//        dqBody1_hat = dqBody1Filt;
//        f_qBody1 = 0.0;
//      }
//
//      // double thR_hat, dthR_hat, f_thR, qBody1_hat, dqBody1_hat, f_qBody1;
//      Eigen::VectorXd xFull(10);
//      xFull << qBody1, thR, dqBody1Filt, dthRFilt, thR_hat, dthR_hat, f_thR, qBody1_hat, dqBody1_hat, f_qBody1;
//
//      Eigen::VectorXd xDesired(4);
//      xDesired << 0,0,0,0;
//
//      Eigen::VectorXd F(4);
//
//      // Break Dance Moves
//      F << -510.449550243191,  -0.244948974282393,  -110.528023245082, -1.14116859943037;
//
//
//      tauL = -0.5*F.transpose()*(xL-xDesired);
//      tauR = -0.5*F.transpose()*(xR-xDesired);
//
//      // Observer Control
//      Eigen::VectorXd ut1(2); ut1 << F[1], F[3];
//      Eigen::VectorXd ut2(2); ut2 << xFull[4]-xDesired[1], xFull[5]-xDesired[3];
//
//      Eigen::VectorXd up1(2); up1 << F[0], F[2];
//      Eigen::VectorXd up2(2); up2 << xFull[7]-xDesired[0], xFull[8]-xDesired[2];
//
//      u_theta = -ut1.transpose()*ut2;
//      u_psi = -up1.transpose()*up2;
//      tau_w = u_psi + u_theta - xFull[9] - xFull[6];
//
//      // Real Dynamics
//      // SIMULATOR HERE
//
//      // Observer Dynamics
//      // x_dot_sys = x
//
//
//      Eigen::VectorXd wh(3); wh << xFull[4], xFull[5], xFull[6];
//
//      Eigen::VectorXd ph(3); ph << xFull[7], xFull[8], xFull[9];
//      // x_dot_sys = Ignore, Just Use Simulation Output
//      x_dot_w = A_*wh + L_1*(xFull[1]-xFull[4]) + B_1*u_theta;
//      x_dot_p = A_*ph + L_2*(xFull[0]-xFull[7]) + B_2*u_psi;
//      x_dot << x_dot_sys, x_dot_w, x_dot_p;
//
//      tauL = 0.5*tau_w;
//      tauR = 0.5*tau_w;
//
//       Apply Control Input

      double dt = m3DOF->getTimeStep();

      // TODO: CHANGE HARD CODED DT DIFFERENCE
      double ddp_step = steps / 10;
      State cur_x = state_traj.col(ddp_step);
      State next_x = state_traj.col(ddp_step + 1);
      Control cur_u = ctl_traj.col(ddp_step);

      param p;
      p.R = 2.500000e-01; p.mw = 0.5; p.Iw = 5.100000e-03; p.L = 0.8; p.g=9.800000e+00;
      p.m_1 = 115.0;
      p.MX_1 = 0; p.MY_1 = 0.8; p.MZ_1 = 0;
      p.XX_1 = 26; p.YY_1 = 0.0832724; p.ZZ_1 = 0.086493;
      p.XY_1 = 2.45462e-05; p.YZ_1 = -0.00131733; p.XZ_1 = 0.00713022;
      p.fric_1 = 15;
      p.XXw = 0.005; p.YYw=0.0025; p.ZZw=0.0025;
      Dynamics krangDynamics(p);

      c_forces dyn_forces = krangDynamics.dynamic_forces(cur_x, cur_u);

      double dd_x = (next_x(3) - cur_x(3)) / dt;
      double dd_psi = (next_x(4) - cur_x(4)) / dt;
      double ddth = cur_u(0);

      double tau_0 =  cur_u(1);
      double tau_1 = dyn_forces.A(2, 0) * dd_x + dyn_forces.A(2, 1) * dd_psi + dyn_forces.A(2, 2) * ddth + dyn_forces.C(2, 0) * cur_x(3) + dyn_forces.C(2, 1) * cur_x(4) + dyn_forces.C(2, 2) + dyn_forces.Q(2) + dyn_forces.Gamma_fric(2);
      double tau_l = (tau_1 + tau_0) / 2;
      double tau_r = (tau_1 - tau_0) / 2;

      mForces << 0, 0, 0, 0, 0, 0, tau_l, tau_r;
      m3DOF->setForces(mForces);
      steps++;



//      printf("%6d ",steps);
//      printf("%12.4f ",tauR);
//      for(int i = 0; i < 10; i++){printf("%12.7f ",xFull[i]);}
//      printf("\n");

      // Integrate Dynamics
//      xFull = xFull + dt*x_dot;
//
//      thR_hat = xFull[4];
//      dthR_hat = xFull[5];
//      f_thR = xFull[6];
//      qBody1_hat = xFull[7];
//      dqBody1_hat = xFull[8];
//      f_qBody1 = xFull[9];

      SimWindow::timeStepping();
    }
    ~MyWindow() {
      outFile.close();     
    }
    

  protected:

    SkeletonPtr m3DOF;

    Eigen::VectorXd qInit;

    Eigen::VectorXd dof1;

    std::unique_ptr<Controller> mController;

    double dist;
    double psi, dpsi, thL, dthL, thR, dthR, qBody1, dqBody1; // State Variables
    double thR_hat, dthR_hat, f_thR, qBody1_hat, dqBody1_hat, f_qBody1; // Observer Variables
    double psiFilt, dpsiFilt, qBody1Filt, dqBody1Filt, dthLFilt, dthRFilt; // Filtered Variables

    double L, R;
    double tauL, tauR;
    
    int steps;
    double ddp_dt;

    Eigen::Matrix<double, 8, 1> mForces;
   
    ofstream outFile; 

    filter *dqFilt, *cFilt;
    ControlTrajectory ctl_traj;
    StateTrajectory state_traj;
};


SkeletonPtr createFloor()
{
  SkeletonPtr floor = Skeleton::create("floor");

  // Give the floor a body
  BodyNodePtr body =
      floor->createJointAndBodyNodePair<WeldJoint>(nullptr).second;
//  body->setFrictionCoeff(1e16);

  // Give the body a shape
  double floor_width = 50;
  double floor_height = 0.05;
  std::shared_ptr<BoxShape> box(
        new BoxShape(Eigen::Vector3d(floor_width, floor_width, floor_height)));
  auto shapeNode
      = body->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
  shapeNode->getVisualAspect()->setColor(dart::Color::Blue());

  // Put the body into position
  Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
  tf.translation() = Eigen::Vector3d(0.0, 0.0, -floor_height / 2.0);
  body->getParentJoint()->setTransformFromParentBodyNode(tf);

  return floor;
}


SkeletonPtr create3DOF_URDF()
{
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  SkeletonPtr threeDOF =
      loader.parseSkeleton("/Users/BakerStreetBakery/Desktop/Krang/10d-3DHighLevelControllerC-2Ddart/examples/3dofddp/3dof.urdf");
  threeDOF->setName("m3DOF");

  // Get it into a useful configuration
  double psiInit = M_PI/4, qBody1Init = 0.01*M_PI;//M_PI;//0;
  Eigen::Transform<double, 3, Eigen::Affine> baseTf = Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  // RotX(pi/2)*RotY(-pi/2+psi)*RotX(-qBody1)
  baseTf.prerotate(Eigen::AngleAxisd(-qBody1Init,Eigen::Vector3d::UnitX())).prerotate(Eigen::AngleAxisd(-M_PI/2+psiInit,Eigen::Vector3d::UnitY())).prerotate(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()));
  Eigen::AngleAxisd aa(baseTf.matrix().block<3,3>(0,0));
  Eigen::Matrix<double, 8, 1> q;
//  q << 1.2092, -1.2092, -1.2092, 0, 0, 0.28, 0, 0;
  q << aa.angle()*aa.axis(), 0, 0, 0.28, 0, 0;
  threeDOF->setPositions(q);

  return threeDOF;
}


int main(int argc, char* argv[])
{

  SkeletonPtr threeDOF = create3DOF_URDF();
  SkeletonPtr floor = createFloor();

  WorldPtr world = std::make_shared<World>();
  world->addSkeleton(threeDOF);
  world->addSkeleton(floor);

  MyWindow window(world);
  glutInit(&argc, argv);
  window.initWindow(1280,720, "3DOF URDF");
  glutMainLoop();
}


