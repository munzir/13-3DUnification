#include <dart/dart.hpp>
#include <dart/gui/gui.hpp>
#include <dart/utils/urdf/urdf.hpp>
#include <iostream>
#include <fstream>
#include <boost/circular_buffer.hpp>
#include <ddp/costs.hpp>
#include <ddp/ddp.hpp>
#include <ddp/mpc.hpp>
#include <ddp/util.hpp>
#include <nlopt.hpp>
#include "krangddp.h"

using namespace std;
using namespace dart::common;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::math;

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

Eigen::Vector3d getBodyCOM(dart::dynamics::SkeletonPtr robot) {
  double fullMass = robot->getMass();
  double wheelMass = robot->getBodyNode("LWheel")->getMass();
  return (fullMass*robot->getCOM() - wheelMass*robot->getBodyNode("LWheel")->getCOM() - wheelMass*robot->getBodyNode("RWheel")->getCOM())/(fullMass - 2*wheelMass);
}

void getSimple(SkeletonPtr& threeDOF, SkeletonPtr& krang) 
{
  // Load the full body with fixed wheel and set the pose q
  // dart::utils::DartLoader loader;
  // SkeletonPtr krangFixedWheel =
  //     loader.parseSkeleton("/home/krang/dart/09-URDF/KrangFixedWheels/krang_fixed_wheel.urdf");
  
  // Body Mass
  double mFull = krang->getMass(); 
  double mLWheel = krang->getBodyNode("LWheel")->getMass();
  double mRWheel = krang->getBodyNode("RWheel")->getMass();
  double mBody = mFull - mLWheel - mRWheel;


  Eigen::Vector3d bodyCOM;
  dart::dynamics::Frame* baseFrame = krang->getBodyNode("Base");
  bodyCOM = (mFull*krang->getCOM(baseFrame) - mLWheel*krang->getBodyNode("LWheel")->getCOM(baseFrame) - mLWheel*krang->getBodyNode("RWheel")->getCOM(baseFrame))/(mFull - mLWheel - mRWheel);

  // Body inertia (axis)
  double m;
  Eigen::Matrix3d iMat;
  Eigen::Matrix3d iBody = Eigen::Matrix3d::Zero();
  double ixx, iyy, izz, ixy, ixz, iyz;  
  Eigen::Matrix3d rot;
  Eigen::Vector3d t;
  Eigen::Matrix3d tMat;
  dart::dynamics::BodyNodePtr b;
  int nBodies = krang->getNumBodyNodes();
  for(int i=0; i<nBodies; i++){
    if(i==1 || i==2) continue; // Skip wheels
    b = krang->getBodyNode(i);
    b->getMomentOfInertia(ixx, iyy, izz, ixy, ixz, iyz);
    rot = b->getTransform(baseFrame).rotation();
    t = krang->getCOM(baseFrame) - b->getCOM(baseFrame) ; // Position vector from local COM to body COM expressed in base frame
    m = b->getMass();
    iMat << ixx, ixy, ixz, // Inertia tensor of the body around its CoM expressed in body frame
            ixy, iyy, iyz,
            ixz, iyz, izz;
    iMat = rot*iMat*rot.transpose(); // Inertia tensor of the body around its CoM expressed in base frame
    tMat << (t(1)*t(1)+t(2)*t(2)), (-t(0)*t(1)),          (-t(0)*t(2)),
            (-t(0)*t(1)),          (t(0)*t(0)+t(2)*t(2)), (-t(1)*t(2)),
            (-t(0)*t(2)),          (-t(1)*t(2)),          (t(0)*t(0)+t(1)*t(1));
    iMat = iMat + m*tMat; // Parallel Axis Theorem
    iBody += iMat;
  }


  // Aligning threeDOF base frame to have the y-axis pass through the CoM
  double th = atan2(bodyCOM(2), bodyCOM(1));
  rot << 1, 0, 0,
         0, cos(th), sin(th),
         0, -sin(th), cos(th);
  bodyCOM = rot*bodyCOM;
  iBody = rot*iBody*rot.transpose();

  // Set the 3 DOF robot parameters
  threeDOF->getBodyNode("Base")->setMomentOfInertia(iBody(0,0), iBody(1,1), iBody(2,2), iBody(0,1), iBody(0,2), iBody(1,2));
  threeDOF->getBodyNode("Base")->setLocalCOM(bodyCOM);
  threeDOF->getBodyNode("Base")->setMass(mBody);

  // Print them out
  cout << "mass: " << mBody << endl;
  cout << "COM: " << bodyCOM(0) << ", " << bodyCOM(1) << ", " << bodyCOM(2) << endl;
  cout << "ixx: " << iBody(0,0) << ", iyy: " << iBody(1,1) << ", izz: " << iBody(2,2) << endl;
  cout << "ixy: " << iBody(0,1) << ", ixz: " << iBody(0,2) << ", iyz: " << iBody(1,2) << endl;

  // Update 3DOF state
  // get positions
  Eigen::Matrix3d baseRot = krang->getBodyNode("Base")->getTransform().rotation();
  baseRot = baseRot*rot.transpose();
  Eigen::AngleAxisd aa(baseRot);
  Eigen::Matrix<double, 8, 1> q, dq;
  q << aa.angle()*aa.axis(), krang->getPositions().segment(3, 5);
  threeDOF->setPositions(q);

  // TODO: When joints are unlocked qBody1 of the 3DOF (= dth = COM angular speed) is not the same as qBody1 of the full robot
  dq << krang->getVelocities().head(8);
  threeDOF->setVelocities(dq);
}


SkeletonPtr create3DOF_URDF(SkeletonPtr krang)
{
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  SkeletonPtr threeDOF = 
      loader.parseSkeleton("/home/krang/dart/09-URDF/3DOF-WIP/3dof.urdf");
  threeDOF->setName("m3DOF");

  // Set parameters of Body that reflect the ones we will actually have 
  // Eigen::Matrix<double, 19, 1> qInit;
  // qInit << -M_PI/4, -4.588, 0.0, 0.0, 0.0548, -1.0253, 0.0, -2.1244, -1.0472, 1.5671, 0.0, -0.0548, 1.0253, 0.0, 2.1244, 1.0472, 0.0037, 0.0;
  // qInit << krang->getPositions();
//   getSimple(threeDOF, krang);   
  
  threeDOF->getJoint(0)->setDampingCoefficient(0, 0.5);
  threeDOF->getJoint(1)->setDampingCoefficient(0, 0.5);

//   // Get it into a useful configuration
//   double psiInit = 0, qBody1Init = 0;
//   Eigen::Transform<double, 3, Eigen::Affine> baseTf = Eigen::Transform<double, 3, Eigen::Affine>::Identity();
//   // RotX(pi/2)*RotY(-pi/2+psi)*RotX(-qBody1)
//   baseTf.prerotate(Eigen::AngleAxisd(-qBody1Init,Eigen::Vector3d::UnitX())).prerotate(Eigen::AngleAxisd(-M_PI/2+psiInit,Eigen::Vector3d::UnitY())).prerotate(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()));
//   Eigen::AngleAxisd aa(baseTf.matrix().block<3,3>(0,0));
//   Eigen::Matrix<double, 8, 1> q;
// //  q << 1.2092, -1.2092, -1.2092, 0, 0, 0.28, 0, 0;
//   q << aa.angle()*aa.axis(), 0, 0, 0.28, 0, 0;
//   threeDOF->setPositions(q);

  return threeDOF;
}



class MyWindow : public dart::gui::SimWindow
{
    using Scalar = double;  
    using Dynamics = Krang3D<Scalar>;
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
      cout << "8" << endl;
  

      setWorld(world);
      mkrang = world->getSkeleton("krang");
      m3DOF = create3DOF_URDF(mkrang);
      world3dof = std::make_shared<World>();
      world3dof->addSkeleton(m3DOF);
      getSimple(m3DOF, mkrang); 
      qInit = m3DOF->getPositions();
      // double psiInit = 0, qBody1Init = 0;
      // Eigen::Transform<double, 3, Eigen::Affine> baseTf = Eigen::Transform<double, 3, Eigen::Affine>::Identity();
      // // RotX(pi/2)*RotY(-pi/2+psi)*RotX(-qBody1)
      // baseTf.prerotate(Eigen::AngleAxisd(-qBody1Init,Eigen::Vector3d::UnitX())).prerotate(Eigen::AngleAxisd(-M_PI/2+psiInit,Eigen::Vector3d::UnitY())).prerotate(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()));
      // Eigen::AngleAxisd aa(baseTf.matrix().block<3,3>(0,0));
      // qInit << aa.angle()*aa.axis(), 0, 0, 0.28, 0, 0;
      // m3DOF->setPositions(qInit);


      psi = 0; // Heading Angle
      steps = 0;
      mpc_steps = -1; 
      mpc_dt = 0.01;
      outFile.open("constraints.csv");
      dqFilt = new filter(8, 50);
      cFilt = new filter(5, 50);
      R = 0.25;
      L = 0.68;//*6;
      mpc_writer.open_file("mpc_traj.csv");
      cout << "9" << endl;

      computeDDPTrajectory();
      cout << "10" << endl;

      
    }

    void computeDDPTrajectory() {
      cout << "9a" << endl;

      param p; 
      double ixx, iyy, izz, ixy, ixz, iyz; 
      Eigen::Vector3d com;
      Eigen::Matrix3d iMat;      
      Eigen::Matrix3d tMat;
            cout << "9b" << endl;

      dart::dynamics::Frame* baseFrame = m3DOF->getBodyNode("Base");
      p.R = 2.500000e-01; p.L = 6.000000e-01; p.g=9.800000e+00;
            cout << "9c" << endl;

      p.mw = m3DOF->getBodyNode("LWheel")->getMass(); 
            cout << "9d" << endl;

      m3DOF->getBodyNode("LWheel")->getMomentOfInertia(ixx, iyy, izz, ixy, ixz, iyz);
            cout << "9e" << endl;

      p.YYw = ixx; p.ZZw = izz; p.XXw = iyy; // Wheel frame of reference in ddp dynamic model is different from the one in DART
      p.m_1 = m3DOF->getBodyNode("Base")->getMass(); 
      com = m3DOF->getBodyNode("Base")->getCOM(baseFrame);
      p.MX_1 = p.m_1*com(0); p.MY_1 = p.m_1*com(1); p.MZ_1 = p.m_1*com(2);
            cout << "9f" << endl;

      m3DOF->getBodyNode("Base")->getMomentOfInertia(ixx, iyy, izz, ixy, ixz, iyz);
      Eigen::Vector3d s = -com; // Position vector from local COM to body COM expressed in base frame
      iMat << ixx, ixy, ixz, // Inertia tensor of the body around its CoM expressed in body frame
              ixy, iyy, iyz,
              ixz, iyz, izz;
      tMat << (s(1)*s(1)+s(2)*s(2)), (-s(0)*s(1)),          (-s(0)*s(2)),
              (-s(0)*s(1)),          (s(0)*s(0)+s(2)*s(2)), (-s(1)*s(2)),
              (-s(0)*s(2)),          (-s(1)*s(2)),          (s(0)*s(0)+s(1)*s(1));
      iMat = iMat + p.m_1*tMat; // Parallel Axis Theorem
      p.XX_1 = iMat(0,0); p.YY_1 = iMat(1,1); p.ZZ_1 = iMat(2,2);
      p.XY_1 = iMat(0,1); p.YZ_1 = iMat(1,2); p.XZ_1 = iMat(0,2);
      p.fric_1 = m3DOF->getJoint(0)->getDampingCoefficient(0); // Assuming both joints have same friction coeff (Please make sure that is true)
            cout << "9g" << endl;

      CSV_writer<Scalar> writer;
      util::DefaultLogger logger;
      bool verbose = true;
      Scalar tf = 20;
      auto time_steps = util::time_steps(tf, mpc_dt);
      int max_iterations = 15;
      

      ddp_dyn = new Dynamics(p);
       // Dynamics ddp_dyn(p);

      // Initial state th, dth, x, dx, desired state, initial control sequence

      State x0 = getCurrentState();
      x0 << 0, 0, 0, 0, 0, 0, 0, 0; 
      cout << x0 << endl;
      Dynamics::State xf; xf << 2, 0, 0, 0, 0, 0, 0.01, 5;
      Dynamics::ControlTrajectory u = Dynamics::ControlTrajectory::Zero(2, time_steps);

      // Costs
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
      DDP_Opt trej_ddp (mpc_dt, time_steps, max_iterations, &logger, verbose);

      cout << "9h" << endl;
      // Get initial trajectory from DDP
      OptimizerResult<Dynamics> DDP_traj = trej_ddp.run(x0, u, *ddp_dyn, cp_cost, cp_terminal_cost);

      cout << "9i" << endl;

      ddp_state_traj = DDP_traj.state_trajectory;
      ddp_ctl_traj = DDP_traj.control_trajectory;

      writer.save_trajectory(ddp_state_traj, ddp_ctl_traj, "initial_traj.csv");

    }

    State getCurrentState() {
      // Read Positions, Speeds, Transform speeds to world coordinates and filter the speeds      
      Eigen::Matrix<double, 4, 4> Tf = m3DOF->getBodyNode(0)->getTransform().matrix();
      psi =  atan2(Tf(0,0), -Tf(1,0));
      qBody1 = atan2(Tf(0,1)*cos(psi) + Tf(1,1)*sin(psi), Tf(2,1));
      Eigen::VectorXd q = m3DOF->getPositions();
      Eigen::VectorXd dq_orig = m3DOF->getVelocities();
      Eigen::Matrix<double, 8, 1> dq;
      dq << (Tf.block<3,3>(0,0) * dq_orig.head(3)) , (Tf.block<3,3>(0,0) * dq_orig.segment(3,3)), dq_orig(6), dq_orig(7);
      dqFilt->AddSample(dq);

      // Calculate the quantities we are interested in
      dpsi = dq(2);
      dpsiFilt = dqFilt->average(2);
      dqBody1 = -dq_orig(0);
      dqBody1Filt = (-dqFilt->average(0)*sin(psi) + dqFilt->average(1)*cos(psi));
      double thL = q(6) + qBody1;
      dthL = dq(6) + dqBody1;
      dthLFilt = dqFilt->average(6) + dqBody1Filt;
      double thR = q(7) + qBody1;
      dthR = dq(7) + dqBody1;
      dthRFilt = dqFilt->average(7) + dqBody1Filt;


      // State: x, psi, theta, dx, dpsi, dtheta, x0, y0
      State cur_state = Dynamics::State::Zero();
      cur_state << R/2 * (thL + thR), psi, qBody1, dq(3) * cos(psi) + dq(4) * sin(psi), dpsi, dqBody1, q(3), q(4); 
      return cur_state;
    }

    void timeStepping() override
    {
      steps++;

      getSimple(m3DOF, mkrang);
      State cur_state = getCurrentState(); 

      // MPC DDP RECEDING HORIZON CALCULATION
      int beginStep = 500; 
      int cur_mpc_steps = ((steps > beginStep) ? ((steps - beginStep) / 10) : -1);

      if (cur_mpc_steps > mpc_steps) {
        mpc_steps = cur_mpc_steps;
        int max_iterations = 15; 
        bool verbose = true; 
        util::DefaultLogger logger;
        int mpc_horizon = 10; 
        
        Dynamics::State target_state;
        target_state = ddp_state_traj.col(mpc_steps + mpc_horizon);
        Dynamics::ControlTrajectory hor_control = Dynamics::ControlTrajectory::Zero(2, mpc_horizon);
        Dynamics::StateTrajectory hor_traj_states = ddp_state_traj.block(0, mpc_steps, 8, mpc_horizon);
        
        DDP_Opt ddp_horizon (mpc_dt, mpc_horizon, max_iterations, &logger, verbose);
        
        Cost::StateHessian Q_mpc, Qf_mpc;
        Cost::ControlHessian ctl_R;
        
        ctl_R.setZero();
        ctl_R.diagonal() << 0.01, 0.01;
        Q_mpc.setZero();
        Q_mpc.diagonal() << 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
        Qf_mpc.setZero();
        Qf_mpc.diagonal() << 0, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4;
        Cost running_cost_horizon(target_state, Q_mpc, ctl_R);
        TerminalCost terminal_cost_horizon(target_state, Qf_mpc);
        
        OptimizerResult<Dynamics> results_horizon;
        results_horizon.control_trajectory = hor_control;
        
        results_horizon = ddp_horizon.run_horizon(cur_state, hor_control, hor_traj_states, *ddp_dyn, running_cost_horizon, terminal_cost_horizon);
        u = results_horizon.control_trajectory.col(0);


        mpc_writer.save_step(cur_state, u);

      }

      double tau_L = 0, tau_R = 0;
      if(mpc_steps > -1) {

        // *************************************** IDEA 2
        double ddth = u(0);
        double tau_0 = u(1);
        State xdot = ddp_dyn->f(cur_state, u);
        double ddx = xdot(3);
        double ddpsi = xdot(4);
        Eigen::Vector3d ddq, dq;
        ddq << ddx, ddpsi, ddth;
        dq = cur_state.segment(3,3);
        c_forces dy_forces = ddp_dyn->dynamic_forces(cur_state, u);
        //double tau_1 = (dy_forces.A.block<1,3>(2,0)*ddq) + (dy_forces.C.block<1,3>(2,0)*dq) + (dy_forces.Q(2)) - (dy_forces.Gamma_fric(2));
        double tau_1 = dy_forces.A.block<1,3>(2,0)*ddq;
        tau_1 += dy_forces.C.block<1,3>(2,0)*dq;
        tau_1 += dy_forces.Q(2);
        tau_1 -= dy_forces.Gamma_fric(2);
        tau_L = -0.5*(tau_1+tau_0);
        tau_R = -0.5*(tau_1-tau_0);

        if(abs(tau_L) > 60 | abs(tau_R) > 60){
          cout << "step: " << steps << ", tau_0: " << tau_0 << ", tau_1: " << tau_1 << ", tau_L: " << tau_L << ", tau_R: " << tau_R << endl;
        }
      }
      mForces(0) = tau_L;
      mForces(1) = tau_R;
      const vector<size_t > index{6, 7};
      mkrang->setForces(index, mForces);

      
      SimWindow::timeStepping();
    }


    ~MyWindow() {
      outFile.close();     
    }
    

  protected:

    SkeletonPtr m3DOF;
    SkeletonPtr mkrang;

    Eigen::VectorXd qInit;

    Eigen::VectorXd dof1;

    double psi, dpsi, qBody1, dqBody1, dthL, dthR;
    double psiFilt, dpsiFilt, qBody1Filt, dqBody1Filt, dthLFilt, dthRFilt;

    double R;
    double L;
    
    int steps;
    int mpc_steps;
    double mpc_dt;

    Eigen::Matrix<double, 2, 1> mForces;   
   
    ofstream outFile; 

    filter *dqFilt, *cFilt;
    ControlTrajectory ddp_ctl_traj;
    StateTrajectory ddp_state_traj;
    Dynamics *ddp_dyn;
    Control u;
    CSV_writer<Scalar> mpc_writer;

    WorldPtr world3dof;


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




struct comOptParams {
  SkeletonPtr robot;
  Eigen::Matrix<double, 25, 1> qInit;
};

double comOptFunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  comOptParams* optParams = reinterpret_cast<comOptParams *>(my_func_data);
  Eigen::Matrix<double, 25, 1> q(x.data());

  if (!grad.empty()) {
    Eigen::Matrix<double, 25, 1> mGrad = q-optParams->qInit;
    Eigen::VectorXd::Map(&grad[0], mGrad.size()) = mGrad;
  }
  return (0.5*pow((q-optParams->qInit).norm(), 2));
}

double comConstraint(const std::vector<double> &x, std::vector<double> &grad, void *com_const_data) {
  comOptParams* optParams = reinterpret_cast<comOptParams *>(com_const_data);
  Eigen::Matrix<double, 25, 1> q(x.data());
  optParams->robot->setPositions(q);
  return (pow(optParams->robot->getCOM()(0)-optParams->robot->getPosition(3), 2) \
    + pow(optParams->robot->getCOM()(1)-optParams->robot->getPosition(4), 2));
}

double wheelAxisConstraint(const std::vector<double> &x, std::vector<double> &grad, void *wheelAxis_const_data) {
  comOptParams* optParams = reinterpret_cast<comOptParams *>(wheelAxis_const_data);
  Eigen::Matrix<double, 25, 1> q(x.data());
  optParams->robot->setPositions(q);
  return optParams->robot->getBodyNode(0)->getTransform().matrix()(2,0);
}

double headingConstraint(const std::vector<double> &x, std::vector<double> &grad, void *heading_const_data) {
  comOptParams* optParams = reinterpret_cast<comOptParams *>(heading_const_data);
  Eigen::Matrix<double, 25, 1> q(x.data());
  optParams->robot->setPositions(q);
  Eigen::Matrix<double, 4, 4> Tf = optParams->robot->getBodyNode(0)->getTransform().matrix();
  double heading = atan2(Tf(0,0), -Tf(1,0));
  optParams->robot->setPositions(optParams->qInit);
  Tf = optParams->robot->getBodyNode(0)->getTransform().matrix();
  double headingInit = atan2(Tf(0,0), -Tf(1,0));
  return heading-headingInit;
}



dart::dynamics::SkeletonPtr createKrang() {
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  dart::dynamics::SkeletonPtr krang =
      loader.parseSkeleton("/home/krang/dart/09-URDF/Krang/KrangOld.urdf");
  krang->setName("krang");

  // Read initial pose from the file
  ifstream file("/home/krang/dart/Unification/10d-3DHighLevelControllerC-2Ddart_Unification/examples/3dofddp_original/defaultInit.txt");
  assert(file.is_open());
  char line [1024];
  file.getline(line, 1024);
  std::istringstream stream(line);
  Eigen::Matrix<double, 24, 1> initPoseParams; // heading, qBase, x, y, z, qLWheel, qRWheel, qWaist, qTorso, qKinect, qLArm0, ... qLArm6, qRArm0, ..., qRArm6
  size_t i = 0; double newDouble;
  while((i < 24) && (stream >> newDouble)) initPoseParams(i++) = newDouble;
  file.close();
  double headingInit; headingInit = initPoseParams(0);
  double qBaseInit; qBaseInit = initPoseParams(1);
  Eigen::Vector3d xyzInit; xyzInit << initPoseParams.segment(2,3);
  double qLWheelInit; qLWheelInit = initPoseParams(5);
  double qRWheelInit; qRWheelInit = initPoseParams(6);
  double qWaistInit; qWaistInit = initPoseParams(7);
  double qTorsoInit; qTorsoInit = initPoseParams(8);
  double qKinectInit; qKinectInit = initPoseParams(9);
  Eigen::Matrix<double, 7, 1> qLeftArmInit; qLeftArmInit << initPoseParams.segment(10, 7);
  Eigen::Matrix<double, 7, 1> qRightArmInit; qRightArmInit << initPoseParams.segment(17, 7);
  
  // Calculating the axis angle representation of orientation from headingInit and qBaseInit: 
  // RotX(pi/2)*RotY(-pi/2+headingInit)*RotX(-qBaseInit)
  Eigen::Transform<double, 3, Eigen::Affine> baseTf = Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  baseTf.prerotate(Eigen::AngleAxisd(-qBaseInit,Eigen::Vector3d::UnitX())).prerotate(Eigen::AngleAxisd(-M_PI/2+headingInit,Eigen::Vector3d::UnitY())).prerotate(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()));
  Eigen::AngleAxisd aa(baseTf.matrix().block<3,3>(0,0));

  // Ensure CoM is right on top of wheel axis
  const int dof = (const int)krang->getNumDofs();
  comOptParams optParams;
  optParams.robot = krang;
  optParams.qInit << aa.angle()*aa.axis(), xyzInit, qLWheelInit, qRWheelInit, qWaistInit, qTorsoInit, qKinectInit, qLeftArmInit, qRightArmInit; 
  nlopt::opt opt(nlopt::LN_COBYLA, dof);
  std::vector<double> q_vec(dof);
  double minf;
  opt.set_min_objective(comOptFunc, &optParams);
  opt.add_equality_constraint(comConstraint, &optParams, 1e-8);
  opt.add_equality_constraint(wheelAxisConstraint, &optParams, 1e-8);
  opt.add_equality_constraint(headingConstraint, &optParams, 1e-8);
  opt.set_xtol_rel(1e-4);
  opt.set_maxtime(10);
  opt.optimize(q_vec, minf);
  Eigen::Matrix<double, 25, 1> q(q_vec.data());
  
  // Initializing the configuration
  // file = ifstream("/home/krang/dart/Unification/10d-3DHighLevelControllerC-2Ddart_Unification/examples/3dofddp_original/overwriteInitPose.txt");
  // assert(file.is_open());
  // file.getline(line, 1024);
  // stream = std::istringstream(line);
  // i = 0;
  // while((i < 25) && (stream >> newDouble)) q(i++) = newDouble;
  // file.close();
  q(3) = 0.0; q(4) = 0.0; q(6) = 0.0; q(7) = 0.0; 
  krang->setPositions(q);
  int joints = krang->getNumJoints();
  for(int i=3; i < joints; i++) {
    krang->getJoint(i)->setActuatorType(dart::dynamics::Joint::ActuatorType::LOCKED);
  }

  return krang;
}


dart::dynamics::SkeletonPtr createTray(dart::dynamics::BodyNodePtr ee) {

  Eigen::Matrix3d EELRot, trayLocalRot, trayRot;
  Eigen::Vector3d EELPos, trayLocalTranslation, trayPos;
  Eigen::Matrix<double, 6, 1> qObject;
  
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  dart::dynamics::SkeletonPtr tray =
      loader.parseSkeleton("/home/krang/dart/09-URDF/scenes/tray.urdf");
  tray->setName("tray");

  // Orientation
  EELRot = ee->getTransform().rotation();
  trayLocalRot << 1,  0, 0,
                  0,  0, 1,
                  0, -1, 0;
  trayRot = EELRot*trayLocalRot;
  Eigen::AngleAxisd aa(trayRot);
  
  // Position
  EELPos = ee->getTransform().translation();
  trayLocalTranslation << 0, -0.286303, 0.04;
  trayPos = EELPos + trayRot*trayLocalTranslation;

  // Set the position
  qObject << (aa.angle()*aa.axis()), trayPos;
  tray->setPositions(qObject);

  return tray;
}

dart::dynamics::SkeletonPtr createCup(dart::dynamics::BodyNodePtr ee) {
  
  Eigen::Matrix3d EELRot, cupLocalRot, cupRot;
  Eigen::Vector3d EELPos, cupLocalTranslation, cupPos;
  Eigen::Matrix<double, 6, 1> qObject;
  
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  dart::dynamics::SkeletonPtr cup =
      loader.parseSkeleton("/home/krang/dart/09-URDF/scenes/cup.urdf");
  cup->setName("cup");

  // Orientation
  EELRot = ee->getTransform().rotation();
  cupLocalRot << 1,  0, 0,
                 0,  0, 1,
                 0, -1, 0;
  cupRot = EELRot*cupLocalRot;
  Eigen::AngleAxisd aa(cupRot);
  
  // Position
  EELPos = ee->getTransform().translation();
  cupLocalTranslation << 0, -0.286303, 0.06;
  cupPos = EELPos + cupRot*cupLocalTranslation;

  // Set the position
  qObject << (aa.angle()*aa.axis()), cupPos;
  cup->setPositions(qObject);

  return cup;
}


int main(int argc, char* argv[])
{

  // create and initialize the world
  cout << "1" << endl;
 // SkeletonPtr threeDOF = create3DOF_URDF();
  cout << "2" << endl;

    // load skeletons
  SkeletonPtr floor = createFloor();
  SkeletonPtr robot = createKrang();
  cout << "3" << endl;

  SkeletonPtr tray = createTray(robot->getBodyNode("lGripper"));
  cout << "4" << endl;

  SkeletonPtr cup = createCup(robot->getBodyNode("lGripper")); //cup->setPositions(tray->getPositions());
  cout << "5" << endl;

  WorldPtr world = std::make_shared<World>();

  world->addSkeleton(floor); //add ground and robot to the world pointer
  world->addSkeleton(robot);
  world->addSkeleton(tray);
  world->addSkeleton(cup);
  cout << "6" << endl;

  // create and initialize the world
  //Eigen::Vector3d gravity(0.0,  -9.81, 0.0);
  //world->setGravity(gravity);
  // world->setTimeStep(1.0/1000);
  cout << "7" << endl; 

  MyWindow window(world);

  glutInit(&argc, argv);
  window.initWindow(1280,720, "3DOF URDF");
  glutMainLoop();

  return 0;
}
