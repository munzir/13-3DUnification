#include "quadrotor.h"
#include "quadrotor_costs.h"
#include "quadrotor_plant.h"
#include "quadrotorwindow.h"
#include <ddp/ddp.hpp>
#include <ddp/mpc.hpp>
#include <ddp/result.hpp>
#include <ddp/util.hpp>
#include <QApplication>
#include <QtConcurrent/QtConcurrent>

int main(int argc, char **argv)
{
    // Convenience typedefs
    using Scalar        = double;
    using Dynamics      = Quadrotor<Scalar>;
    using RunningCost   = QuadrotorCost<Dynamics>;
    using TerminalCost  = QuadrotorTerminalCost<Dynamics>;
    using Plant         = QuadrotorPlant;

    // Logger
    util::DefaultLogger logger;
    bool verbose = true;

    // Dynamics
    Scalar m = 1.0;
    Scalar l = 0.24;
    Scalar Jx = 8.1e-3, Jy = 8.1e-3, Jz = 14.2e-3;
    Scalar g = 9.81;
    Dynamics dynamics(m, l, Jx, Jy, Jz, g);

    // Receding horizon
    Scalar tf = 1.0;
    Scalar dt = 0.01;
    auto time_steps = util::time_steps(tf, dt);

    // Total time steps
    Scalar tN = 7.5;
    auto total_time_steps = util::time_steps(tN, dt);

    // Iterations per time step
    int iterations = 2;

    // Initial state
    Dynamics::State x0;
    x0.setZero();
    x0(0) = 1.0;
    x0(1) = 0.0;
    x0(2) = 1.0;

    // Trajectory to track (same size as total time steps)
    Dynamics::StateTrajectory xf;
    xf.setZero(Dynamics::StateSize, total_time_steps);
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> steps = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>::LinSpaced(total_time_steps, 0.0, 2.0 * M_PI);
    xf.row(0) = steps.array().cos();
    xf.row(1) = steps.array().sin();
    xf.row(2) = 1.0 * Eigen::Matrix<Scalar, 1, Eigen::Dynamic>::Ones(total_time_steps);

    // Initial control trajectory (same size as receding horizon)
    Dynamics::ControlTrajectory u = Dynamics::ControlTrajectory::Zero(Dynamics::ControlSize, time_steps);

    // Running state, running control, and terminal state cost
    RunningCost::StateCostWeight Q;
    Q.setIdentity();
    Q(0,0) = 250.0;
    Q(1,1) = 250.0;
    Q(2,2) = 250.0;

    RunningCost::ControlCostWeight R;
    R.setIdentity();
    R *= 0.015;

    TerminalCost::StateCostWeight Qf;
    Qf.setIdentity();
    Qf(0,0) = 1000.0;
    Qf(1,1) = 1000.0;
    Qf(2,2) = 1000.0;

    RunningCost run_cost(Q, R);
    TerminalCost term_cost(Qf);

    // Simulate noisy plant
    Scalar sv = 0.10;
    Scalar cv = 0.20;
    auto state_var = sv * Plant::StateNoiseVariance::Identity();
    auto control_var = cv * Plant::ControlNoiseVariance::Identity();
    Plant plant(dynamics, run_cost, dt, state_var, control_var);

    // Initialize MPC
    ModelPredictiveController<Dynamics, optimizer::DDP> rhc(dt, time_steps, iterations, &logger, verbose,
                                                            /*DDP options*/ &logger, verbose);

    // Configure and display plot window
    QApplication a(argc, argv);
    MainWindow win(dt);
    QObject::connect(&plant, &Plant::update, &win, &MainWindow::update, Qt::QueuedConnection);
    win.show();

    // Define a termination condition
    // Update the target point at each time step
    using Result = ModelPredictiveController<Dynamics, optimizer::DDP>::Result;
    using StateRef = Eigen::Ref<const Dynamics::State>;
    auto termination =
        [&](int i, Result &r, const StateRef &x, RunningCost &rc, TerminalCost &tc, Scalar true_cost)
        {
            rc.target() = xf.col(i % total_time_steps);
            tc.target() = xf.col(i % total_time_steps);
            return (QThreadPool::globalInstance()->activeThreadCount() == 0);
        };

    // Run in receding horizon mode
    QtConcurrent::run([&](){ rhc.run(x0, u, termination, dynamics, plant, run_cost, term_cost); });

    return a.exec();
}
