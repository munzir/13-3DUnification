#include "krang3d.h"
#include "krang3dwindow.h"
#include <ddp/costs.hpp>
#include <ddp/ddp.hpp>
#include <ddp/mpc.hpp>
#include <ddp/util.hpp>
#include <QApplication>
#include <QtConcurrent/QtConcurrent>
#include <iostream>
#include <vector>
#include <fstream>

int main(int argc, char **argv)
{
    using Scalar = double;
    using Dynamics = Krang3D<Scalar>;
    using DDP_Opt = optimizer::DDP<Dynamics>;
    using Plant = Krang3DPlant;
    using Cost = Krang3DCost<Scalar>;
    using TerminalCost = Krang3DTerminalCost<Scalar>;
    using StateTrajectory = typename Dynamics::StateTrajectory ;
    using ControlTrajectory= typename Dynamics::ControlTrajectory ;

    // Log output
    util::DefaultLogger logger;

    CSV_writer<Scalar> writer;
    bool verbose = true;

    // Timing
    Scalar tf = 20;
    Scalar dt = 0.01;
    auto time_steps = util::time_steps(tf, dt);
    int max_iterations = 15;
//    std::cout << "1" << std::endl;
    // Parameters for the dynamics
    param p;
    p.R = 2.500000e-01; p.mw = 5.100000e-01; p.Iw = 5.100000e-03; p.L = 6.000000e-01; p.g=9.800000e+00;
    p.m_1 = 1.342160e+02;
    p.MX_1 = 0; p.MY_1 = 8.052960e+01; p.MZ_1 = 0;
    p.XX_1 = 6.487107e+01; p.YY_1 = 4.473867e+00; p.ZZ_1 = 6.845016e+01;
    p.XY_1 = 0; p.YZ_1 = 0; p.XZ_1 = 0;
    p.fric_1 = 15;
    p.XXw = 1.673438e-02; p.YYw=3.283125e-02; p.ZZw=1.673438e-02;
//    std::cout << "2" << std::endl;

    // Dynamics
    Dynamics cp_dynamics(p);
//    std::cout << "3" << std::endl;

    // Noisy Plant
    Scalar ssigma = 0.00;
    Scalar csigma = 0.00;
    Plant cp_plant(cp_dynamics, dt, ssigma, csigma);
//    std::cout << "4" << std::endl;

    // Initial state th, dth, x, dx, desired state, initial control sequence
    Dynamics::State x0 = Dynamics::State::Zero();
    Dynamics::State xf; xf << 2, 0, 0, 0, 0, 0, 0.01, 5;
    Dynamics::ControlTrajectory u = Dynamics::ControlTrajectory::Zero(2, time_steps);
//    std::cout << "5" << std::endl;

    // Costs
    Cost::StateHessian Q;
    Q.setZero();
    Q.diagonal() << 0,0.1,0.1,0.1,0.1,0.1,0.1,0.1;
//    std::cout << "6" << std::endl;

    Cost::ControlHessian R;
    R.setZero();
    R.diagonal() << 0.01, 0.01;
//    std::cout << "7" << std::endl;

    TerminalCost::Hessian Qf;
    Qf.setZero();
    Qf.diagonal() << 0,1e4,1e4,1e4,1e4,1e4,1e4,1e4;
//    std::cout << "8" << std::endl;

    Cost cp_cost(xf, Q, R);
    TerminalCost cp_terminal_cost(xf, Qf);
//    std::cout << "9" << std::endl;

    // initialize DDP for trajectory planning
    DDP_Opt trej_ddp (dt, time_steps, max_iterations, &logger, verbose);
//    std::cout << "10" << std::endl;

    // Get initial trajectory from DDP
    OptimizerResult<Dynamics> traj_results = trej_ddp.run(x0, u, cp_dynamics, cp_cost, cp_terminal_cost);
//    std::cout << "11" << std::endl;

    StateTrajectory xs = traj_results.state_trajectory;
    ControlTrajectory us = traj_results.control_trajectory;

    writer.save_trajectory(xs, us, "initial_traj.csv");

    logger.info("Obtained initial state trajectory");
    for (int m = 0; m < xs.cols(); ++m) {
        logger.info("\n");
        for (int n = 0; n < xs.rows(); ++n) {
            logger.info("%f ", xs(n, m));
        }
    }

    logger.info("Obtained initial control trajectory");
    for (int m = 0; m < us.cols(); ++m) {
        logger.info("\n");
        for (int n = 0; n < us.rows(); ++n) {
            logger.info("%f ", us(n, m));
        }
    }
    // Connect signals and configure plot window
    QApplication a(argc, argv);
    // final angle, then position
    MainWindow win(xf(6), xf(7), xf(2), dt);
    QObject::connect(&cp_plant, &Plant::update, &win, &MainWindow::update, Qt::QueuedConnection);
    win.show();

    // Create new MPC
    using Result = ModelPredictiveController<Dynamics, optimizer::DDP>::Result;
    using StateRef = Eigen::Ref<const Dynamics::State>;

    // Define a termination condition (when no more threads alive)
    auto termination =
            [&](int i, Result &r, const StateRef &x, Cost &rc, TerminalCost &tc, Scalar true_cost)
            {
                return (QThreadPool::globalInstance()->activeThreadCount() == 0);
            };


    // ORIGINAL MPC

    // Initialize receding horizon controller
    // ModelPredictiveController<Dynamics, optimizer::DDP> mpc(dt, time_steps, max_iterations, &logger, verbose, /*DDP options*/ &logger, verbose);

    // Run in receding horizon mode
    // QtConcurrent::run([&](){ mpc.run(x0, u, termination, cp_dynamics, cp_plant, cp_cost, cp_terminal_cost); });

    // END ORIGINAL MPC

    // Run in receding horizon mode
    QtConcurrent::run([&](){

        int horizon = 10;
        Dynamics::State x = x0;
        Dynamics::State x_old = x0;

        StateTrajectory mpc_xs = xs.replicate(1, 1);
        ControlTrajectory mpc_us = us.replicate(1, 1);

        mpc_xs.col(0) = x0;
        // For Each Time Step
        for (int k = 0; k < (time_steps - 1); ++k) {
            int look_ahead;
            if ((time_steps - k) > horizon) {
               look_ahead = horizon;
            } else {
                look_ahead = (time_steps - 1) - k;
            }

            logger.info("\n");
            logger.info("iteration: %d", k);
            logger.info("look ahead: %d", look_ahead);

            // New Goal State
            Dynamics::State x_target = xs.col(k + look_ahead);

            Dynamics::ControlTrajectory u_horizon = Dynamics::ControlTrajectory::Zero(2, look_ahead);
            Dynamics::StateTrajectory x_horizon_ref = xs.block(0, k, 8, look_ahead);

//            logger.info("Obtained horizon state trajectory from ddp: ");
//            for(int m = 0; m < x_horizon_ref.cols(); ++m) {
//                logger.info("\n");
//                for (int n = 0; n < x_horizon_ref.rows(); ++n) {
//                    logger.info("%f ", x_horizon_ref(n, m));
//                }
//            }
//
//            logger.info("Obtained horizon control trajectory from ddp: ");
//            for(int m = 0; m <u_horizon.cols(); ++m) {
//                logger.info("\n");
//                for (int n = 0; n < u_horizon.rows(); ++n) {
//                    logger.info("%f ", u_horizon(n, m));
//                }
//            }
            // initialize DDP for trajectory planning
            DDP_Opt ddp_horizon (dt, look_ahead, max_iterations, &logger, verbose);

            Cost::StateHessian Q_mpc, Qf_mpc;
            Q_mpc.setZero();
            Q_mpc.diagonal() << 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
            Qf_mpc.setZero();
            Qf_mpc.diagonal() << 0, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4;
            Cost cp_cost_horizon(x_target, Q_mpc, R);
            TerminalCost cp_terminal_cost_horizon(x_target, Qf_mpc);

//
//            logger.info("\n");
//            logger.info("current state");
//            for (int n=0; n < 4; ++n) {
//                logger.info("\n");
//                logger.info("%f", x(n));
//                logger.info("\n");
//            }

//            logger.info("\n");
//            logger.info("new target");
//            for (int m =0; m < 4; ++m) {
//                logger.info("\n");
//                logger.info("%f", x_target(m));
//                logger.info("\n");
//            }

            OptimizerResult<Dynamics> results_horizon;
            results_horizon.control_trajectory = u_horizon;
            // Get initial trajectory from DDP

            results_horizon = ddp_horizon.run_horizon(x_old, u_horizon, x_horizon_ref, cp_dynamics, cp_cost_horizon, cp_terminal_cost_horizon);
//            results_horizon = ddp_horizon.run(x_old, u_horizon, cp_dynamics, cp_cost_horizon, cp_terminal_cost_horizon);

            Dynamics::Control u_new = results_horizon.control_trajectory.col(0);
            Dynamics::StateTrajectory x_horizon = results_horizon.state_trajectory;

//            logger.info("Obtained horizon state trajectory from optimizer: ");
//            for(int m = 0; m < xs_horizon.cols(); ++m) {
//                logger.info("\n");
//                for (int n = 0; n < xs_horizon.rows(); ++n) {
//                    logger.info("%f ", xs_horizon(n, m));
//                }
//            }
//            logger.info("\n");
//
//            logger.info("\n Obtained control trajectory from optimizer: \n");
//            for(int k =0; k < u_new.cols(); ++k) {
//                logger.info("%f", u_new(k));
//                logger.info("\n");
//            }
//
//
//            logger.info("Obtained control from optimizer: ");
//            for(int m = 0; m < u_new.rows(); ++m) { logger.info("%f ", u_new(m)); }
//            logger.info("\n");

            logger.info("new control %f", u_new(0));
            logger.info("reference state %f", xs(2, k));
            x = cp_plant.f_with_ref(x_old, u_new, xs.col(k), us.col(k));
            x_old = x;
            mpc_us.col(k) = u_new;
            mpc_xs.col(k + 1) = x;
        }

        writer.save_trajectory(mpc_xs, mpc_us, "mpc_trajectory.csv");

    });

    return a.exec();
}

void save_traj() {

}
