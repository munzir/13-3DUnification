#include "autorallyrt.h"
#include "autorallyrtwindow.h"
#include "autorally_bfm.h"
#include "autorally_costs.h"
#include <ddp/ddp.hpp>
#include <ddp/mpc.hpp>
#include <ddp/util.hpp>
#include <QApplication>
#include <QtConcurrent/QtConcurrent>

int main(int argc, char **argv)
{
    using Dynamics = AutoRally<double>;
    using Scalar = Dynamics::Scalar;
    using Plant = AutoRallyPlotter;
    using RunningCost = AutoRallyCost<Dynamics>;
    using TerminalCost = AutoRallyTerminalCost<Dynamics>;

    // Create a logger for output
    util::DefaultLogger logger;
    bool verbose = true;

    // Dynamics
    Dynamics car_dynamics("THETA.txt", &logger);

    // Noise
    Scalar sn_multiplier = 0.05;
    Scalar cn_multiplier = 0.05;
    Plant::StateNoiseVariance state_noise = sn_multiplier * Plant::StateNoiseVariance::Identity();
    Plant::ControlNoiseVariance control_noise = cn_multiplier * Plant::ControlNoiseVariance::Identity();

    // Plant
    Scalar dt = 1.0 / 100.0;
    Plant car_plant(car_dynamics, state_noise, control_noise, dt);

    // Initial state and desired states
    Dynamics::State x0;
    x0.setZero();
    x0(2) = static_cast<Scalar>(2.35);

    Dynamics::State xf;
    xf.setZero();
    xf(4) = static_cast<Scalar>(11.0);

    // Costs
    Eigen::Matrix<Scalar, 16, 1> simulated_track;
    simulated_track << -1.31621629e-01,  4.93592932e-01,  4.22996447e-02, -1.82646619e-03,
                        5.92302211e-01,  7.99102512e-03,  1.46501357e-03,  3.50704914e-06,
                        3.37724237e-02, -8.54907425e-03,  5.88551838e-05,  3.54671893e-05,
                       -6.84942760e-05, -3.26088526e-04, -1.40257190e-06,  1.25879429e-06;
    Scalar  speed_coef =       1.75,    // speed cost weight
            crash_coef =       10000.0, // crash cost weight
            stability_cost =   100.0,   // (in)stability cost weight
            max_slip =         0.75,    // max allowable slip angle
            steering_coef =    5.0,    // steering cost weight
            throttle_coef =    5.0,    // throttle cost weight
            track_weight =     20.0;    // track coefficient MULTIPLIER
    RunningCost car_cost(xf, simulated_track, speed_coef, crash_coef, stability_cost, max_slip, steering_coef, throttle_coef, track_weight);
    TerminalCost car_term_cost(xf, simulated_track, speed_coef, crash_coef, stability_cost, max_slip, track_weight);

    // Time and iterations
    Scalar tf = 1.50;
    int time_steps = util::time_steps(tf, dt);
    int iterations = 10;

    // Initial control trajectory
    Dynamics::ControlTrajectory u;
    u.setZero(Dynamics::ControlSize, time_steps);

    // Control limits
    Dynamics::Control u_min = Dynamics::Control::Constant(static_cast<Scalar>(-1.0));
    Dynamics::Control u_max = Dynamics::Control::Constant(static_cast<Scalar>(1.0));

    // Initialize receding horizon controller
    ModelPredictiveController<Dynamics, optimizer::DDP> rhc(dt, time_steps, iterations, &logger, verbose, /*DDP options*/ &logger, verbose);

    // Connect signal and configure plot window
    QApplication a(argc, argv);
    AutoRallyRTWindow win(dt);
    QObject::connect(&car_plant, &Plant::signal_update, &win, &AutoRallyRTWindow::update_graph, Qt::QueuedConnection);

    // Define a termination condition
    using Result = ModelPredictiveController<Dynamics, optimizer::DDP>::Result;
    using StateRef = Eigen::Ref<const Dynamics::State>;
    auto termination =
        [=](int i, Result &r, const StateRef &x, RunningCost &rc, TerminalCost &tc, Scalar true_cost)
        {
            return (QThreadPool::globalInstance()->activeThreadCount() == 0);
        };

    // Run in receding horizon mode until the user quits the program
    QtConcurrent::run([&](){ rhc.run(x0, u, termination, car_dynamics, car_plant, car_cost, car_term_cost, /*DDP options*/ u_min, u_max); });

    // Show plot window
    win.show();

    return a.exec();
}
