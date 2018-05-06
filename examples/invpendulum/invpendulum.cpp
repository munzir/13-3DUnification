#include "invpendulum.h"
#include "invpendulumwindow.h"
#include <ddp/ddp.hpp>
#include <QApplication>

int main(int argc, char **argv)
{
    using Dynamics      = InvPendulum<double>;
    using Cost          = InvPendulumCost<double>;
    using TerminalCost  = InvPendulumTerminalCost<double>;

    double tf           = 2.0;
    double dt           = 0.01;
    auto time_steps     = util::time_steps(tf, dt);
    int max_iterations  = 10;

    Dynamics::State x0 = Dynamics::State::Zero();
    Dynamics::State xf; xf << M_PI, 0.0;
    Dynamics::ControlTrajectory u = Dynamics::ControlTrajectory::Zero(time_steps);

    Cost::StateHessian Q;
    Q.setZero(); /*Q.diagonal() << 1.0, 1.0;*/
    Cost::ControlHessian R;
    R << 0.1;
    TerminalCost::Hessian Qf;
    Qf.setZero(); Qf.diagonal() << 100.0, 0.0;

    Dynamics ip_dynamics(1.0, 0.25, 0.01, 9.81);
    Cost ip_cost(xf, Q, R);
    TerminalCost ip_terminal_cost(xf, Qf);

    DDP<Dynamics> solver(dt, time_steps, max_iterations);
    auto result = solver.run(x0, u, ip_dynamics, ip_cost, ip_terminal_cost);

    // Plot targets
    std::vector<double> target_angle(result.state_trajectory.cols(), xf(0));
    std::vector<double> target_velocity(result.state_trajectory.cols(), xf(1));

    // Plot trajectory over time
    std::vector<double> angle(result.state_trajectory.cols());
    std::vector<double> velocity(result.state_trajectory.cols());
    std::vector<double> time(result.state_trajectory.cols(), 0.0);
    for(int t = 0; t < result.timesteps; ++t)
    {
        angle.at(t) = result.state_trajectory(0, t);
        velocity.at(t) = result.state_trajectory(1, t);
        if(t > 0) time.at(t) = time.at(t - 1) + dt;
    }

    // Plot cost per iteration
    std::vector<double> cost(result.iterations);
    Eigen::VectorXd cpi = result.cost.rowwise().sum();
    Eigen::VectorXd::Map(&cost[0], result.iterations) = cpi;

    // Setup and show the plot window
    QApplication a(argc, argv);
    InvPendulumWindow win(angle, target_angle, velocity, target_velocity, time, cost);
    win.show();
    return a.exec();
}
