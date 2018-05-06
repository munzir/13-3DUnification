#ifndef QUADROTOR_COSTS_H
#define QUADROTOR_COSTS_H

#include <ddp/costs.hpp>
#include <Eigen/Dense>

template <typename DynamicsT>
struct QuadrotorCost: public CostFunction<DynamicsT>
{
    using Dynamics  = DynamicsT;
    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Control   = typename Dynamics::Control;
    using Gradient  = typename CostFunction<Dynamics>::Gradient;
    using Hessian   = typename CostFunction<Dynamics>::Hessian;

    static const int N = Dynamics::StateSize;
    static const int M = Dynamics::ControlSize;
    using StateCostWeight = Eigen::Matrix<Scalar, N, N>;
    using ControlCostWeight = Eigen::Matrix<Scalar, M, M>;

public:
    QuadrotorCost(const Eigen::Ref<const StateCostWeight> &Q, const Eigen::Ref<const ControlCostWeight> &R)
    : Q_(Q), R_(R)
    {
        QR_.setZero();
        QR_.topLeftCorner(Dynamics::StateSize, Dynamics::StateSize) = Q;
        QR_.bottomRightCorner(Dynamics::ControlSize, Dynamics::ControlSize) = R;
    }

    Scalar c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        State error = x - this->target();
        return (error.transpose() * Q_ * error).value() + (u.transpose() * R_ * u).value();
    }

    Gradient dc(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        return (Gradient() << Q_ * (x - this->target()), R_ * u).finished();
    }

    Hessian d2c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        return QR_;
    }

private:
    StateCostWeight Q_;
    ControlCostWeight R_;
    Hessian QR_;
};

template <typename DynamicsT>
struct QuadrotorTerminalCost: public TerminalCostFunction<DynamicsT>
{
    using Dynamics  = DynamicsT;
    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Gradient  = typename TerminalCostFunction<Dynamics>::Gradient;
    using Hessian   = typename TerminalCostFunction<Dynamics>::Hessian;

    static const int N = Dynamics::StateSize;
    using StateCostWeight = Eigen::Matrix<Scalar, N, N>;

public:
    QuadrotorTerminalCost(const Eigen::Ref<const StateCostWeight> &Qf)
    : Qf_(Qf) {}

    Scalar c(const Eigen::Ref<const State> &x)
    {
        return static_cast<Scalar>(0.5) * ((x - this->target()).transpose() * Qf_ * (x - this->target())).value();
    }

    Gradient dc(const Eigen::Ref<const State> &x)
    {
        return Qf_ * (x - this->target());
    }

    Hessian d2c(const Eigen::Ref<const State> &x)
    {
        return Qf_;
    }

private:
    StateCostWeight Qf_;
};

#endif // QUADROTOR_COSTS_H
