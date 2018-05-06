#ifndef INVPENDULUM_H
#define INVPENDULUM_H

#include <ddp/costs.hpp>
#include <ddp/dynamics.hpp>

template <class T>
struct InvPendulum: public Dynamics<T, 2, 1>
{
    using Scalar = T;
    using State = typename Dynamics<Scalar, 2, 1>::State;
    using Control = typename Dynamics<Scalar, 2, 1>::Control;

    InvPendulum(Scalar m, Scalar l, Scalar b, Scalar g)
    : m_(m), l_(l), b_(b), g_(g) {}

    State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) const
    {
        State dx;
        dx(0) = x(1);
        dx(1) = ((-b_ / (m_ * l_)) * x(1) - (g_ / l_) * sin(x(0))) + u(0);
    }

private:
    Scalar m_, l_, b_, g_;
};

template <class T>
struct InvPendulumCost: public CostFunction<InvPendulum<T>>
{
    using Scalar = T;
    using Dynamics = InvPendulum<T>;
    using State = typename CostFunction<InvPendulum<T>>::State;
    using Control = typename CostFunction<InvPendulum<T>>::Control;
    using Gradient = typename CostFunction<InvPendulum<T>>::Gradient;
    using Hessian = typename CostFunction<InvPendulum<T>>::Hessian;
    using StateHessian = Eigen::Matrix<Scalar, InvPendulum<T>::StateSize, InvPendulum<T>::StateSize>;
    using ControlHessian = Eigen::Matrix<Scalar, InvPendulum<T>::ControlSize, InvPendulum<T>::ControlSize>;

    InvPendulumCost(const Eigen::Ref<const State> &xf, const Eigen::Ref<const StateHessian> &Q, const Eigen::Ref<const ControlHessian> &R)
    : CostFunction<InvPendulum<T>>(xf), Q_(Q), R_(R)
    {
        QR_.setZero();
        QR_.topLeftCorner(Dynamics::StateSize, Dynamics::StateSize) = Q;
        QR_.bottomRightCorner(Dynamics::ControlSize, Dynamics::ControlSize) = R;
    }

    Scalar c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) const
    {
        return ((x - CostFunction<InvPendulum<T>>::xf).transpose() * Q_ * (x - CostFunction<InvPendulum<T>>::xf) +
                u.transpose() * R_ * u).value();
    }

    Gradient dc(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) const
    {
        Gradient g;
        g.head(Dynamics::StateSize) = Q_ * (x - CostFunction<InvPendulum<T>>::xf);
        g.tail(Dynamics::ControlSize) = R_ * u;
        return g;
    }

    Hessian d2c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) const
    {
        return QR_;
    }

private:
    StateHessian Q_;
    ControlHessian R_;
    Hessian QR_;
};

template <class T>
struct InvPendulumTerminalCost: public TerminalCostFunction<InvPendulum<T>>
{
    using Scalar = T;
    using Dynamics = InvPendulum<T>;
    using State = typename TerminalCostFunction<InvPendulum<T>>::State;
    using Gradient = typename TerminalCostFunction<InvPendulum<T>>::Gradient;
    using Hessian = typename TerminalCostFunction<InvPendulum<T>>::Hessian;

    InvPendulumTerminalCost(const Eigen::Ref<const State> &xf, const Eigen::Ref<const Hessian> &Q)
    : TerminalCostFunction<InvPendulum<T>>(xf), Q_(Q) {}

    Scalar c(const Eigen::Ref<const State> &x) const
    {
        return (x - TerminalCostFunction<InvPendulum<T>>::xf).transpose() *
                Q_ * (x - TerminalCostFunction<InvPendulum<T>>::xf);
    }

    Gradient dc(const Eigen::Ref<const State> &x) const
    {
        return Q_ * (x - TerminalCostFunction<InvPendulum<T>>::xf);
    }

    Hessian d2c(const Eigen::Ref<const State> &x) const
    {
        return Q_;
    }

private:
    Hessian Q_;
};

#endif // INVPENDULUM_H
