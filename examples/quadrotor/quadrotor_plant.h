#ifndef QUADROTOR_PLANT_H
#define QUADROTOR_PLANT_H

#include "quadrotor.h"
#include "quadrotor_costs.h"
#include <ddp/plant.hpp>
#include <ddp/eigenmvn.hpp>
#include <Eigen/Dense>
#include <QtCore>
#include <QVector>

struct QuadrotorPlant: public QObject, public Plant<typename Quadrotor<double>::Scalar, Quadrotor<double>::StateSize, Quadrotor<double>::ControlSize>
{
    using Dynamics  = Quadrotor<double>;
    using Cost      = QuadrotorCost<Dynamics>;
    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Control   = typename Dynamics::Control;

    static const int N = Dynamics::StateSize;
    static const int M = Dynamics::ControlSize;
    using StateNoiseVariance    = Eigen::Matrix<Scalar, N, N>;
    using ControlNoiseVariance  = Eigen::Matrix<Scalar, M, M>;

    QuadrotorPlant(Dynamics &dynamics, Cost &cost, Scalar dt,
                   const Eigen::Ref<const StateNoiseVariance> &sn,
                   const Eigen::Ref<const ControlNoiseVariance> &cn)
    : dynamics_(dynamics), cost_(cost), dt_(dt),
      sn_(Eigen::Matrix<Scalar, N, 1>::Zero(), sn),
      cn_(Eigen::Matrix<Scalar, M, 1>::Zero(), cn)
    {

    }

    State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        static QVector<Scalar> sx(StateSize);
        static QVector<Scalar> su(ControlSize);
        static QVector<Scalar> sxf(StateSize);

        Control u_noisy = u + cn_.samples(1);
        State xnew = x + (dynamics_.f(x, u_noisy) + sn_.samples(1)) * dt_;
        Eigen::Map<State>(sx.data(), StateSize) = xnew;
        Eigen::Map<Control>(su.data(), ControlSize) = u_noisy;
        Eigen::Map<State>(sxf.data(), StateSize) = cost_.target();
        update(sx, su, sxf, cost_.c(xnew, u_noisy));
        return xnew;
    }

signals:
    void update(const QVector<Scalar> &x, const QVector<Scalar> &u, const QVector<Scalar> &xf, Scalar true_cost);

private:
    Q_OBJECT

    Dynamics &dynamics_;
    Cost &cost_;
    Scalar dt_;
    Eigen::EigenMultivariateNormal<Scalar, N> sn_;
    Eigen::EigenMultivariateNormal<Scalar, M> cn_;
};

#endif // QUADROTOR_PLANT_H
