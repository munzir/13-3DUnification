#ifndef AUTORALLYRT_H
#define AUTORALLYRT_H

#include "autorally_bfm.h"
#include <ddp/plant.hpp>
#include <ddp/eigenmvn.hpp>
#include <QtCore>
#include <QVector>

struct AutoRallyPlotter: public QObject, public Plant<double, 7, 2>
{
    using Base                  = Plant<double, 7, 2>;
    using Scalar                = typename Base::Scalar;
    using State                 = typename Base::State;
    using Control               = typename Base::Control;
    using StateNoiseVariance    = Eigen::Matrix<Scalar, StateSize, StateSize>;
    using ControlNoiseVariance  = Eigen::Matrix<Scalar, ControlSize, ControlSize>;

    AutoRallyPlotter(AutoRally<Scalar> &dynamics, const StateNoiseVariance &state_noise, const ControlNoiseVariance &control_noise, Scalar dt,
                     QObject *parent = nullptr)
    : QObject(parent),
      dynamics_(dynamics),
      state_noise_dist_(State::Zero(), state_noise), control_noise_dist_(Control::Zero(), control_noise),
      dt_(dt)
    {
    }

    inline State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        static QVector<Scalar> sx(StateSize);
        static QVector<Scalar> su(ControlSize);

        Control u_noisy = u + control_noise_dist_.samples(1);
        saturate(u_noisy);

        State xnew = x + (dynamics_.f(x, u_noisy) + state_noise_dist_.samples(1)) * dt_;
        Eigen::Map<State>(sx.data(), StateSize) = xnew;
        Eigen::Map<Control>(su.data(), ControlSize) = u_noisy;
        signal_update(sx, su);
        return xnew;
    }

signals:
    void signal_update(const QVector<Scalar> &state, const QVector<Scalar> &control) const;

private:
    Q_OBJECT
    AutoRally<Scalar> dynamics_;
    Eigen::EigenMultivariateNormal<Scalar, StateSize> state_noise_dist_;
    Eigen::EigenMultivariateNormal<Scalar, ControlSize> control_noise_dist_;
    Scalar dt_;

    void saturate(Eigen::Ref<Control> u) const
    {
        for(int i = 0; i < Plant<double, 7, 2>::ControlSize; ++i)
        {
            if(u(i) < -1.0)
                u(i) = -1.0;
            if(u(i) > 1.0)
                u(i) = 1.0;
        }
    }
};

#endif // AUTORALLYRT_H
