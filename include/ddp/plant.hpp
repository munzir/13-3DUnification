#ifndef TRAJOPT_PLANT_HPP
#define TRAJOPT_PLANT_HPP

#include <Eigen/Dense>

template <class T, int S, int C>
struct Plant
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    enum { StateSize = S, ControlSize = C };
    using Scalar            = T;
    using State             = Eigen::Matrix<T, StateSize, 1>;
    using Control           = Eigen::Matrix<T, ControlSize, 1>;
    using StateTrajectory   = Eigen::Matrix<T, StateSize, Eigen::Dynamic>;
    using ControlTrajectory = Eigen::Matrix<T, ControlSize, Eigen::Dynamic>;

    Plant() = default;
    Plant(const Plant &other) = default;
    Plant(Plant &&other) = default;
    Plant& operator=(const Plant &other) = default;
    Plant& operator=(Plant &&other) = default;
    virtual ~Plant() = default;

    /**
     * @brief   Pass information to the Plant in order to apply controls and obtain the new state.
     *
     * The user must provide an implementation that takes in the system state and the control calculated
     * by the optimizer, applies one control to the system, performs any desired calculations or updates,
     * and returns the new state of the system after the application of the control.
     *
     * Note that unlike the Dynamics functor, this function must return a STATE rather than a state transition.
     * This is in the interests of compatibility with real robotic systems whose state estimators will usually
     * return a state rather than a transition.
     *
     * @param x The state calculated by the optimizer for the current time window.
     * @param u The control calculated by the optimizer for the current time window.
     * @return  The new state of the system.
     */
    virtual State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) = 0;
};

#endif // TRAJOPT_PLANT_HPP
