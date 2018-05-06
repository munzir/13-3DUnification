#ifndef TRAJOPT_AUTORALLY_COSTS_H
#define TRAJOPT_AUTORALLY_COSTS_H

#include <ddp/costs.hpp>
#include <Eigen/Dense>

template <class Dynamics>
struct AutoRallyCost: public CostFunction<Dynamics>
{
    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Control   = typename Dynamics::Control;
    using Gradient  = typename CostFunction<Dynamics>::Gradient;
    using Hessian   = typename CostFunction<Dynamics>::Hessian;

    AutoRallyCost(const State &xf,
                  const Eigen::Ref<Eigen::Matrix<Scalar, 16, 1>> &track_coef,
                  Scalar speed_coef, Scalar crash_coef, Scalar stability_cost, Scalar max_slip_angle, Scalar steering_coef, Scalar throttle_coef, Scalar track_weight)
            : CostFunction<Dynamics>(xf),
              tc_(track_coef), ms_(max_slip_angle), spc_(speed_coef), cc_(crash_coef), sc_(stability_cost), steering_cost_(steering_coef), throttle_cost_(throttle_coef), track_weight_(track_weight)
    {
        grad_.setZero();
        hess_.setZero();
    }

    inline void update_speed(const Scalar desired_speed)
    {
        if(desired_speed < static_cast<Scalar>(0.0))
        {
            this->target(4) = static_cast<Scalar>(0.0);
        }
        else
        {
            this->target(4) = desired_speed;
        }
    }

    inline bool crashed(const Eigen::Ref<const State> &state) const
    {
        bool crash = false;
        if(std::abs(state(3)) > 1.57) // Don't roll too much
        {
            crash = true;
        }
        else if(std::abs(state(6)) > 5.0) // Don't yaw too rapidly
        {
            crash = true;
        }
        else if((state(4) > 0.01) && (std::abs(-std::atan(state(5)/std::abs(state(4)))) > 0.785))
        {
            // If the slip angle is above 45 degrees ...
            crash = true;
        }
        return crash;
    }

    inline bool unstable(const Eigen::Ref<const State> &state) const
    {
        bool high_slip = false;
        if(state(4) > static_cast<Scalar>(1.0))
        {
            Scalar slip = -std::atan(state(5) / std::abs(state(4)));
            if(std::abs(slip) > ms_)
            {
                high_slip = true;
            }
        }
        return high_slip;
    }

    inline Scalar c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {   
        return (crashed(x) ? cc_ * std::pow(static_cast<Scalar>(1.0) + x(4), 2) : static_cast<Scalar>(0.0))
                        + (unstable(x) ? sc_ : static_cast<Scalar>(0.0))
                        + spc_*std::pow((x(4) - this->target(4)), 2)
                        + steering_cost_*std::pow(u(0), 2) + throttle_cost_*std::pow(u(1), 2)
                        + track_weight_ *
                            ( std::pow((tc_(0) + tc_(1)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(4)*x(0) + tc_(5)*x(0)*x(1) + tc_(6)*x(0)*std::pow(x(1), 2) +
                              tc_(7)*x(0)*std::pow(x(1), 3) + tc_(8)*std::pow(x(0), 2) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) +
                              tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(12)*std::pow(x(0), 3) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) +
                              tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3)), 2)
                            );
    }

    inline Gradient dc(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        grad_(0) = 2*(tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) + tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) +
                 tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        grad_(1) = 2*(tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) + 3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) + tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) +
                 tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        grad_(0) *= track_weight_;
        grad_(1) *= track_weight_;
        grad_(4) = 2 * spc_* (x(4) - this->target(4));
        grad_(7) = 2 * steering_cost_ * u(0);
        grad_(8) = 2 * throttle_cost_ * u(1);

        return grad_;
    }

    inline Hessian d2c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        hess_(0, 0) = 2*std::pow((tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3)), 2) +
                2*(2*tc_(8) + 6*tc_(12)*x(0) + 2*tc_(9)*x(1) + 6*tc_(13)*x(0)*x(1) + 2*tc_(10)*std::pow(x(1), 2) + 6*tc_(14)*x(0)*std::pow(x(1), 2) + 2*tc_(11)*std::pow(x(1), 3) + 6*tc_(15)*x(0)*std::pow(x(1), 3))*(tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) +
                  tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) + tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(0, 1) = 2*(tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) + 3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2))*
                (tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3)) +
               2*(tc_(5) + 2*tc_(9)*x(0) + 3*tc_(13)*std::pow(x(0), 2) + 2*tc_(6)*x(1) + 4*tc_(10)*x(0)*x(1) + 6*tc_(14)*std::pow(x(0), 2)*x(1) + 3*tc_(7)*std::pow(x(1), 2) + 6*tc_(11)*x(0)*std::pow(x(1), 2) + 9*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 2))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) +
                 tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(1, 0) = 2*(tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) + 3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2))*
                (tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3)) +
               2*(tc_(5) + 2*tc_(9)*x(0) + 3*tc_(13)*std::pow(x(0), 2) + 2*tc_(6)*x(1) + 4*tc_(10)*x(0)*x(1) + 6*tc_(14)*std::pow(x(0), 2)*x(1) + 3*tc_(7)*std::pow(x(1), 2) + 6*tc_(11)*x(0)*std::pow(x(1), 2) + 9*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 2))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) +
                 tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(1, 1) = 2*std::pow((tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) +
                        3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2)), 2) + 2*(2*tc_(2) + 2*tc_(6)*x(0) + 2*tc_(10)*std::pow(x(0), 2) + 2*tc_(14)*std::pow(x(0), 3) + 6*tc_(3)*x(1) + 6*tc_(7)*x(0)*x(1) + 6*tc_(11)*std::pow(x(0), 2)*x(1) + 6*tc_(15)*std::pow(x(0), 3)*x(1))*
                      (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) +
                       tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(0, 0) *= track_weight_;
        hess_(0, 1) *= track_weight_;
        hess_(1, 0) *= track_weight_;
        hess_(1, 1) *= track_weight_;
        hess_(4, 4) = 2 * spc_;
        hess_(7, 7) = 2 * steering_cost_;
        hess_(8, 8) = 2 * throttle_cost_;

        return hess_;
    }

    Eigen::Matrix<Scalar, 16, 1> tc_;
    Scalar spc_, cc_, sc_, ms_, steering_cost_, throttle_cost_, track_weight_;
    Gradient grad_;
    Hessian hess_;
};

template <typename Dynamics>
struct AutoRallyTerminalCost: public TerminalCostFunction<Dynamics>
{
    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Control   = typename Dynamics::Control;
    using Gradient  = typename TerminalCostFunction<Dynamics>::Gradient;
    using Hessian   = typename TerminalCostFunction<Dynamics>::Hessian;

    AutoRallyTerminalCost(const State &xf, const Eigen::Ref<Eigen::Matrix<Scalar, 16, 1>> &track_coef, Scalar speed_coef, Scalar crash_coef, Scalar stability_cost, Scalar max_slip_angle, Scalar track_weight)
            : TerminalCostFunction<Dynamics>(xf),
              tc_(track_coef), ms_(max_slip_angle), spc_(speed_coef), cc_(crash_coef), sc_(stability_cost), track_weight_(track_weight)
    {
        grad_.setZero();
        hess_.setZero();
    }

    inline void update_speed(const Scalar desired_speed)
    {
        if(desired_speed < static_cast<Scalar>(0.0))
        {
            this->target(4) = static_cast<Scalar>(0.0);
        }
        else
        {
            this->target(4) = desired_speed;
        }
    }

    inline bool crashed(const Eigen::Ref<const State> &state) const
    {
        bool crash = false;
        if(std::abs(state(3)) > 1.57) // Don't roll too much
        {
            crash = true;
        }
        else if(std::abs(state(6)) > 5.0) // Don't yaw too rapidly
        {
            crash = true;
        }
        else if((state(4) > 0.01) && (std::abs(-std::atan(state(5)/std::abs(state(4)))) > 0.785))
        {
            // If the slip angle is above 45 degrees ...
            crash = true;
        }
        return crash;
    }

    inline bool unstable(const Eigen::Ref<const State> &state) const
    {
        bool high_slip = false;
        if(state(4) > static_cast<Scalar>(1.0))
        {
            Scalar slip = -std::atan(state(5) / std::abs(state(4)));
            if(std::abs(slip) > ms_)
            {
                high_slip = true;
            }
        }
        return high_slip;
    }

    inline Scalar c(const Eigen::Ref<const State> &x)
    {
        return (crashed(x) ? cc_ * std::pow(static_cast<Scalar>(1.0) + x(4), 2) : static_cast<Scalar>(0.0))
                + (unstable(x) ? sc_ : static_cast<Scalar>(0.0))
                + spc_*std::pow((x(4) - this->target(4)), 2)
                + track_weight_ *
                    ( std::pow((tc_(0) + tc_(1)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(4)*x(0) + tc_(5)*x(0)*x(1) + tc_(6)*x(0)*std::pow(x(1), 2) +
                      tc_(7)*x(0)*std::pow(x(1), 3) + tc_(8)*std::pow(x(0), 2) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) +
                      tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(12)*std::pow(x(0), 3) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) +
                      tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3)), 2)
                    );
    }

    inline Gradient dc(const Eigen::Ref<const State> &x)
    {
        grad_(0) = 2*(tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) + tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) +
                 tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        grad_(1) = 2*(tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) + 3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) + tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) +
                 tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        grad_(4) = 2*spc_*(x(4) - this->target(4));
        grad_(0) *= track_weight_;
        grad_(1) *= track_weight_;
        return grad_;
    }

    inline Hessian d2c(const Eigen::Ref<const State> &x)
    {
        hess_(0, 0) = 2*std::pow((tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3)), 2) +
                2*(2*tc_(8) + 6*tc_(12)*x(0) + 2*tc_(9)*x(1) + 6*tc_(13)*x(0)*x(1) + 2*tc_(10)*std::pow(x(1), 2) + 6*tc_(14)*x(0)*std::pow(x(1), 2) + 2*tc_(11)*std::pow(x(1), 3) + 6*tc_(15)*x(0)*std::pow(x(1), 3))*(tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) +
                  tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) + tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(0, 1) = 2*(tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) + 3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2))*
                (tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3)) +
               2*(tc_(5) + 2*tc_(9)*x(0) + 3*tc_(13)*std::pow(x(0), 2) + 2*tc_(6)*x(1) + 4*tc_(10)*x(0)*x(1) + 6*tc_(14)*std::pow(x(0), 2)*x(1) + 3*tc_(7)*std::pow(x(1), 2) + 6*tc_(11)*x(0)*std::pow(x(1), 2) + 9*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 2))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) +
                 tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(1, 0) = 2*(tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) + 3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2))*
                (tc_(4) + 2*tc_(8)*x(0) + 3*tc_(12)*std::pow(x(0), 2) + tc_(5)*x(1) + 2*tc_(9)*x(0)*x(1) + 3*tc_(13)*std::pow(x(0), 2)*x(1) + tc_(6)*std::pow(x(1), 2) + 2*tc_(10)*x(0)*std::pow(x(1), 2) + 3*tc_(14)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(7)*std::pow(x(1), 3) + 2*tc_(11)*x(0)*std::pow(x(1), 3) + 3*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 3)) +
               2*(tc_(5) + 2*tc_(9)*x(0) + 3*tc_(13)*std::pow(x(0), 2) + 2*tc_(6)*x(1) + 4*tc_(10)*x(0)*x(1) + 6*tc_(14)*std::pow(x(0), 2)*x(1) + 3*tc_(7)*std::pow(x(1), 2) + 6*tc_(11)*x(0)*std::pow(x(1), 2) + 9*tc_(15)*std::pow(x(0), 2)*std::pow(x(1), 2))*
                (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) +
                 tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(1, 1) = 2*std::pow((tc_(1) + tc_(5)*x(0) + tc_(9)*std::pow(x(0), 2) + tc_(13)*std::pow(x(0), 3) + 2*tc_(2)*x(1) + 2*tc_(6)*x(0)*x(1) + 2*tc_(10)*std::pow(x(0), 2)*x(1) + 2*tc_(14)*std::pow(x(0), 3)*x(1) + 3*tc_(3)*std::pow(x(1), 2) + 3*tc_(7)*x(0)*std::pow(x(1), 2) + 3*tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 2) +
                        3*tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 2)), 2) + 2*(2*tc_(2) + 2*tc_(6)*x(0) + 2*tc_(10)*std::pow(x(0), 2) + 2*tc_(14)*std::pow(x(0), 3) + 6*tc_(3)*x(1) + 6*tc_(7)*x(0)*x(1) + 6*tc_(11)*std::pow(x(0), 2)*x(1) + 6*tc_(15)*std::pow(x(0), 3)*x(1))*
                      (tc_(0) + tc_(4)*x(0) + tc_(8)*std::pow(x(0), 2) + tc_(12)*std::pow(x(0), 3) + tc_(1)*x(1) + tc_(5)*x(0)*x(1) + tc_(9)*std::pow(x(0), 2)*x(1) + tc_(13)*std::pow(x(0), 3)*x(1) + tc_(2)*std::pow(x(1), 2) + tc_(6)*x(0)*std::pow(x(1), 2) + tc_(10)*std::pow(x(0), 2)*std::pow(x(1), 2) + tc_(14)*std::pow(x(0), 3)*std::pow(x(1), 2) + tc_(3)*std::pow(x(1), 3) + tc_(7)*x(0)*std::pow(x(1), 3) +
                       tc_(11)*std::pow(x(0), 2)*std::pow(x(1), 3) + tc_(15)*std::pow(x(0), 3)*std::pow(x(1), 3));
        hess_(4, 4) = 2 * spc_;
        hess_(0, 0) *= track_weight_;
        hess_(0, 1) *= track_weight_;
        hess_(1, 0) *= track_weight_;
        hess_(1, 1) *= track_weight_;
        return hess_;
    }

    Eigen::Matrix<Scalar, 16, 1> tc_;
    Scalar spc_, cc_, sc_, ms_, track_weight_;
    Gradient grad_;
    Hessian hess_;
};

#endif //TRAJOPT_AUTORALLY_COSTS_H
