#ifndef TRAJOPT_AUTORALLY_BFM_H
#define TRAJOPT_AUTORALLY_BFM_H

#include <ddp/dynamics.hpp>
#include <ddp/util.hpp>
#include <fstream>
#include <iostream>
#include <string>

template <class T>
struct AutoRally: public Dynamics<T, 7, 2>
{
    static const int NUM_BFS = 25;
    using Scalar                    = typename Dynamics<T, 7, 2>::Scalar;
    using State                     = typename Dynamics<T, 7, 2>::State;
    using Control                   = typename Dynamics<T, 7, 2>::Control;
    using StateTrajectory           = typename Dynamics<T, 7, 2>::StateTrajectory;
    using ControlTrajectory         = typename Dynamics<T, 7, 2>::ControlTrajectory;
    using FeedbackGainTrajectory    = typename Dynamics<T, 7, 2>::FeedbackGainTrajectory;
    using FeedforwardGainTrajectory = typename Dynamics<T, 7, 2>::FeedforwardGainTrajectory;

    static const int AUTORALLY_BAD_FILE_PATH = -2;

    AutoRally(const std::string &filename, util::Logger *logger)
    {
        std::ifstream theta_file(filename);
        if(!theta_file)
        {
            logger->error("Unable to open theta file");
            std::exit(AUTORALLY_BAD_FILE_PATH);
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < NUM_BFS; j++)
                {
                    theta_file >> theta_(i, j);
                }
            }
        }

        theta_file.close();
    }

    inline State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        State dx;
        static Eigen::Matrix<T, NUM_BFS, 1> bf_vec;

        // Compute basis functions
        for (int i = 0; i < NUM_BFS; i++)
        {
            bf_vec(i) = basisFuncX(i, x.data(), u.data());
        }

        // Kinematics
        dx(0) = std::cos(x(2)) * x(4) - std::sin(x(2)) * x(5);
        dx(1) = std::sin(x(2)) * x(4) + std::cos(x(2)) * x(5);
        dx(2) = -x(6);

        // Dynamics
        dx.tail(4) = theta_ * bf_vec;

        return dx;
    }

private:
    Eigen::Matrix<T, 4, NUM_BFS, Eigen::RowMajor> theta_;

    T basisFuncX(int idx, const T *state, const T *control) const
    {
      T phi = 0;
      switch(idx) {
        case 0: phi = control[1];
                break;
        case 1: phi = state[4]/10.0;
                break;
        case 2: phi = (state[4] > .1) ?
                      std::sin(control[0])*std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0])/1200.0 :
                      std::sin(control[0])*std::tan(-control[0])/1200.0;
                break;
        case 3: phi = (state[4] > .1) ?
                      std::sin(control[0])*std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0])*std::abs(std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0]))/1440000.0:
                      std::sin(control[0])*std::tan(-control[0])*std::abs(std::tan(-control[0]))/1440000.0;
                break;
        case 4: phi = (state[4] > .1) ?
                      std::sin(control[0])*std::pow(std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0]), 3)/1728000000.0 :
                      std::sin(control[0])*std::pow(std::tan(-control[0]), 3)/1728000000.0;
                break;
        case 5: phi = state[6]*state[5]/25.0;
                break;
        case 6: phi = state[6]/10.0;
                break;
        case 7: phi = state[5]/10.0;
                break;
        case 8: phi = std::sin(control[0]);
                break;
        case 9: phi = (state[4] > .1) ?
                       state[5]/state[4]/40.0 :
                       0;
                break;
        case 10: phi = (state[4] > .1) ?
                        std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0])/1400.0 :
                        std::tan(-control[0])/1400.0;
                break;
        case 11: phi = (state[4] > .1) ?
                        std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0])*std::abs(std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0]))/1960000:
                        std::tan(-control[0])*std::abs(std::tan(-control[0]))/1960000;
                break;
        case 12: phi = (state[4] > .1) ?
                        std::pow(std::tan(std::atan(state[5]/state[4] + .45*state[6]/state[4]) - control[0]), 3)/2744000000:
                        std::pow(std::tan(-control[0]), 3)/2744000000;
                break;
        case 13: phi = (state[4] > .1) ?
                       (state[5]/state[4] - .35*state[6]/state[4])/40.0 :
                        0;
                break;
        case 14: phi = (state[4] > .1) ?
                       (state[5]/state[4] - .35*state[6]/state[4])*std::abs(state[5]/state[4] - .35*state[6]/state[4])/1600.0:
                        0;
                break;
        case 15: phi = (state[4] > .1) ?
                        std::pow(state[5]/state[4] - .35*state[6]/state[4], 3)/64000.0 :
                        0;
                break;
        case 16: phi = state[6]*state[4]/50.0;
                break;
        case 17: phi = state[3];
                break;
        case 18: phi = state[3]*state[6];
                break;
        case 19: phi = state[3]*state[4]/3.0;
                break;
        case 20: phi = state[3]*state[4]*state[6]/5.0;
                break;
        case 21: phi = std::pow(state[4], 2)/100.0;
                break;
        case 22: phi = std::pow(state[4], 3)/1000.0;
                break;
        case 23: phi = std::pow(control[1], 2);
                break;
        case 24: phi = std::pow(control[1], 3);
                break;
      }
      return phi;
    }
};

#endif // TRAJOPT_AUTORALLY_BFM_H
