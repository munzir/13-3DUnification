#ifndef TRAJOPT_QUADROTOR_H
#define TRAJOPT_QUADROTOR_H

#include <ddp/costs.hpp>
#include <ddp/dynamics.hpp>
#include <ddp/eigenmvn.hpp>
#include <cmath>

template <typename T>
Eigen::Matrix<T, 4, 1> forces(const Eigen::Ref<const Eigen::Matrix<T, 4, 1>> &u)
{
    static constexpr T constant = static_cast<T>(0.5 * std::sqrt(2.0));

    return (Eigen::Matrix<T, 4, 1>() <<
                u.sum(),
                constant * (u(0) + u(2) - u(1) - u(3)),
                constant * (u(2) + u(3) - u(0) - u(1)),
                u(0) + u(3) - u(1) - u(2)
            ).finished();
}

template <typename T>
Eigen::Matrix<T, 3, 3> inertial_to_body(T roll, T pitch, T yaw)
{
    using std::sin; using std::cos;

    return (Eigen::Matrix<T, 3, 3>() <<
                cos(pitch) * cos(yaw),                                      cos(pitch) * sin(yaw),                                      -sin(pitch),
                sin(pitch) * sin(roll) * cos(yaw) - cos(roll) * sin(yaw),   sin(pitch) * sin(roll) * sin(yaw) + cos(roll) * cos(yaw),   sin(roll) * cos(pitch),
                sin(roll) * sin(yaw) + cos(roll) * sin(pitch) * cos(yaw),   sin(pitch) * sin(yaw) * cos(roll) - sin(roll) * sin(yaw),   cos(pitch) * cos(roll)
            ).finished().transpose();
}

template <typename T>
Eigen::Matrix<T, 3, 3> body_to_euler(T roll, T pitch, T yaw)
{
    using std::sin; using std::cos; using std::tan;

    return (Eigen::Matrix<T, 3, 3>() <<
                static_cast<T>(1.0),    tan(pitch) * sin(roll),   tan(pitch) * cos(roll),
                static_cast<T>(0.0),    cos(roll),                -sin(roll),
                static_cast<T>(0.0),    sin(roll) / cos(pitch),   cos(roll) / cos(pitch)
            ).finished();
}

template <typename T>
struct Quadrotor: public Dynamics<T, 12, 4>
{
    using Scalar    = T;
    using State     = typename Dynamics<Scalar, 12, 4>::State;
    using Control   = typename Dynamics<Scalar, 12, 4>::Control;

    Quadrotor(Scalar mass = 0.5, Scalar arm_length = 0.17,
              Scalar Jx = 3.2e-3, Scalar Jy = 3.2e-3, Scalar Jz = 5.5e-3,
              Scalar kt = 1.69e-2, Scalar g = 9.81)
    : minv(1.0 / mass), L(arm_length), Jx(Jx), Jy(Jy), Jz(Jz), kt(kt), g(g) {}

    inline State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        auto force = forces(u);
        Scalar torque_roll = L * force(1);
        Scalar torque_pitch = L * force(2);
        Scalar torque_yaw = L * force(3);

        auto rot_ib = inertial_to_body(x(6), x(7), x(8));
        auto er = body_to_euler(x(6), x(7), x(8));

        State dx;
        dx(0) = x(3);
        dx(1) = x(4);
        dx(2) = x(5);

        dx(3) = minv * rot_ib(0, 2) * force(0);
        dx(4) = minv * rot_ib(1, 2) * force(0);
        dx(5) = minv * rot_ib(2, 2) * force(0) - g;

        dx(6) = er(0, 0) * x(9) + er(0, 1) * x(10) + er(0, 2) * x(11);
        dx(7) = er(1, 1) * x(10) + er(1, 2) * x(11);
        dx(8) = er(2, 1) * x(10) + er(2, 2) * x(11);

        dx(9) = (1.0 / Jx) * ((Jy - Jz) * x(10) * x(11) + torque_roll);
        dx(10) = (1.0 / Jy) * ((Jz - Jx) * x(9) * x(11) + torque_pitch);
        dx(11) = (1.0 / Jz) * ((Jx - Jy) * x(9) * x(10) + torque_yaw);

        return dx;
    }

    Scalar minv;
    Scalar L;
    Scalar Jx, Jy, Jz;
    Scalar kt;
    Scalar g;
};

template <typename T>
struct QuadrotorSimple: public Dynamics<T, 12, 4>
{
    using Scalar    = typename Dynamics<T, 12, 4>::Scalar;
    using State     = typename Dynamics<T, 12, 4>::State;
    using Control   = typename Dynamics<T, 12, 4>::Control;

    QuadrotorSimple(Scalar m, Scalar l, Scalar Jx, Scalar Jy, Scalar Jz, Scalar g)
    : minv_(1.0 / m), l_(l), Ixx_(Jx), Iyy_(Jy), Izz_(Jz),
      Ixxinv_(static_cast<T>(1.0) / Ixx_), Iyyinv_(static_cast<T>(1.0) / Iyy_), Izzinv_(static_cast<T>(1.0) / Izz_),
      g_(g) {}

    inline State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        using std::sin; using std::cos;
        State state_dot;

        Scalar u1 = u(0) + u(1) + u(2) + u(3);
        Scalar u2 = u(3) - u(1);
        Scalar u3 = u(0) - u(2);
        Scalar u4 = 0.05 * (u(1) + u(3) - u(0) - u(2));

        state_dot(0) = x(3);
        state_dot(1) = x(4);
        state_dot(2) = x(5);

        state_dot(3) = minv_ * (cos(x(6)) * sin(x(7)) * cos(x(8)) + sin(x(6)) * sin(x(8))) * u1;
        state_dot(4) = minv_ * (cos(x(6)) * sin(x(7)) * sin(x(8)) - sin(x(6)) * cos(x(8))) * u1;
        state_dot(5) = -g_ + minv_ * (cos(x(6)) * cos(x(7))) * u1;

        state_dot(6) = x(9);
        state_dot(7) = x(10);
        state_dot(8) = x(11);

        state_dot(9)  = Ixxinv_ * (x(10)*x(11) * (Iyy_ - Izz_) + l_ * u2);
        state_dot(10) = Iyyinv_ * (x(9)*x(11) * (Izz_ - Ixx_) + l_ * u3);
        state_dot(11) = Izzinv_ * (x(9)*x(10) * (Ixx_ - Iyy_) + u4);

        return state_dot;
    }

    T minv_;
    T l_;
    T Ixx_;
    T Iyy_;
    T Izz_;
    T Ixxinv_;
    T Iyyinv_;
    T Izzinv_;
    T g_;
};

#endif // TRAJOPT_QUADROTOR_H
