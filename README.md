# TrajOpt
Real-time trajectory optimization.

## Setup (not required for use; required only to build this project with examples)

* Install [CMake](https://cmake.org/download/)
* Install [Qt](https://www.qt.io/download-open-source/)
* (optional) Install Intel MKL
* Clone this repository

## Use

Simply drop the contents of the `include/` directory into your project.
Note that you will also need to include [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) if it is not accessible system-wide.
A copy is included with this project in the `eigen/` directory.

To use this software for your own systems:
* Implement an object representing your system's dynamics by deriving from the `Dynamics` class and implementing `f` (for the system dynamics)
and `df` (for the Jacobian of dynamics with respect to both the state and controls).  `df` will be determined using numerical differentiation
if no implementation is provided.

* Implement an object representing a plant to control by deriving from the `Plant` class and implementing `f`.
No derivative is required.  Note that the `Dynamics` object representing system dynamics may be used as a `Plant` as well with some minor changes.
See the `Autorally` and `Quadrotor` examples.

* Create objects representing the running and terminal costs by deriving from the `CostFunction` and `TerminalCostFunction` classes,
respectively, and implementing `c` (for cost), `dc` (for the cost gradient with respect to both the state and the control),
and `d2c` (for the second derivative of the cost with respect to both the state and the control).
Note that the cost function derivatives ***will not*** be determined numerically in the current implementation -- derivatives must be provided,
though the user's provided implementation is free to use numerical derivatives.

Finally, run the optimizer by providing the dynamics, costs, and other required parameters (optimizer, time step, time horizon, etc.).
 See the example applications which demonstrate how to use the library in vanilla or receding horizon modes.

## Examples
A few example applications are included (cart pole, rally car, quadrotor).  These are built by default when **cmake** is used to build the project.
 **Qt5** is required to build the examples; see the **Setup** section above.

    git clone https://github.gatech.edu/ACDS/DDP.git
    cd DDP
    mkdir build
    cd build
    cmake -DCMAKE_PREFIX_PATH=/path/to/qt5 ..
    make (for GCC or Clang users on Linux or Mac -- Windows users should use the appropriate command depending on their compiler of choice)

The example applications will be in the `build/examples` folder after building.

### Cart Pole ###

The cart pole example demonstrates how to use the library in vanilla mode, i.e., single-shot optimization of a trajectory with a provided optimizer.
The Differential Dynamic Programming (DDP) optimizer is the currently the only optimizer provided.

### Autorally ###

The Autorally example demonstrates how to use the library in receding horizon mode to continuously attempt to minimize the cost function.
No trajectory is provided in this case; the path to follow is encoded in the cost function.

### Quadrotor ###

The quadrotor example demonstrates how to use the library to track a user-provided trajectory.
