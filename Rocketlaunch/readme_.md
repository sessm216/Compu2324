
Simulation of the trajectory of a rocket traveling within the gravitational system defined by the Earth and the Moon. 4 coupled first-order ODE. 2D orbital mechanics. 
Using the 4th order Runge-Kutta (RK) numerical integration method: solves non-linear eq of motion in polar coordinates. Variables conjugate momenta are used to calculate H=V+K to verify numerical stability. 

Suppositions:
-Earth's referential frame is fixed, does not move. Phi and theta angles wrt the Eart.
-Energy needed for launch is considered an instant pulse. 
-Angular v of Moon is constant. 

About the eq of motion: Polar coordinates $(r, \phi)$, corresponding canonical momenta $(p_r, p_\phi)$

Changing the launch point: changing phi0 and theta0. 
theta0: launch angle from a specific point. Compute as: theta0=wishedangle*(Ì/180.0)
phi0: starting angle relative to the reference axis. PI/2 is north pole. 


**State Vector $f$**:
$f[0] = r$: Radial distance from Earth's center (normalized by Earth-Moon distance $d_{TL}$).
$f[1] = \phi$: Polar angle of the rocket.
$f[2] = p_r$: Radial momentum component.
$f[3] = p_\phi$: Angular momentum component.
$f[4] = t$: Simulation time (seconds).

### System Constants

| **$G$** | `G` | $6.67 \times 10^{-11}$ | $\text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$ | Gravitational Constant |
| **$M_t$** | `Mt` | $5.9736 \times 10^{24}$ | $\text{kg}$ | Mass of Earth |
| **$M_l$** | `Ml` | $0.07349 \times 10^{24}$ | $\text{kg}$ | Mass of Moon |
| **$R_T$** | `RT` | $6.37816 \times 10^6$ | $\text{m}$ | Mean radius of Earth |
| **$d_{TL}$** | `dTl` | $3.844 \times 10^8$ | $\text{m}$ | Earth-Moon distance |
| **$\omega$** | `w` | $2.6617 \times 10^{-6}$ | $\text{rad/s}$ | Moon's orbital angular velocity |


3 output files:
`naveangulo_6.dat`: Contains Cartesian coordinates $(x, y)$ of both the rocket and the Moon. Downsampled to record every 60 steps to optimize output file size. Necessary for video visualization.
`Hnormalizado_6.dat: Stores normalized energy ratio $H / H_0$. Useful for verifying energy conservation (ideal value $\approx 1.0$).
'Hangulo_6'.dat: stores the hamiltonian (total energy) at every iteration step. 

Names can be updated. 
> ⚠️ **Note**: Check the hardcoded output directories in `main()` before compiling (default path: `/home/compu2324/cohete/`). Ensure the target directory exists or update the paths as needed.

About compilation: 

Requirements:standard C++ compiler (e.g., `g++` or `clang++`).

Build: compile the source code using standard optimization flags:

```bash
g++ -O2 simulation.cpp -o version3.cpp
```

Final notes: check directory and Hnormalized computation. 

