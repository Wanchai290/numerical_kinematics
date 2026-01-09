# Numerical Kinematics
## How to run ?
The file `newton_2d.m` is an exemple usage of the Newton-Raphson method for
linear optimization.

The project in itself is contained in the `main.m` file. There are in total 3
sections of the code described below

### First section : Direct geometric symbolic model
Defines the model, as well as the target point `P_target` for the robot to attain. The condition for each joint to be at `2 * Obs_r` away from the obstacle is written at line 83.

|Line|Parameter|Description|
|-----|----------|----------------|
| 132 |    Lv    | Length of each joint
| 132 |    q0    | Initial solution for Newton-Raphson solver
| 144 | P_target | Column vector of the desired position for the end-effector
----------------------------------

### Second section : Inverse geometric model solver
Uses the Newton-Raphson method to solve the inverse geometric model of the robot, along with obstacle avoidance (if enabled by passing a parameter).

|Line|Parameter|Description|
|-----|-------|----------------|
| 184 | Obs_P | Column vector of the obstacle's center
| 185 | Obs_r | Obstacle radius
--------------------------------

### Third section : Resulting 3D plot
Displays a 3D plot with both the simple solution and the solution that avoids obstacles.

## About pseudo-inverse computation
We used a result of the SVD to perform it, but really we should've changed the totality of the code to use the SVD decomposition instead, so computations run faster. But there was no time left for this project.
