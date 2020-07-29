#ifndef GRADIENT_PLANNER_LIB_H
#define GRADIENT_PLANNER_LIB_H

#include <iostream>
#include <eigen3/Eigen/Dense>	
#include <eigen3/Eigen/Sparse>
#include <math.h>

using namespace Eigen;
using namespace std;

/**
 * @brief definition of data type for easier sparse matrix initialization
 */
typedef Eigen::Triplet<double> T;

/**
 * @brief Structure of constants for algorithm
 * 
 * @param N number of time samples [-]
 * @param T sample time [s]
 * @param v car speed [m/s]
 * @param WB car wheelbase [m]
 * @param wr ratio between car position measuring point to rear axle and wheelbase \f$ l_r/WB \f$
 */

typedef struct params_s{
    int N; // number of samples
    double T; // sample time [s]
    double v; // constant speed [m/s]
    double WB; // wheelbase [m]
    double wr; // ratio between rear axle and wheelbase
} params;

/**
 * @brief Compute clothoid trajectory between two points
 * @details 
 * 
 * ### Algorithm description
 * 
 * This function implements Newton-type algorithm for solving of constrained optimization problem by following the gradient descent.
 * 
 * Optimization is constrained by discretized simplified car bicycle model:
 * 
 * \f$
 *  \begin{array}{ccc}
 *      x_{k+1} & = & v T + x_k
 *      \\ y_{k+1} & = & v \sin(\psi_k + c_k w_r W)\cdot T + y_k
 *      \\ \psi_{k+1} & = & v c_k \cdot T + \psi_k
 *      \\ c_{k+1} & = & u_k T + c_k
 *        \end{array}
 * \f$
 * 
 * Where
 *      - \f$ x,\ y\f$ is car position in m
 *      - \f$ \psi \f$ is car heading in rad
 *      - \f$ c \f$ is trajectory curvature in 1/m
 *      - \f$ v \f$ is car speed in m/s
 *      - \f$ T \f$ is discretization sampe time in s
 *      - \f$ W \f$ is car wheelbase in m
 *      - \f$ w_r \f$ is ratio between point where we measure the car position to the rear axle and wheelbase
 *      - \f$ u \f$ is curvature derivative [1/m/s], which is control input to the system
 *
 * This model is valid under the assumption that steering radius \f$ R \f$ is much greater than car wheelbase and speed is constant during the plan.
 * 
 * Relationship between steering angle \f$ \delta \f$ and curvature is \f$ \tan(\delta) = W/R \f$
 *
 * General optimization problem is defined as minimization of control input (curvature derivative) over the trajectory. By weighting of inputs, we can obtain clothoid shape of curve.
 * 
 * Inside the function, the Newton iteration is implemented
 * \f$ \left[\begin{array}{c} \mathbf{X}_{k+1} \\ \lambda_{k+1}\end{array}\right] = 
 *       \left[\begin{array}{c} \mathbf{X}_{k} \\ 0\end{array}\right] -
 *       \left[\begin{array}{cc} B_k & \nabla \mathbf{g}^T(\mathbf{X}_k) \\ \nabla \mathbf{g}(\mathbf{X}_k) & 0\end{array}\right]^{-1}
 *       \left[\begin{array}{c} \nabla L^T(\mathbf{X}_{k}) \\ \mathbf{g}(\mathbf{X}_{k}) \end{array}\right] \f$
 * and solved as \f$ x=A^{-1}b \f$
 * 
 * ### Implementation notes
 * 
 * This code uses [eigen](http://eigen.tuxfamily.org/dox/) library for matrix computation. Mainly sparse matrix structures are used.
 * 
 * @see [M. Diehl, S. Gros. Numerical Optimal Control](https://www.syscop.de/files/users/Jesus.lago/book.pdf)
 * @see pdf document, if you want more functionality, which is described in the document, ask Jiri Zahora to implement it.
 * @todo add link to pdf document
 * 
 * @param p         Structure with algorithm parameters
 * @param start_pt  Initial state of car in whatever coordinate frame.
 *                  State is defined by 4 numbers:
 *                  - x[m]
 *                  - y[m]
 *                  - heading[rad]
 *                  - curvature[1/m] 
 *                  
 * @param end_pt    Final state
 *                  in the same coordinate frame as initial state
 * @return Computed trajectory (in the same coordinate frame as initial state input) is 6 x (N+1) matrix with columns.
 *          - x [m]
 *          - y[ m]
 *          - heading [rad]
 *          - curvature [1/m]
 *          - curvature derivative [1/m/s]
 *          - time [s]
 *          
 */
MatrixXd compute_trajectory(params& p, Vector4d& start_pt, Vector4d& end_pt);

/**
 * @brief State space equation of simplified bicycle model. x[k+1] = f(x[k],u[k])
 * @param X car state [x,y,heading,curvature]
 * @param u control input [curvature derivative]
 * @param p Ssructure with algorithm parameters
 * @return car state in the next timestep
 */
Vector4d F(Vector4d& X, double u, params& p);

/**
 * @brief Compute vector of constraints g(X) = 0
 * @details Here is filled vector
 * \f$ \mathbf{g}(\mathbf{X}) = \left[\begin{array}{c} \mathbf{f}(\mathbf{x}_i,\mathbf{u}_i)-\mathbf{x}_{i+1} \\ \mathbf{x}_N - \mathbf{x}_f \end{array}\right],\ i\in\{0,1,\dots,N-1\} \f$
 * in the vector \f$ b \f$
 * @param X vector of all states, inputs and costates
 * @param x0 inital state (just [y,heading,cuvature] is relevant for optimization)
 * @param xf final state (just [y,heading,cuvature] is relevant for optimization)
 * @param p Ssructure with algorithm parameters
 * @param b output of function, vector which this function is filling
 */
void G(VectorXd& X, Vector3d& x0, Vector3d& xf, params& p, VectorXd& b);

/**
 * @brief Compute gradient of constraints
 * @details This function fill in the gradients of constraints in the matrix A. For timestep i the matrix is
 * \f$ \nabla \mathbf{g}(\mathbf{X})_i = 
 *       \left[\begin{array}{cccccc}
 *              1 & vT\cos(\psi_i+c_iw_rW) & vTw_rW\cos(\psi_i+c_iw_rW) & -1 & 0 & 0
 *           \\ 0 & 1 & vT & 0 & -1 & 0
 *           \\ 0 & 0 & 1 & 0 & 0 & -1
 *        \end{array}\right]
 * \f$
 * @param X vector of all states, inputs and costates
 * @param p Ssructure with algorithm parameters
 * @param gg output of function, part of matrix A which is filled in
 */
// gradient of constraints grad(g(x))
//MatrixXd grad_G(VectorXd& X, params& p);
void grad_G(VectorXd& X, params& p, std::vector<T>& gg);

// gradient of cost function
//VectorXd grad_f(VectorXd& X, VectorXd& W, params& p);

/**
 * @brief Set angle into the +-Pi range
 * @param phi angle which we want to modify
 */
void in_pi_range(double *phi);
/**
 * @brief Transformation into the car frame (init_state = [0,0,0,curvature])
 * @param x0 initial state
 * @param xf final state
 */
void transform_to_car_frame(Vector4d& x0, Vector4d& xf);

#endif
