#ifndef __wavebox_hpp__  
#define __wavebox_hpp__


#include <armadillo>
#include <random>
#include <chrono>
#include <math.h>
#include <assert.h>

class Wavebox
{
  private:

  public:
    arma::mat V;  //Stores the potential matrix
    int M, L, m, N_steps;  //Matrix size, Matrix-vector length, internal matrix, simulation steps
    double T, dt, h, x_c, y_c, p_x, p_y, s_x, s_y, v_0;  //Sim time interval, time step, space step, parameters for wave packet and potential
    arma::cx_vec u; // wave function vector
    arma::cx_vec p_sum; //probability conservation vector
    arma::cx_cube u_cube; //storage for wavefunction matrices at each timestep

    //constructor
    Wavebox(double h_in, double dt_in, double T_in);

    //Converts from vector to matrix coordinates
    arma::ivec coords_vec_to_mat(int k, int m);
    
    //Converts from matrix to vector coordinates
    int coords_mat_to_vec(int i, int j, int m);

    //Converts a matrix into a vector
    arma::cx_vec convert_to_vector(arma::cx_mat mat, int old_dim);

    //Converts a vector into a matrix
    arma::cx_mat convert_to_matrix(arma::cx_vec vec, int old_dim);

    //Makes A and B matrices for the Crank-Nicolson method
    void make_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B);

    //Initializes a gaussian wave packet 
    void initialize_packet(double c_x, double c_y, double p_x, double p_y, double s_x, double s_y);

    //Initializes a slit potential in a box
    void generate_slit_potential(double x_thick, double x_center, int n_slits, std::vector<double> slit_widths , std::vector<double> wall_widths, double v0);

    //Simulates the propagation of the Schroedinger equation through time and stores the results
    void simulate(bool store_u=false, bool store_p_sum=false);


};
#endif  