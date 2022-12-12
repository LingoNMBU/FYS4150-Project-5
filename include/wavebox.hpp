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
    arma::mat V;
    int M, L, m, N_steps;
    double T, dt, h, x_c, y_c, p_x, p_y, s_x,s_y, v_0;
    arma::cx_vec u;
    arma::cx_vec p_sum;
    arma::cx_cube u_cube;

    //constructor
    Wavebox(double h_in, double dt_in, double T_in);

    arma::ivec coords_vec_to_mat(int k, int m);

    int coords_mat_to_vec(int i, int j, int m);

    arma::cx_vec convert_to_vector(arma::cx_mat mat, int old_dim);

    arma::cx_mat convert_to_matrix(arma::cx_vec vec, int old_dim);

    void make_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B);

    void initialize_packet(double c_x, double c_y, double p_x, double p_y, double s_x, double s_y);

    void generate_slit_potential(double x_thick, double x_center, int n_slits, std::vector<double> slit_widths , std::vector<double> wall_widths, double v0);

    void simulate(bool store_u=false, bool store_p_sum=false);


};
#endif  