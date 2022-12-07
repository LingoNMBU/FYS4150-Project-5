#ifndef __wavebox_hpp__  
#define __wavebox_hpp__


#include <armadillo>
#include <random>
#include <chrono>
#include <algorithm>

class Wavebox
{
  private:

  public:
    arma::mat V;
    int M, L, m;
    double T, dt, h, x_c, y_c, p_x, p_y, s_x,s_y, v_0;
    arma::cx_vec u;

    //constructor
    Wavebox(int M_in, arma::mat V_in, double dt_in, double T_in);

    arma::ivec coords_vec_to_mat(int k, int L);

    int coords_mat_to_vec(int i, int j, int L);

    void make_matrices(arma::cx_mat& A, arma::cx_mat& B);

    void initialize_packet(double c_x, double c_y, double p_x, double p_y, double s_x, double s_y);

    void generate_slit_potential(double x_thick, double x_center, int n_slits, std::vector<double> slit_widths , std::vector<double> wall_widths);

    void simulate(void);


};
#endif  