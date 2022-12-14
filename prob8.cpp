#include "wavebox.hpp"
using namespace std;

int main()
{
    //sim params
    double h   = 0.005; //Stepsize space
    double dt  = 2.5e-5; //Stepsize time
    double T   = 0.002; //Max time

    //Packet params
    double c_x = 0.25; //packet start center x
    double s_x = 0.05; //packet start std x
    double p_x = 200.; //packet start momentum x
    double c_y = 0.5; //packet start center y
    double s_y = 0.2; //packet start std y
    double p_y = 0.0; //packet start momentum y

    //potential params
    double v_0 = 1.e10; //potential
    double x_thick = 0.02; //thicknes of barrier
    double x_center = 0.5; // symmetry center of barrier
    double slit_width1 = 0.05; //width slit
    double wall_width1 = 0.05; //width of wall separating slits
    int n_slits1 = 1;
    int n_slits2 = 2;
    int n_slits3 = 3;

    std::vector<double> slit_widths;
    std::vector<double> wall_widths;
    slit_widths.push_back(slit_width1);
    slit_widths.push_back(slit_width1);
    slit_widths.push_back(slit_width1);
    wall_widths.push_back(wall_width1);
    wall_widths.push_back(wall_width1);



    bool prob8_noV = true;
    if (prob8_noV)
    {
        Wavebox catbox8 = Wavebox(h, dt, T);

        catbox8.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);

        catbox8.simulate(true, false);

        arma::cx_cube u_cube = catbox8.u_cube;

        u_cube.save("prob8_u_cube_noV");
    }

    bool prob8_single = true;
    if (prob8_single)
    {
        Wavebox catbox8 = Wavebox(h, dt, T);

        catbox8.generate_slit_potential(x_thick, x_center, n_slits1 ,slit_widths ,wall_widths ,v_0);

        catbox8.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);

        catbox8.simulate(true, false);

        arma::cx_cube u_cube = catbox8.u_cube;

        u_cube.save("prob8_u_cube_single");
    }

    bool prob8_double = true;
    if (prob8_double)
    {
        Wavebox catbox8 = Wavebox(h, dt, T);

        catbox8.generate_slit_potential(x_thick, x_center, n_slits2 ,slit_widths ,wall_widths ,v_0);

        catbox8.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);

        catbox8.simulate(true, false);

        arma::cx_cube u_cube = catbox8.u_cube;

        u_cube.save("prob8_u_cube_double");
    }

    bool prob8_triple = true;
    if (prob8_triple)
    {
        Wavebox catbox8 = Wavebox(h, dt, T);

        catbox8.generate_slit_potential(x_thick, x_center, n_slits3 ,slit_widths ,wall_widths ,v_0);

        catbox8.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);

        catbox8.simulate(true, false);

        arma::cx_cube u_cube = catbox8.u_cube;

        u_cube.save("prob8_u_cube_triple");
    }
    


}