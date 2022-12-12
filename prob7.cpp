#include "wavebox.hpp"
using namespace std;

int main()
{
    //sim params
    double h   = 0.005; //Stepsize space
    double dt  = 2.5e-5; //Stepsize time
    double T   = 0.008; //Max time

    //Packet params
    double c_x = 0.25; //packet start center x
    double s_x = 0.05; //packet start std x
    double p_x = 200.; //packet start momentum x
    double c_y = 0.5; //packet start center y
    double s_y = 0.05; //packet start std y
    double s_y2 = 0.1; //packet start std y
    double p_y = 0.0; //packet start momentum y

    //potential params
    double v_0 = 1.e10; //potential
    double x_thick = 0.02;
    double x_center = 0.5;
    double slit_width1 = 0.05;
    double slit_width2 = 0.05;
    double wall_width1 = 0.05;
    int n_slits = 2;
    std::vector<double> slit_widths;
    std::vector<double> wall_widths;
    slit_widths.push_back(slit_width1);
    slit_widths.push_back(slit_width2);
    wall_widths.push_back(wall_width1);

    bool prob7a = true;
    if (prob7a)
    {
        Wavebox catbox7a = Wavebox(h, dt, T);

        catbox7a.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);

        catbox7a.simulate(false,true);

        arma::cx_vec p_sum = catbox7a.p_sum;

        p_sum.save("prob7a");
    }

    bool prob7b = true;
    if (prob7b)
    {
        Wavebox catbox7b = Wavebox(h, dt, T);

        catbox7b.generate_slit_potential(x_thick, x_center, n_slits ,slit_widths ,wall_widths ,v_0);

        catbox7b.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y2);

        catbox7b.simulate(false,true);

        arma::cx_vec p_sum = catbox7b.p_sum;

        p_sum.save("prob7b_psum");
    }
    


}