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
    double p_y = 0.0; //packet start momentum y

    //potential params
    double v_0 = 0.0; //potential
    double x_thick = 0.02;
    double x_center = 0.5;
    double slit_width1 = 0.05;
    double wall_width1 = 0.05;
    int n_slits = 1;
    std::vector<double> slit_widths;
    std::vector<double> wall_widths;
    slit_widths.push_back(0.05);
    wall_widths.push_back(0.05);

    bool prob7 = true;
    if (prob7)
    {
        Wavebox catbox7 = Wavebox(h, dt, T);

        catbox7.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);
        arma::cx_vec u = catbox7.u;

        catbox7.simulate(false,true);

        arma::cx_vec p_sum = catbox7.p_sum;

        std::cout << std::endl << std::endl;
        std::cout << p_sum << std::endl;
        std::cout << std::endl;

        p_sum.save("prob7");
    }

    


}