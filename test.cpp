#include "wavebox.hpp"
#include "functions.hpp"
using namespace std;

int main()
{

    bool test1 = false;

    if (test1)
    {
        double dt = 0.1;
        double h = 0.1;
        int M = 5;
        int m = M-2;
        int L = m * m;

        arma::mat V = arma::mat(m,m, arma::fill::randu);
        arma::cx_mat A = arma::cx_mat(L,L).fill(arma::cx_double(0, 0));
        arma::cx_mat B = arma::cx_mat(L,L).fill(arma::cx_double(0, 0));

        make_matrices( M, h, dt, V, A, B);

        cout << "A" << endl;
        cout << A << endl;
        cout << endl;
        cout << "B" << endl;
        cout << B << endl;
        cout << endl;
    }

    bool test2 = false;
    if (test2)
    {
        double dt1 = 0.1;
        double h1 = 0.1;
        int M1 = 5;
        double T1 = 10.0;
        int L1 = (M1-2) * (M1-2);

        arma::mat V1 = arma::mat(M1-2,M1-2, arma::fill::randu);

        arma::cx_mat A1 = arma::cx_mat(L1,L1).fill(arma::cx_double(0, 0));
        arma::cx_mat B1 = arma::cx_mat(L1,L1).fill(arma::cx_double(0, 0));

        Wavebox catbox1 = Wavebox(M1, V1, dt1, T1);

        catbox1.make_matrices(A1,B1);

        cout << "A1" << endl;
        cout << A1 << endl;
        cout << endl;
        cout << "B1" << endl;
        cout << B1 << endl;
        cout << endl;
    }


        bool test3 = true;
    if (test3)
    {
        double dt2 = 0.01;
        double h2 = 0.01;
        int M2 = 5;
        double T2 = 1.0;
        int L2 = (M2-2) * (M2-2);

        double x_c2 = 0.5; 
        double y_c2 = 0.5;
        double p_x2 = 0.1;
        double p_y2 = 0.1;
        double s_x2 = 0.05;
        double s_y2 = 0.05;

        arma::mat V2 = arma::mat(M2-2,M2-2, arma::fill::randu);

        Wavebox catbox2 = Wavebox(M2, V2, dt2, T2);

        catbox2.initialize_packet(x_c2, y_c2, p_x2, p_y2, s_x2, s_y2);

        catbox2.simulate();
    }

    bool test4 = true;
    if (test4)
    {
        double dt3 = 0.01;
        int M3 = 100;
        double T3 = 1.0;
        int L3 = (M3-2) * (M3-2);

        arma::mat V3 = arma::mat(M3-2,M3-2, arma::fill::randu);

        Wavebox catbox3 = Wavebox(M3, V3, dt3, T3);

        double x_thick = 0.02;
        double x_center = 0.5;
        int n_slits = 1;
        std::vector<double> slit_widths;
        std::vector<double> wall_widths;

        slit_widths.push_back(0.05);
        wall_widths.push_back(0.05);

       catbox3.generate_slit_potential(x_thick ,x_center ,n_slits, slit_widths ,wall_widths);
       arma::mat V = catbox3.V;

        std::cout << std::endl << std::endl;
        std::cout << V << std::endl;
        std::cout << std::endl;
    }

    bool test5 = false;
    if (test5)
    {
        double dt5 = 0.01;
        int M5 = 100;
        double T5 = 1.0;
        int L5 = (M5-2) * (M5-2);

        double x_thick = 0.02;
        double x_center = 0.5;
        int n_slits = 1;
        std::vector<double> slit_widths;
        std::vector<double> wall_widths;

        slit_widths.push_back(0.05);
        wall_widths.push_back(0.05);


        arma::mat V5 = arma::mat(M5-2,M5-2, arma::fill::randu);

        Wavebox catbox5 = Wavebox(M5, V5, dt5, T5);
        
        catbox5.generate_slit_potential(x_thick ,x_center ,n_slits, slit_widths ,wall_widths);
        arma::mat V = catbox5.V;

        std::cout << std::endl << std::endl;
        std::cout << V << std::endl;
        std::cout << std::endl;
    }
return 0;

}