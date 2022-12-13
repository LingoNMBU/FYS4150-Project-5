#include "wavebox.hpp"
using namespace std;

int main()
{

    bool test2 = true;
    if (test2)
    {
        //sim params
        double h1  = 0.1; //Stepsize space
        double dt1  = 2.5e-5; //Stepsize time
        double T1   = 0.002; //Max time      

        Wavebox catbox1 = Wavebox(h1, dt1, T1);

        int L1 = catbox1.L;
        arma::sp_cx_mat A1 = arma::sp_cx_mat(L1,L1);
        arma::sp_cx_mat B1 = arma::sp_cx_mat(L1,L1);

        catbox1.make_matrices(A1,B1);

        arma::cx_mat A2(A1);
        arma::cx_mat B2(B1);

        A2.save("A2");
        B2.save("B2");

        cout << "A1" << endl;
        cout << A2 << endl;
        cout << endl;
        cout << "B1" << endl;
        cout << B2 << endl;
        cout << endl;
    }

    bool test_coords = false;
    if (test_coords)
    {
        double h1   = 0.25; //Stepsize space
        double dt1  = 2.5e-5; //Stepsize time
        double T1 = 2.5e-5; //Stepsize time 

        int dim = 100;  

        arma::cx_mat mat1 = arma::cx_mat(dim,dim, arma::fill::randu);

        Wavebox catbox3 = Wavebox(h1, dt1, T1);


        arma::cx_vec vec1 = catbox3.convert_to_vector(mat1, dim);

        arma::cx_mat mat2 = catbox3.convert_to_matrix(vec1, dim*dim);

        arma::cx_vec vec2 = catbox3.convert_to_vector(mat2, dim);

        assert(arma::approx_equal(mat1, mat2, "absdiff", 0.0001));
        assert(arma::approx_equal(vec1, vec2, "absdiff", 0.0001));
    }

    bool test4 = false;
    if (test4)
    {
        double dt3 = 0.01;
        int M3 = 100;
        double h3 = 1./M3;
        double T3 = 1.0;
        double v_03 = 100. ;

        Wavebox catbox3 = Wavebox(h3, dt3, T3);

        double x_thick = 0.02;
        double x_center = 0.5;
        int n_slits = 1;
        std::vector<double> slit_widths;
        std::vector<double> wall_widths;

        slit_widths.push_back(0.05);
        wall_widths.push_back(0.05);

       catbox3.generate_slit_potential(x_thick ,x_center ,n_slits, slit_widths ,wall_widths, v_03);
       arma::mat V = catbox3.V;

        std::cout << std::endl << std::endl;
        std::cout << V << std::endl;
        std::cout << std::endl;
    }

    bool test_double_slit = false;
    if (test_double_slit)
    {
        //sim params
        double h   = 0.005; //Stepsize space
        double dt  = 2.5e-5; //Stepsize time
        double T   = 0.002; //Max time
        double v_03 = 1.e10 ;

        Wavebox catbox3 = Wavebox(h, dt, T);

        double x_thick = 0.02;
        double x_center = 0.5;
        int n_slits = 2;
        std::vector<double> slit_widths;
        std::vector<double> wall_widths;

        slit_widths.push_back(0.05);
        slit_widths.push_back(0.05);
        slit_widths.push_back(0.05);
        wall_widths.push_back(0.05);
        wall_widths.push_back(0.05);

        catbox3.generate_slit_potential(x_thick ,x_center ,1, slit_widths ,wall_widths, v_03);
        arma::mat V1 = catbox3.V;
        V1.save("V_single");

        catbox3.generate_slit_potential(x_thick ,x_center ,2, slit_widths ,wall_widths, v_03);
        arma::mat V2 = catbox3.V;
        V2.save("V_double");

        catbox3.generate_slit_potential(x_thick ,x_center ,3, slit_widths ,wall_widths, v_03);
        arma::mat V3 = catbox3.V;
        V3.save("V_triple");

    }


    bool test5 = false;
    
    if (test5)
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

        Wavebox catbox5 = Wavebox(h, dt, T);

        catbox5.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);
        arma::cx_vec u = catbox5.u;

        std::cout << std::endl << std::endl;
        std::cout << arma::size(u) << std::endl;
        std::cout << std::endl;

        std::cout << std::endl << std::endl;
        std::cout << arma::accu(arma::conj(u)%u)<< std::endl;
        std::cout << std::endl;

    }
    bool test6 = false;
    
    if (test6)
    {
        //sim params
        double h   = 0.005; //Stepsize space
        double dt  = 2.5e-5; //Stepsize time
        double T   = 0.0008; //Max time

        //Packet params
        double c_x = 0.25; //packet start center x
        double s_x = 0.05; //packet start std x
        double p_x = 200.; //packet start momentum x
        double c_y = 0.5; //packet start center y
        double s_y = 0.05; //packet start std y
        double p_y = 0.0; //packet start momentum y

        Wavebox catbox6 = Wavebox(h, dt, T);

        catbox6.initialize_packet(c_x, c_y, p_x, p_y, s_x, s_y);
        arma::cx_vec u = catbox6.u;

        std::cout << std::endl << std::endl;
        std::cout << arma::size(u) << std::endl;
        std::cout << std::endl;

        std::cout << std::endl << std::endl;
        std::cout << arma::accu(arma::conj(u)%u)<< std::endl;
        std::cout << std::endl;

        catbox6.simulate(false,true);

        arma::cx_vec p_sum = catbox6.p_sum;

        std::cout << std::endl << std::endl;
        std::cout << p_sum << std::endl;
        std::cout << std::endl;



    }
return 0;

}