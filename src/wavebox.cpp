#ifndef __wavebox_cpp__  
#define __wavebox_cpp__

#include "wavebox.hpp"

// Constructor
Wavebox::Wavebox(int M_in, arma::mat V_in, double dt_in, double T_in)
{

  V = V_in; //potential

  M = M_in; //number of points along each axis
  m = M-2; //Number of internal points
  L = m*m; //number of elements in internal matrices
  h = 1./(M-1); //stepsize

  dt = dt_in; //timestepsize
  T = T_in; //max time

  u = arma::cx_vec(L); // init wavefunction vector



}

// Change the spin of a single particle in the lattice
int Wavebox::coords_mat_to_vec(int i, int j, int L)
{
    //Translates two matrix coordinates into one vector coordinate
    // 0,0  0,1  0,2
    // 1,0  1,1  1,2
    // 2,0  2,1  2,2

    //0    1    2    3    4    5    6    7    8
    //0,0  1,0  2,0  0,1  1,1  2,1  0,2  1,2  2,2

    //L=3

    int k = i + j*L;
    return k;
}

arma::ivec Wavebox::coords_vec_to_mat(int k, int L)
{
    //Transaltes one vector coordinate into two matrix coordinates
    // 0,0  0,1  0,2
    // 1,0  1,1  1,2
    // 2,0  2,1  2,2

    //0    1    2    3    4    5    6    7    8
    //0,0  1,0  2,0  0,1  1,1  2,1  0,2  1,2  2,2

    //L=3

    int i = k%L;
    int j = k/L;

    arma::ivec coords = { i, j };

    return coords;
}


void Wavebox::make_matrices(arma::cx_mat& A, arma::cx_mat& B)
{
    //Constructs the A and B matrices
    double r = dt/(2*h*h);

    arma::cx_vec a = arma::cx_vec(L);
    arma::cx_vec b = arma::cx_vec(L);

    for (int k = 0 ; k < L ; k++)
    {
        arma::ivec coords = coords_vec_to_mat(k, m);
        
        double real = 1;
        double imag_a = 4*r + (dt/2)*V(coords(0), coords(1));
        double imag_b = 4*r - (dt/2)*V(coords(0), coords(1));

        a(k) = arma::cx_double(real, imag_a);
        b(k) = arma::cx_double(real, imag_b);
    }

    //THis is mad, there must be a better way if extra time or super slow
    for (int i=0 ; i<L ; i++)
    {
        for(int j=0 ; j < L ; j++)
        {
            //If on main diagonal
            if (i==j)
            {
                A(i,j) = a(i);
                B(i,j) = b(i);
            }

            //If on third diagonal 
            if(abs(j-i) == 3)
            {
                A(i,j) = -r;
                B(i,j) = r;
            }
            //If on first diagonal
            if(abs(i-j) == 1)
            {   
                //If on the element between the M-2 submatrices, aka boundary condition
                if((i%m == 0 and j%m == m-1) or (i%m == m-1 and j%m == 0))
                {

                    A(i,j) = arma::cx_double(0, 0);
                    B(i,j) = arma::cx_double(0, 0);
                }
                //Else inside the submatrices
                else
                {
                    A(i,j) =-r;
                    B(i,j) = r;
                }
            }
        }
    }
} 

void Wavebox::initialize_packet(double x_c, double y_c, double p_x, double p_y, double s_x, double s_y)
{
    
    for (int i=0; i < m; i++)
    {
        for (int j=0; j < m; j++)
        {
            //Boundary condition by only iterating over internal points -> border points are always zero
            //if (y*x == 0)
            //{
            //    u(k) = arma::cx_double(0.0, 0.0)
            //}
            //if (i == m or j == m)

            double x = i*h;
            double y = j*h;


            double pos = -(x-x_c)*(x-x_c)/(2*s_x*s_x) - (y-y_c)*(y-y_c)/(2*s_y*s_y);
            double mom = p_x*(x-x_c) + p_y*(y-y_c);

            int k = coords_mat_to_vec(i, j, m);

            //std::cout <<k<< std::endl;
            //std::cout <<i<< std::endl;
            //std::cout <<j<< std::endl;
            //std::cout <<m<< std::endl;
            //std::cout <<L<< std::endl;


            
            u(k) = exp(arma::cx_double(pos, mom));
        }
    }
    //Probability ums to 1
    std::complex<double> probsum = arma::accu(arma::conj(u)%u);
    u = u/probsum;




}

void Wavebox::generate_slit_potential(double x_thick, double x_center, int n_slits, std::vector<double> slit_widths , std::vector<double> wall_widths)
{
    arma::mat V_new = arma::mat(m,m).fill(0);
    int mid_ind = round(m/2);
    double v0 = 1.0 * 10e6;

    if (n_slits == 1)
    {
     int slit_start =  round((mid_ind*h - slit_widths[0] / 2.) / h);
     int slit_end   =  round((mid_ind*h + slit_widths[0] / 2.) / h);

     int slit_w_start   =  round((mid_ind*h - x_thick / 2.) / h);
     int slit_w_end     =  round((mid_ind*h + x_thick / 2.) / h);

     //std::cout << std::endl;
     //std::cout << slit_start << std::endl;
     //std::cout << std::endl;
     //std::cout << slit_end << std::endl;
     //std::cout << std::endl;
     //std::cout << slit_w_start << std::endl;
     //std::cout << std::endl;
     //std::cout << slit_w_end << std::endl;
     //std::cout << std::endl;
     //std::cout << mid_ind*h - x_thick/2.0 << std::endl;
     //std::cout << std::endl;
     //std::cout << (mid_ind*h - x_thick/2.0)/h << std::endl;     
     //std::cout << std::endl;
     //std::cout << (mid_ind*h + x_thick/2.0)/h << std::endl;     
     
     
     for(int i = 0 ; i < m ; i++)
        if (i < slit_start or i > slit_end)
        {
            for (int j = slit_w_start; j < slit_w_end; j++)
            {
                V_new(i, j) = v0;
            }           
        }
    }
    V = V_new;
    //return V;
}

void Wavebox::simulate(void)
{
    //1. Create potential

    //TODO, find elegant way to initialize potential from all dem variables
    //for now, assume already set up

    //2. Initialize wavepacket

    //TODO, find elegant way to initialize packet from all dem variables
    //for now, assume already set up 

    //2. Initialize and create Crank-Nicholson matrices
    arma::cx_mat A = arma::cx_mat(L,L).fill(arma::cx_double(0, 0));
    arma::cx_mat B = arma::cx_mat(L,L).fill(arma::cx_double(0, 0));
    make_matrices(A,B);


    double T_curr = 0;
    while (T_curr < T)
    {
        //step 1, matrix multiplication
        arma::cx_vec b = B*u;

        //step 2, solve equations
        arma::cx_vec u1 = arma::solve(A,b);

        //Update values
        u = u1;

        //Move time
        T_curr += dt;
    }

}





#endif 