#ifndef __wavebox_cpp__  
#define __wavebox_cpp__

#include "wavebox.hpp"

// Constructor
Wavebox::Wavebox(double h_in, double dt_in, double T_in)
{
    dt = dt_in; //timestepsize
    T = T_in; //max time
    h = h_in; //stepsize
    M = round(1. + 1./h);
    m = M-2; //Number of internal points
    L = m*m; //number of elements in internal matrices
    N_steps = round(T/dt);
    //std::cout << std::endl << std::endl;
    //std::cout << round(1. + 1./h_in) << std::endl;
    //std::cout << std::endl;

    u = arma::cx_vec(L); // init wavefunction vector
    V = arma::mat(M,M).fill(0.);

    //Initialize datastoring
    p_sum = arma::cx_vec(N_steps).fill(arma::cx_double(0,0));
    //us = arma::cx_mat(N_steps, L); //might cause memory issues

}

// Change the spin of a single particle in the lattice
int Wavebox::coords_mat_to_vec(int i, int j, int m)
{
    //Translates two matrix coordinates into one vector coordinate
    // 0,0  0,1  0,2
    // 1,0  1,1  1,2
    // 2,0  2,1  2,2

    //0    1    2    3    4    5    6    7    8
    //0,0  1,0  2,0  0,1  1,1  2,1  0,2  1,2  2,2

    //m=3

    //TODO test this function
    //TODO remove L as variable
    int k = j*m + i;
    return k;
}

arma::ivec Wavebox::coords_vec_to_mat(int k, int m)
{
    //Transaltes one vector coordinate into two matrix coordinates
    // 0,0  0,1  0,2
    // 1,0  1,1  1,2
    // 2,0  2,1  2,2

    //0    1    2    3    4    5    6    7    8
    //0,0  1,0  2,0  0,1  1,1  2,1  0,2  1,2  2,2

    //m=3

    //TODO test this function
    //TODO remove L as variable

    int i = k%m;
    int j = k/m;

    arma::ivec coords = { i, j };

    return coords;
}

arma::cx_mat Wavebox::convert_to_matrix(arma::cx_vec vec, int old_dim)
{
    int new_dim = int(sqrt(old_dim));
    arma::cx_mat mat = arma::cx_mat(new_dim, new_dim);

    for (int i=0 ; i < old_dim; i++)
    {
        arma::ivec coords = coords_vec_to_mat(i, new_dim);
        mat(coords(0),coords(1)) = vec(i);
    }
    return mat;
}

arma::cx_vec Wavebox::convert_to_vector(arma::cx_mat mat, int old_dim)
{
    int new_dim = (old_dim*old_dim);
    arma::cx_vec vec = arma::cx_vec(new_dim);

    for (int i=0 ; i < old_dim; i++)
    {
        for (int j=0; j < old_dim; j++)
        {
            int coord = coords_mat_to_vec(i, j, old_dim);
            //std::cout << std::endl;
            //std::cout << i << std::endl;
            //std::cout << std::endl;

            //std::cout << std::endl;
            //std::cout << j << std::endl;
            //std::cout << std::endl;

            //std::cout << std::endl;
            //std::cout << coord << std::endl;
            //std::cout << std::endl;

            vec(coord) = mat(i,j);
        }
    }
    return vec;
}


//void Wavebox::make_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B)
//{
//    //Constructs the A and B matrices
//    double r = dt/(2.*h*h);
//
//    arma::cx_vec a = arma::cx_vec(L);
//    arma::cx_vec b = arma::cx_vec(L);
//
//    for (int k = 0 ; k < L ; k++)
//    {
//        arma::ivec coords = coords_vec_to_mat(k, m);
//        
//        double real = 1.;
//        double imag_a =  4.*r + (dt/2)*V(coords(0), coords(1));
//        double imag_b = -4.*r - (dt/2)*V(coords(0), coords(1));
//
//        a(k) = arma::cx_double(real, imag_a);
//        b(k) = arma::cx_double(real, imag_b);
//    }
//
//    //THis is mad, there must be a better way if extra time or super slow
//    for (int i=0 ; i<L ; i++)
//    {
//        for(int j=0 ; j < L ; j++)
//        {
//            //If on main diagonal
//            if (i==j)
//            {
//                A(i,j) = a(i);
//                B(i,j) = b(i);
//            }
//
//            //If on third diagonal 
//            if(abs(j-i) == 3)
//            {
//                A(i,j) = -r;
//                B(i,j) = r;
//            }
//            //If on first diagonal
//            if(abs(i-j) == 1)
//            {   
//                //If on the element between the M-2 submatrices, aka boundary condition
//                if((i%m == 0 and j%m == m-1) or (i%m == m-1 and j%m == 0))
//                {
//                    //Assumed 0 from start by sparse matrix
//                    //A(i,j) = arma::cx_double(0, 0);
//                    //B(i,j) = arma::cx_double(0, 0);
//                }
//                //Else inside the submatrices
//                else
//                {
//                    A(i,j) =-r;
//                    B(i,j) = r;
//                }
//            }
//        }
//    }
//} 

void Wavebox::make_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B)
{
    arma::cx_double r = arma::cx_double(0, dt/(2.*h*h));
    
    arma::cx_vec a = arma::cx_vec(L);
    arma::cx_vec b = arma::cx_vec(L);

    for (int k=0; k < L; k++)
    {
        arma::ivec coords = coords_vec_to_mat(k, m);

        //plus one since we are only iterating over internal points and V is full size
        arma::cx_double a_imag =  4.*r + arma::cx_double(0, (dt/2)*V(coords(0)+1,coords(1)+1));
        arma::cx_double b_imag = -4.*r - arma::cx_double(0, (dt/2)*V(coords(0)+1,coords(1)+1));

        a(k) = arma::cx_double(1., 0.) + a_imag;
        b(k) = arma::cx_double(1., 0.) + b_imag;
    }

        //diagonals
    for (int i=0; i<L; i++)
    {

        A(i,i) = a(i);
        B(i,i) = b(i);
    }
    A.diag( 1).fill(-r);
    A.diag(-1).fill(-r);
    A.diag( m).fill(-r);
    A.diag(-m).fill(-r);

    B.diag( 1).fill(r);
    B.diag(-1).fill(r);
    B.diag( m).fill(r);
    B.diag(-m).fill(r);

    for (int i = m-1; i < L-2; i+=3)
    {
        //Setting boundary condition zeroes
        A.diag( 1)(i) = 0;
        A.diag(-1)(i) = 0;
        B.diag( 1)(i) = 0;
        B.diag(-1)(i) = 0;
    }
    std::cout << "made dem matrices" << std::endl;
    std::cout << "made dem matrices" << std::endl;
    std::cout << "made dem matrices" << std::endl;

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

            if (i == 0 || i == m - 1 || j == 0 || j == m - 1)
            {
                continue;
            }

            double x = i*h;
            double y = j*h;


            double u_real = -pow(x-x_c,2) / (2.*s_x*s_x) - pow(y-y_c,2) / (2.*s_y*s_y);
            double u_imag = p_x*(x-x_c) + p_y*(y-y_c);

            int k = coords_mat_to_vec(i, j, m);

            u(k) = exp(arma::cx_double(u_real, u_imag));
        }
    }
    //Probability ums to 1
    arma::cx_double probsum = arma::accu(arma::conj(u)%u); 
    //arma::cx_double probsum = sqrt(arma::accu(arma::real(u)%arma::real(u) + arma::imag(u)%arma::imag(u))); 

    u = u*sqrt(1./probsum); // A^2 *probsum = 1 -> A = sqrt(1/probsum) -> u* sqrt(1/probsum)

    //std::cout << std::endl << std::endl;
    //std::cout << arma::accu(arma::conj(u)%u)<< std::endl;
    //std::cout << std::endl;


        


}

void Wavebox::generate_slit_potential(double x_thick, double x_center, int n_slits, std::vector<double> slit_widths , std::vector<double> wall_widths, double v0)
{
     //std::cout << std::endl;
     //std::cout << M << std::endl;
     //std::cout << std::endl;


    arma::mat V_new = arma::mat(M,M).fill(0);
    int mid_ind = round(M/2);

    if (n_slits == 1)
    {
     int slit_start =  round((mid_ind*h - slit_widths[0] / 2.) / h);
     int slit_end   =  round((mid_ind*h + slit_widths[0] / 2.) / h);

     int slit_w_start   =  round((mid_ind*h - x_thick / 2.) / h);
     int slit_w_end     =  round((mid_ind*h + x_thick / 2.) / h);

     
     
     for(int i = 0 ; i < M ; i++)
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



void Wavebox::simulate(bool store_u, bool store_p_sum)
{
    //1. Create potential

    //TODO, find elegant way to initialize potential from all dem variables
    //for now, assume already set up

    //2. Initialize wavepacket

    //TODO, find elegant way to initialize packet from all dem variables
    //for now, assume already set up 

    //3. Initialize and create Crank-Nicholson matrices

    arma::sp_cx_mat A = arma::sp_cx_mat(L,L);
    arma::sp_cx_mat B = arma::sp_cx_mat(L,L);
    make_matrices(A,B);


    if (store_p_sum)
    {
        p_sum(0) = arma::accu(arma::conj(u)%u);
    }

    //double T_curr = 0;
    for (int t=1; t < N_steps; t++)
    {
        //step 1, matrix multiplication
        arma::cx_vec b = B*u;

        //step 2, solve equations
        arma::cx_vec u1 = arma::spsolve(A,b);


        //TODO INVESTIGATE A AND B, check that eigenvalues are less or equal to 1, check that the matrices are correct

        //Update values
        u = u1;

        //Store values
        if (store_p_sum)
        {
            p_sum(t) = arma::accu(arma::conj(u1)%u1);
        }

        //if (store_u)
        //{
        //    us(t, arma::span::all) = u;
        //}

        //Move time
        //T_curr += dt;
    }

}





#endif 