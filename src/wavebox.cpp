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

    u = arma::cx_vec(L); // init wavefunction vector
    V = arma::mat(M,M).fill(0.);

    //Initialize datastoring
    p_sum = arma::cx_vec(N_steps).fill(arma::cx_double(0,0));
    u_cube = arma::cx_cube(m, m, N_steps); 
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

    int i = k%m;
    int j = k/m;

    arma::ivec coords = { i, j };

    return coords;
}

arma::cx_mat Wavebox::convert_to_matrix(arma::cx_vec vec, int old_dim)
{

    int new_dim = sqrt(old_dim);
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

            vec(coord) = mat(i,j);
        }
    }
    return vec;
}

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
    std::cout << std::endl;
    std::cout << "made A and B matrices" << std::endl;
    std::cout << std::endl;

} 



void Wavebox::initialize_packet(double x_c, double y_c, double p_x, double p_y, double s_x, double s_y)
{
    
    for (int i=0; i < m; i++)
    {
        for (int j=0; j < m; j++)
        {
            //Boundary condition by only iterating over internal points -> border points are always zero

            double x = i*h;
            double y = j*h;


            double u_real = -pow(x-x_c,2) / (2.*s_x*s_x) - pow(y-y_c,2) / (2.*s_y*s_y);
            double u_imag = p_x*(x-x_c) + p_y*(y-y_c);

            int k = coords_mat_to_vec(j, i, m);

            u(k) = exp(arma::cx_double(u_real, u_imag));
        }
    }
    //Probability should sum to 1
    arma::cx_double probsum = arma::accu(arma::conj(u)%u); 
    //arma::cx_double probsum = sqrt(arma::accu(arma::real(u)%arma::real(u) + arma::imag(u)%arma::imag(u))); 

    u = u*sqrt(1./probsum); // A^2 *probsum = 1 -> A = sqrt(1/probsum) -> u* sqrt(1/probsum)

    std::cout << std::endl;
    std::cout << "initialized gaussian wavepacket at x, y = " << std::endl;
    std::cout << x_c << std::endl;
    std::cout << y_c << std::endl;
    std::cout << std::endl;

}

void Wavebox::generate_slit_potential(double x_thick, double x_center, int n_slits, std::vector<double> slit_widths , std::vector<double> wall_widths, double v0)
{
    arma::mat V_new = arma::mat(M,M).fill(0.);
    int mid_ind = round(M/2.);

    if (n_slits == 1)
    {
     int slit_start =  round((mid_ind*h - slit_widths[0] / 2.) / h);
     int slit_end   =  round((mid_ind*h + slit_widths[0] / 2.) / h);

     int slit_w_start   =  round((mid_ind*h - x_thick / 2.) / h);
     int slit_w_end     =  round((mid_ind*h + x_thick / 2.) / h);     
     
        for(int i = 0 ; i < M ; i++)
            if (i < slit_start or i > slit_end)
            {
                for (int j = slit_w_start; j < slit_w_end+1; j++)
                {
                    V_new(i, j) = v0;
                }           
            }
        std::cout << std::endl;
        std::cout << "mid" << std::endl;
        std::cout << mid_ind << std::endl;
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << "start1" << std::endl;
        std::cout << slit_start << std::endl;
        std::cout << std::endl;

        std::cout << "end1" << std::endl;
        std::cout << slit_end << std::endl;
        std::cout << std::endl;
    }

    if (n_slits == 2)
    {
        int slit1_start =  round((mid_ind*h - wall_widths[0] / 2. - slit_widths[0]) / h);
        int slit1_end   =  round((mid_ind*h - wall_widths[0] / 2.) / h);

        int slit2_start =  round((mid_ind*h + wall_widths[0] / 2.) / h);
        int slit2_end   =  round((mid_ind*h + wall_widths[0] / 2. + slit_widths[1]) / h);

        int slit_w_start   =  round((mid_ind*h - x_thick / 2.) / h);
        int slit_w_end     =  round((mid_ind*h + x_thick / 2.) / h);     

        std::cout << std::endl;
        std::cout << "mid" << std::endl;
        std::cout << mid_ind << std::endl;
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << "start1" << std::endl;
        std::cout << slit1_start << std::endl;
        std::cout << std::endl;

        std::cout << "end1" << std::endl;
        std::cout << slit1_end << std::endl;
        std::cout << std::endl;

        std::cout << "start2" << std::endl;
        std::cout << slit2_start << std::endl;
        std::cout << std::endl;

        std::cout << "end2" << std::endl;
        std::cout << slit2_end << std::endl;
        std::cout << std::endl;
     
        for(int i = 0 ; i < M ; i++)
        {
            if ((i < slit1_start) or (i > slit1_end and i < slit2_start) or (i > slit2_end))
            {
                for (int j = slit_w_start; j < slit_w_end+1; j++)
                {
                    V_new(i, j) = v0;
                }           
            }
        }
    }        


    if (n_slits == 3)
    {
        int slit1_start =  round((mid_ind*h - slit_widths[1] / 2. - wall_widths[0] - slit_widths[0]) / h);
        int slit1_end   =  round((mid_ind*h - slit_widths[1] / 2. - wall_widths[0]) / h);

        int slit2_start =  round((mid_ind*h - slit_widths[1] / 2) / h);
        int slit2_end   =  round((mid_ind*h + slit_widths[1] / 2) / h);

        int slit3_start =  round((mid_ind*h + slit_widths[1] / 2. + wall_widths[1]) / h);
        int slit3_end   =  round((mid_ind*h + slit_widths[1] / 2. + wall_widths[1] + slit_widths[2]) / h);

        //int wall1_start  =  round((mid_ind*h - slit_widths[1] / 2. - wall_widths[0]) / h);
        //int wall1_end    =  round((mid_ind*h - slit_widths[1] / 2.) / h);     

        //int wall2_start  =  round((mid_ind*h + slit_widths[1] / 2.) / h);
        //int wall2_end    =  round((mid_ind*h + slit_widths[2] / 2. + wall_widths[1]) / h);     

        int slit_w_start   =  round((mid_ind*h - x_thick / 2.) / h);
        int slit_w_end     =  round((mid_ind*h + x_thick / 2.) / h);

        std::cout << std::endl;
        std::cout << "mid" << std::endl;
        std::cout << mid_ind << std::endl;
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << "start1" << std::endl;
        std::cout << slit1_start << std::endl;
        std::cout << std::endl;

        std::cout << "end1" << std::endl;
        std::cout << slit1_end << std::endl;
        std::cout << std::endl;

        std::cout << "start2" << std::endl;
        std::cout << slit2_start << std::endl;
        std::cout << std::endl;

        std::cout << "end2" << std::endl;
        std::cout << slit2_end << std::endl;
        std::cout << std::endl;

        std::cout << "start3" << std::endl;
        std::cout << slit3_start << std::endl;
        std::cout << std::endl;

        std::cout << "end3" << std::endl;
        std::cout << slit3_end << std::endl;
        std::cout << std::endl;

        for(int i = 0 ; i < M ; i++)
        {
            if ((i < slit1_start) or (i > slit1_end and i < slit2_start) or (i > slit2_end and i < slit3_start) or (i > slit3_end))
            {
                for (int j = slit_w_start; j < slit_w_end+1; j++)
                {
                    V_new(i, j) = v0;
                }           
            }
        }
    }

    //Store V
    V = V_new;

    std::cout << std::endl;
    std::cout << "set up potential with n_slits =" << std::endl;
    std::cout << n_slits << std::endl;
    std::cout << std::endl;

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

    std::cout << std::endl;
    std::cout << "Started simulations, number of steps is" << std::endl;
    std::cout << N_steps << std::endl;
    std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    //double T_curr = 0;
    for (int t=1; t < N_steps; t++)
    {
        if (t%10 == 0)
        {
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

            std::cout << "duration [s]" << std::endl;
            std::cout << duration.count() << std::endl;
            std::cout << std::endl;

            double ETA = (duration.count()/t)*(N_steps-t);

            std::cout << "ETA [s]" << std::endl;
            std::cout << ETA << std::endl;
            std::cout << std::endl;
        }
        //step 1, matrix multiplication
        arma::cx_vec b = B*u;

        //step 2, solve equations
        arma::cx_vec u1 = arma::spsolve(A,b);

        //Update values
        u = u1;

        //Store values
        if (store_p_sum)
        {
            p_sum(t) = arma::accu(arma::conj(u1)%u1);
        }

        if (store_u)
        {
            arma::cx_mat u_mat = convert_to_matrix(u,L);
            u_cube.slice(t) = u_mat;
        }
    }
}

#endif 