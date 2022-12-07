#include <armadillo>
#include <string>
#include <iostream>

int coords_mat_to_vec(int i, int j, int L)
{

    // 0,0  0,1  0,2
    // 1,0  1,1  1,2
    // 2,0  2,1  2,2

    //0    1    2    3    4    5    6    7    8
    //0,0  1,0  2,0  0,1  1,1  2,1  0,2  1,2  2,2

    //L=3

    int k = i + j*L;
    return k;
}

arma::ivec coords_vec_to_mat(int k, int L)
{

    // 0,0  0,1  0,2
    // 1,0  1,1  1,2
    // 2,0  2,1  2,2

    //0    1    2    3    4    5    6    7    8
    //0,0  1,0  2,0  0,1  1,1  2,1  0,2  1,2  2,2

    //L=3

    int i = k%L;
    int j = k/L;
    //int j2 = (k - i)/L;

    //std::cout << std::endl;
    //std::cout << "i" <<std::endl;
    //std::cout << i << std::endl;
    //std::cout << std::endl;
    //std::cout << "j" <<std::endl;
    //std::cout << j << std::endl;
    //std::cout << "j2" << std::endl;
    //std::cout << j2 << std::endl;
    //std::cout << std::endl;

    arma::ivec coords = { i, j };


    return coords;
}



void make_matrices(int M, double h, double dt, arma::mat V, arma::cx_mat& A, arma::cx_mat& B)
{
    int L = (M-2) * (M-2); //M-2 * M-2
    int m = M-2;
    double r = dt/(2*h*h);

    //std::cout << L << std::endl;
    //std::cout << std::endl;
    //std::cout << m << std::endl;
    //std::cout << std::endl;

    arma::cx_vec a = arma::cx_vec(L);
    arma::cx_vec b = arma::cx_vec(L);

    for (int k = 0 ; k < L ; k++)
    {
        arma::ivec coords = coords_vec_to_mat(k, m);
        
        double real = 1;
        double imag_a = 4*r + (dt/2)*V(coords(0), coords(1));
        double imag_b = 4*r - (dt/2)*V(coords(0), coords(1));

        //std::cout << std::endl;
        //std::cout << "__LINE__" << std::endl;
        //std::cout << std::endl;
        //std::cout << std::endl;
        //std::cout << k << std::endl;
        //std::cout << std::endl;

        a(k) = arma::cx_double(real, imag_a);
        b(k) = arma::cx_double(real, imag_b);
    }

    //std::cout << std::endl;
    //std::cout << "__LINE__" << std::endl;
    //std::cout << std::endl;

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
