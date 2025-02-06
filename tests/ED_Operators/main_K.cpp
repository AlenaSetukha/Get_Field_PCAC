#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "Kernel_Par.h"
#include "K.h"
#include "K0.h"
#include "R.h"
#include "Num_Par.h"
#include "element_geom.h"
#include "kernel_lib.h"


/**
 * Тест для проверки вычисления основных интегральных операторов задач
 * рассеяния: K = rot rot
 * Порядок обхода ячеек важен: либо по, либо против часовой
 * стрелки - все одновременно.
 * Ток должен лежить в плоскости ячейки(если в операторе присутствует
 * контурный интеграл)
 */



int main(int argc, char **argv)
{
    //---------Тестирование операторов--------------
    /*
    const double rut1[4][3] = {{0., 0., 0.},
                         {2., 0., 0.},
                         {3., 2., 0.},
                         {0., 2., 0.}};

    const double rut2[3][3] = {{0., 0., 0.},
                         {2., 0., 0.},
                         {0., 2., 0.}};
    
    
    const double rut3[3][3] = {{3., 2., 0.},
                         {0., 2., 0.},
                         {2., 0., 0.}};


    const double x[3] = {3., 3., 0.};
    const std::complex<double> j[3] = {1., 1., 0.,};
*/

    const double rut1[4][3] = {{1., 1., 0.},
                         {2., 0., 0.},
                         {2., 2., 0.},
                         {0., 2., 0.}};

    const double rut2[3][3] = {{1., 1., 0.},
                         {2., 0., 0.},
                         {2., 2., 0.}};
    
    
    const double rut3[3][3] = {{1., 1., 0.},
                         {2., 2., 0.},
                         {0., 2., 0.}};


    const double x[3] = {4., 1., 1.};


    const std::complex<double> j[3] = {1., 0., 0.,};

    //-------------------Test1---------------------
    std::cout << "Test1. K при наборе матрицы, сглаживание." << std::endl;
    const std::complex<double> k = std::complex<double>(1., 0.);

    Num_Par num_param(0.00001, 0.5, 0.5, 8, 10, 3, 3);


    std::complex<double> res[3];
    K_rot_rot(j, x, rut1, num_param, k, res);
    std::cout << "Res1: " << res[0] << " " << res[1] << " " << res[2] << std::endl;


    std::complex<double> res1[3], res2[3];
    K_rot_rot(j, x, rut2, num_param, k, res1);
    K_rot_rot(j, x, rut3, num_param, k, res2);

    std::cout << "Res2: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;





    //-------------------Test2---------------------
    std::cout << "Test2. K в дальней зоне" << std::endl;
    K_rot_rot_Far(j, x, rut1, num_param, k, res);
    std::cout << "Res1: " << res[0] << " " << res[1] << " " << res[2] << std::endl;

    K_rot_rot_Far(j, x, rut2, num_param, k, res1);
    K_rot_rot_Far(j, x, rut3, num_param, k, res2);

    std::cout << "Res2: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;



    //-------------------Test3---------------------
    std::cout << "Test3. K в ближней зоне с выделением особенности" << std::endl;
    std::complex<double> res3[3];
    K_rot_rot_Near(j, x, rut1, num_param, k, res3);
    std::cout << "Res1: " << res3[0] << " " << res3[1] << " " << res3[2] << std::endl;



    std::complex<double> res13[3], res23[3];
    K_rot_rot_Near(j, x, rut2, num_param, k, res13);
    K_rot_rot_Near(j, x, rut3, num_param, k, res23);

    std::cout << "Res2: " << res13[0] + res23[0] << " " << res13[1] + res23[1] << " " << res13[2] + res23[2] << std::endl;
    
    return 0;
}


