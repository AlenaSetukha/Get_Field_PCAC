#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "integral_universal_pnt.h"
#include "integral_universal_seg_pnt.h"
#include "Kernel_Par.h"
#include "Integral_Par.h"
#include "element_geom.h"
#include "kernel_lib.h"


/**
 * Тест для проверки вычисления основных интеграторов по 
 * четырехугольным и треугольным ячейкам.
 * Порядок обхода ячеек важен: либо по, либо против часовой
 * стрелки - все одновременно.
 */

void f0(const double *x, const double *y,
                const Kernel_Par &param, double *ff) 
{
    ff[0] = x[0] * x[1] + 15.9;
    ff[1] = x[2] + x[1];
    ff[2] = x[0];
}




int main(int argc, char **argv)
{
    //-------Тестирование интеграторов-------------
    const double rut1[4][3] = {{0., 0., 0.},
                         {2., 0., 0.},
                         {3., 2., 0.},
                         {0., 2., 0.}};

    const double rut2[3][3] = {{0., 2., 0.},
                         {0., 0., 0.},
                         {2., 0., 0.}};
    
    
    const double rut3[3][3] = {{3., 2., 0.},
                         {0., 2., 0.},
                         {2., 0., 0.}};


    const double x[3] = {10., 10., 10.};


    //-----Тестирование поверхностного интеграла----
    //-------------------Test1---------------------+
    std::cout << "Test1" << std::endl;
    Kernel_Par param(0.5, std::complex<double>(1., 0.));
    Integral_Par int_param(1, 8, 3, 0.000001);

    std::complex<double> res31_c[1];
    integral_universal_pnt(x, rut1, f_simple_pot_G, param, int_param, res31_c);
    std::cout << "     Интеграл [4][3]: " << res31_c[0] << std::endl;
    

    std::complex<double> res32_c[1], res33_c[1];
    integral_universal_pnt(x, rut2, f_simple_pot_G, param, int_param, res32_c);
    integral_universal_pnt(x, rut3, f_simple_pot_G, param, int_param, res33_c);
    std::cout << "     Интеграл [3][3] + [3][3]: " << res32_c[0] + res33_c[0] << std::endl;

    //-------------------Test2---------------------
    std::cout << "Test2" << std::endl;
    Kernel_Par param2(0.01, std::complex<double>(1., 0.));
    Integral_Par int_param2(3, 4, 1, 0.000001);


    double res31_c2[3];
    integral_universal_pnt(x, rut1, f0, param2, int_param2, res31_c2);
    std::cout << "     Интеграл [4][3]: " << res31_c2[0] << " " << res31_c2[1] << " " << res31_c2[2] << std::endl;
    

    double res32_c2[3], res33_c2[3];
    integral_universal_pnt(x, rut2, f0, param2, int_param2, res32_c2);
    integral_universal_pnt(x, rut3, f0, param2, int_param2, res33_c2);
    std::cout << "     Интеграл [3][3] + [3][3]: " << res32_c2[0] + res33_c2[0] << " " << 
                    res32_c2[1] + res33_c2[1] << " " << res32_c2[2] + res33_c2[2] << std::endl;





    //---Тестирование криволинейного интеграла------
    //-------------------Test3---------------------+
    std::cout << "Test3" << std::endl;
    Kernel_Par param3(0.01, std::complex<double>(1., 0.));
    Integral_Par int_param3(1, 8, 1, 0.000001);
    double res_seg[1];
    integral_universal_seg_pnt(rut1[0], rut1[1], x, f0, param3, int_param3, res_seg);
    std::cout << "Интеграл по отрезку: " << res_seg[0] << std::endl;
    integral_universal_seg_pnt(rut1[1], rut1[0], x, f0, param3, int_param3, res_seg);
    std::cout << "Интеграл по отрезку: " << res_seg[0] << std::endl;
    
    return 0;
}