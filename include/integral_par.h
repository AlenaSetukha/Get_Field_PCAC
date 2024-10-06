#ifndef _INTEGRAL_PAR_H_
#define _INTEGRAL_PAR_H_

#include <iostream>
//===============================================================================================
//-------------------------Класс параметров интегрирования---------------------------------------
//===============================================================================================
/**
 * Поля:
 *      idim - размерность подынтегральной функции
 *      n_start - стартовое разбиение ячейки/отрезка 
 *      p_max - предельное число разбиений(итераций)- 2^p_max
 *      eps - точность вычисления интеграла(имеет смысл при p_max != 1)
 */

class integral_par
{
public:
    int idim, n_start, p_max;
    double eps;

    integral_par();
    integral_par(const int idim_in, const int n_start_in, const int p_max_in, const double eps_in);
    integral_par(const integral_par& obj);
};
#endif

