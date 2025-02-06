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

class Integral_Par {
public:
    int idim;
    int n_start, p_max;
    double eps;

    Integral_Par() = default;
    Integral_Par(const int idim_in, const int n_start_in, const int p_max_in, const double eps_in);
    Integral_Par(const Integral_Par& obj);
    ~Integral_Par() = default;
};
#endif // _INTEGRAL_PAR_H_

