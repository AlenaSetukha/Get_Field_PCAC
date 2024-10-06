#ifndef _NUM_PAR_H_
#define _NUM_PAR_H_

#include <string>
//===============================================================================================
//-----------------------Класс численных параметров задачи---------------------------------------
//===============================================================================================
/**
 * Поля:
 *      eps - точность вычисления интегралов
 *      rs - радиус сглаживания для интеграла по ячейке
 *      rs_seg - радиус сглаживания для интеграла по отрезку
 *      n_start_seg - стартовое разбиение на отрезке
 *      n_start - стартовое разбиение на ячейке
 *      p_max - предельное разбиение на ячейке
 *      p_max_seg - предельное разбиение на отрезке
 *      k - число шагов по времени
 *      T - временной отрезок
 *      dt - шаг по времени
 *      kappa - параметр задачи(дополнительный)
 *      M - параметр задачи(дополнительный)
 * Радиус сглаживания может подаваться в долях от каждой ячейки/шага сетки,
 * в зависимости от применимости.
 */

class Num_Par
{
public:
    double eps, rs, rs_seg;
    int n_start_seg, n_start, p_max, p_max_seg, k;
    double T, dt;
    double kappa, M;

    Num_Par();
    Num_Par(const std::string &filename);
    Num_Par(const Num_Par& obj);
};
#endif
