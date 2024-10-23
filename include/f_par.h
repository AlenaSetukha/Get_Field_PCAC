#ifndef _F_PAR_H_
#define _F_PAR_H_

#include <complex>
//===============================================================================================
//---------------------Класс параметров ядра интегрирования--------------------------------------
//===============================================================================================
/**
 * Поля:
 *      rs - радиус сглаживания подинтегральной функции
 *      k - волновое число функции источника(если есть)
 *      n - вектор внешней нормали
 *      e0 - направляющий вектор падающей волны
 *      a/ort - дополнительный вектор
 *      eps_edge - расстояние до края(для криволинейного интеграла)
 *      calc_dist - параметр аналитического/численного вычисления интеграла
 */

class f_par {
public:
  double rs, calc_dist, eps_edge;
  std::complex<double> k, vec_cmplx[3];
  double vec_dbl[3] = {0, 0, 0}, n[3], ort[3], e0[3];

  f_par();
  f_par(const double rs_in);
  f_par(const double rs_in, const std::complex<double> k_in);
  f_par(const double rs_in, const std::complex<double> k_in,
        const double *n_in);

  f_par(const f_par &obj);
};
#endif
