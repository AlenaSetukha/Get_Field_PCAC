#ifndef _KERNEL_PAR_H_
#define _KERNEL_PAR_H_

#include <complex>
//===============================================================================================
//---------------------Класс параметров ядра интегрирования--------------------------------------
//===============================================================================================
/**
 * Поля:
 *      rs - радиус сглаживания подинтегральной функции
 *      k - волновое число функции источника(если есть)
 *      calc_dist - параметр аналитического/численного вычисления интеграла(если есть)
 *      eps_edge - расстояние до края(для криволинейного интеграла)(если есть)
 *      vec_cmplx/vec_dbl - дополнительные векторы(если есть)
 *      n - вектор нормали(если есть)
 *      e0 - направляющий вектор падающей волны(если есть)
 *      ort - дополнительный вектор(если есть)
 */

class Kernel_Par {
public:
  double rs;
  std::complex<double> k;
  double calc_dist, eps_edge;

  std::complex<double> vec_cmplx[3];
  double vec_dbl[3], n[3], ort[3], e0[3];

  Kernel_Par() = default;
  Kernel_Par(const double rs_in);
  Kernel_Par(const double rs_in, const std::complex<double> k_in);
  Kernel_Par(const double rs_in, const std::complex<double> k_in,
        const double *n_in);

  Kernel_Par(const Kernel_Par& obj);
  ~Kernel_Par() = default;
};
#endif // _KERNEL_PAR_H_
