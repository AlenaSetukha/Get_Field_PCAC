#ifndef _GET_J_H_
#define _GET_J_H_

#include <complex>
#include <string>
#include "common_type.h"


//===============================================================================================
//-----------Подсчет токов в точках коллокации и запись токов, разложенных по базису-------------
//===============================================================================================
// Параметры:
//      a - четырехугольная сетка TGrid_DC_Full
//      b - коэффициенты системы, ответ решения СЛАУ
//      filename - имя для записи токов
//      j_vec - токи, результат
void  get_j_basis(const TGrid_DC_Full& a, const std::complex<double>* b, const std::string filename,
        std::complex<double>**  j_vec);



//-------------------------Дополнительная функция считывания готовых токов-----------------------
void get_j_from_files(const std::string &filename_real, const std::string &filename_image, std::vector<std::complex<double>[3]> &j_vec);
#endif

