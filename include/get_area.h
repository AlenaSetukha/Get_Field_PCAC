#ifndef _GET_AREA_H_
#define _GET_AREA_H_

#include <iostream>
#include <vector>

#include "element_geom.h"

//===============================================================================================
//-------------------Определение области, где находится точка------------------------------------
//===============================================================================================
/**
 * Определение, где находится точка относительно заданного тела
 * Параметры:
 *      x[3] - точка
 *      cells[num_frm][4/3][3] - координаты ячеек тела
 *      norm[num_frm][3] - координаты нормалей
 * Результат: 0 - внутри объекта, -1 - снаружи.
 */

template<size_t CellPoints>
int get_area(const double* x,
        const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm) {
    double teta, vec[3];
    double sum = 0.;

    int num_frm = cells.size();
    for (int i = 0; i < num_frm; i++) {
        vec[0] = cells[i][0][0] - x[0];
        vec[1] = cells[i][0][1] - x[1];
        vec[2] = cells[i][0][2] - x[2];

        //123
        teta = solid_angle(cells[i][0], cells[i][1], cells[i][2], x);

        if (scal_prod(vec, norm[i]) > 0) {
            sum -= fabs(teta);
        } else {
            sum += fabs(teta);
        }

        //134
        teta = solid_angle(cells[i][0], cells[i][2], cells[i][3], x);

        if (scal_prod(vec, norm[i]) > 0) {
            sum -= fabs(teta);
        } else {
            sum += fabs(teta);
        }

    }

    //Result
    if (sum > -0.5) {
        return 0;
    } else return -1;
}
#endif // _GET_AREA_H_