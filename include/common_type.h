#ifndef _COMMON_TYPE_H_
#define _COMMON_TYPE_H_

#include <string>
#include <vector>
#include <array>

#include "cell.h"

//=====================================================================================
//--------------------------Класс четырехугольной сетки--------------------------------
//=====================================================================================
/**
 * Класс четырехугольной сетки с различными ячейками. Сквозная сетка
 * Поля:
 *      cell_list[frm] - сквозной список всех ячеек
 *      num_obj/num_mod/num_frm - число объектов/модулей/ячеек
 *      num_bl/num_blSeg - число линий отрыва/отрезков на всех линиях отрыва 
 *      grid_step - шаг сетки
 *      x_bound[num_blSeg][2][3] - координаты концов отрезков линий отрыва
 *      beg_endMod[mod][2] - номер 1 и посл ячейки на модуле
 *      beg_endObj[obj][2] - номер 1 и посл модуля на каждом объекте
 * Методы:
 *      check_step - вычисление шага сетки
 */

class TGrid_DC_Full
{
public:
    int num_obj, num_mod, num_frm;
    int num_bl, num_blSeg;
    double grid_step;

    std::vector<Cell> cell_list;

    std::vector<std::array<int, 2>> beg_endMod, beg_endObj;
    std::vector<std::array<std::array<double, 3>, 2>> x_bound;



    void check_step(const double (&cur_root)[4][3], double cur_grid_step);


    TGrid_DC_Full();

    TGrid_DC_Full(const std::string& filename);

    TGrid_DC_Full(const TGrid_DC_Full& obj);

    ~TGrid_DC_Full();
};
#endif