#ifndef _COMMON_TYPE_H_
#define _COMMON_TYPE_H_

#include <string>
#include <vector>
#include <array>

#include "cell.h"

//===============================================================================================
//--------------------------Класс четырехугольной сетки------------------------------------------
//===============================================================================================
// Все ячейки различны

class TGrid_DC_Full
{
public:

    std::vector<Cell> cell_list;                 // [frm] - список всех ячеек 

    int num_obj;                                 // число объектов
    int num_mod;                                 // число модулей на всех объектах суммарно
    int num_frm;                                 // число ячеек

    std::vector<std::array<int, 2>> beg_end;     // [mod][2] номер 1 и посл ячейки на модуле с номером mod
    std::vector<std::array<int, 2>> beg_end_obj; // [obj][2] 1 и посл модуль на каждом объекте с нoмером obj

    double max_diag;                             // шаг сетки

    int num_bl;                                  // число линий отрыва
    int num_seg;                                 // число отрезков на всех линиях отрыва
    
    std::vector<std::array<std::array<double, 3>, 2>> x_bound;// [num_seg][2][3] - точки линий отрыва




    void check_step(const double (&root_tmp)[4][3], double max_diag_tmp);//поиск шага сетки

//------------------------------------------Constructor------------------------------------------
    TGrid_DC_Full(const std::string filename);//сквозная сетка на все объекты и все модули вместе

//----------------------------------------Copy Constructor---------------------------------------
    TGrid_DC_Full(const TGrid_DC_Full& obj);

//------------------------------------------Destructor-------------------------------------------
    ~TGrid_DC_Full();
};
#endif