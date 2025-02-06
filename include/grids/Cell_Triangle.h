#ifndef _CELL_TRIANGLE_H_
#define _CELL_TRIANGLE_H_

//===============================================================================================
//----------------------------Класс треугольной ячейки-------------------------------------------
//===============================================================================================
/**
 * Поля:
 *      cell - координаты 3-х вершин
 *      tau - координаты двух локальных базисных векторов
 *      rkt - координаты точки коллокации
 *      norm - координаты вектора нормали
 *      s - площадь ячейки
 * Методы:
 *      cell_fill - заполнение одной ячейки
 */

class Cell_Triangle {
public:
    double cell[3][3];
    double tau[2][3];
    double rkt[3];
    double norm[3];
    double s;

    Cell_Triangle() = default;
    Cell_Triangle(const Cell_Triangle& obj);
    ~Cell_Triangle() = default;

    void cell_fill(const double (&cell_in)[3][3]);
};

#endif // _CELL_TRIANGLE_H_