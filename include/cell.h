#ifndef _CELL_H_
#define _CELL_H_

/**
 * Класс четырехугольной ячейки
 * Поля:
 *      cell - координаты 4-х вершин
 *      tau - координаты локального базиса
 *      rkt - координаты точки коллокации
 *      norm - координаты нормали
 *      s - площадь ячейки
 * Методы:
 *      Cell_fill - заполнение одной ячейки
 */

class Cell
{
public:
    double cell[4][3];
    double tau[2][3];
    double rkt[3];
    double norm[3];
    double s;


    Cell();
    Cell(const Cell& obj);
    ~Cell(){};

    void Cell_fill(const double (&cell_in)[4][3]);
};
#endif