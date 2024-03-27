#ifndef _POLAR_TYPE_H_
#define _POLAR_TYPE_H_

class E0_Polar_Type
{
public:
    int type;//тип поляризации: H - 0, V - 1, своя - 2
    double e0[3];//вектор поляризации

    E0_Polar_Type(const int type_in, const double* k_in); //Горизонтальная(k_in = 0) / Вертикальная(k_in = 1) поляризации
    E0_Polar_Type(const double* e0_in); //Не типичная поляризация, в направлении некоторого вектора

};
#endif

