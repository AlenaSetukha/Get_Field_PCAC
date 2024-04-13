#ifndef _NUM_PAR_H_
#define _NUM_PAR_H_

class Num_Par
{
public:
    double eps, rs;// точность расчета интегралов, радиус сглаживания
    int n_start_seg, n_start, p_max, p_max_seg, k;//стартовое разбиение на отрезке, на ячейке, предельное разбиение на ячейке, отрезке, число шагов по времени
    double T, dt; //Временной отрезок, шаг по времени 
    double kappa, M;

    Num_Par(const std::string filename);
};
#endif
