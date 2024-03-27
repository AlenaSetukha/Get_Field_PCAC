#ifndef _ED_PAR_H_
#define _ED_PAR_H_

#include <complex>
#include <string>

//===============================================================================================
//---------------------Класс электродинамических параметров--------------------------------------
//===============================================================================================
class ED_Par
{
public:
    int k_medium;                // число областей
    std::complex<double>* eps_d; // если eps действ, то просто задать (10., 0.) в файле
    std::complex<double>* m_d;
    std::complex<double>* k;     // k[i] = omega * sqrt(eps[i] * mu[i] * mu0 * eps0)
    std::complex<double>* lambda;

    double k_vec[3];             // волновой вектор
    double alpha, omega, alpha_start, alpha_end; //угол падения, omega, угол начала/конца расчета ЭПР

    ED_Par(const std::string filename);
    ED_Par(const ED_Par& ed_obj);
    ~ED_Par();

};
#endif
