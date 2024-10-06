#ifndef _ED_PAR_H_
#define _ED_PAR_H_

#include <complex>
#include <string>

//=====================================================================================
//---------------------Класс электродинамических параметров----------------------------
//=====================================================================================
/**
 * Поля:
 *      k_medium - число областей, нумерация начинается с внешней области.
 *      eps_d[k_medium] - относительные диэлектрические проницаемости каждой области
 *      mu_d[k_medium] - относительные магнитные проницаемости каждой области
 *      k[k_medium] - волновые числа каждой области
 *      lambda[k_medium] - длина падающей волны в каждой области
 * 
 *      k_vec - волновой вектор внешнего поля в вакууме
 *      phi/teta - углы падающей волны: "x" -> "y", "z" -> "y". Задаются в градусах
 *      omega0 - круговая частота падающей волны в вакууме
 *      lambda0 - длины падающей волны в вакууме [мкм]
 *      epr_phiStart/End - угол в плоскости "oxy" начала/конца расчета ЭПР/Обратной ЭПР 
 *      epr_step - шаг расчета ЭПР
 * 
 * 
 * Все углы/шаги хранятся в радианах.
 * Если величина действительная, в файле задавать как (10., 0.)
 */
class ED_Par
{
public:
    int k_medium;
    std::complex<double> *eps_d, *mu_d, *k, *lambda;

    double k_vec[3], e0_H[3], e0_V[3];
    double phi, teta, omega0, lambda0;
    double epr_phiStart, epr_phiEnd, epr_step;

    ED_Par();
    ED_Par(const std::string& filename);
    ED_Par(const double lambda0,
        const double real_radius, const double program_radius,
        const double ne, const double ni, const double ki,
        const double phi_in, const double teta_in,
        const double epr_phiStartIn, const double epr_phiEndIn);
    ED_Par(const ED_Par& ed_obj);
    ~ED_Par();
};
#endif
