#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "ED_Par.h"
#include "constants.h"
#include "element_geom.h"

using namespace Constants;
//===============================================================================================
//---------------------------------Constructor1--------------------------------------------------
//===============================================================================================
// ЭД параметры сразу задаются реальными
ED_Par::ED_Par(const std::string &filename) {
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    throw std::runtime_error("Read ED_Param.txt error");
  }

  fin >> phi >> teta >> omega0;
  fin >> e0_H[0] >> e0_H[1] >> e0_H[2];
  fin >> e0_V[0] >> e0_V[1] >> e0_V[2];
  fin >> k_medium;
  phi *= M_PI / 180.;
  teta *= M_PI / 180.;
  lambda0 = 2. * M_PI * Constants::c_vacuum / omega0;

  mu_d = new std::complex<double>[k_medium];
  eps_d = new std::complex<double>[k_medium];
  k = new std::complex<double>[k_medium];
  lambda = new std::complex<double>[k_medium];

  for (int i = 0; i < k_medium; i++) {
    fin >> eps_d[i];
    fin >> mu_d[i];
    k[i] = omega0 * sqrt(eps_d[i] * mu_d[i] * Constants::eps0 * Constants::m0);
    lambda[i] = M_PI * 2. / k[i];
  }

  // Волновой вектор внешней среды
  k_vec[0] = -std::abs(k[0]) * cos(phi) * sin(teta);
  k_vec[1] = -std::abs(k[0]) * sin(phi) * sin(teta);
  k_vec[2] = -std::abs(k[0]) * cos(teta);

  fin >> epr_phiStart >> epr_phiEnd >> epr_step;
  epr_phiStart *= M_PI / 180.;
  epr_phiEnd *= M_PI / 180.;
  epr_step *= M_PI / 180.;
  fin.close();
}

//===============================================================================================
//---------------------------------Constructor2--------------------------------------------------
//===============================================================================================
/**
 * ЭД параметры для диэлектрика. Две среды - внещняя и внутренняя.
 * Задаются реальными, происходит пересчет к единичному радиусу.
 * Внешняя среда немагнитная
 * real_radius задается в мкм
 */
ED_Par::ED_Par(const double lambda0, const double real_radius,
               const double program_radius, const double ne, const double ni,
               const double ki, const double phi_in, const double teta_in,
               const double epr_phiStartIn, const double epr_phiEndIn) {
  phi = phi_in * M_PI / 180.;
  teta = teta_in * M_PI / 180.;
  epr_phiStart = epr_phiStartIn * M_PI / 180.;
  epr_phiEnd = epr_phiEndIn * M_PI / 180.;

  k_medium = 2;
  mu_d = new std::complex<double>[k_medium];
  eps_d = new std::complex<double>[k_medium];
  k = new std::complex<double>[k_medium];
  lambda = new std::complex<double>[k_medium];

  // 0 - внешняя среда, 1 - внутренняя(облучаемая частица)
  omega0 = 2. * M_PI * Constants::c_vacuum * 1e6 / lambda0; // реальная омега

  // Внешняя область
  mu_d[0] = std::complex<double>(1., 0.);
  eps_d[0] = std::complex<double>(ne * ne, 0.);
  k[0] = omega0 * sqrt(eps_d[0] * mu_d[0] * Constants::eps0 *
                       Constants::m0);               // реальное k
  k[0] = k[0] * real_radius * 1e-6 / program_radius; // программное k
  lambda[0] = lambda0 * 1e-6 / ne;                   // реальное lambda [м]
  lambda[0] =
      lambda[0] / (real_radius * 1e-6 / program_radius); // программное lambda

  std::cout << "Реальное k внешней среды: " << k[0] << std::endl;
  std::cout << "Программное k внешней среды: " << k[0] << std::endl;
  std::cout << "Реальное lambda внешней среды: " << lambda[0] << std::endl;
  std::cout << "Программное lambda внешней среды: " << lambda[0] << std::endl;

  // Внутренняя область
  mu_d[1] = std::complex<double>(1., 0.);
  eps_d[1] = std::complex<double>(ni * ni - ki * ki, 2. * ni * ki);
  k[1] = omega0 * sqrt(eps_d[1] * mu_d[1] * Constants::eps0 *
                       Constants::m0);               // реальное k
  k[1] = k[1] * real_radius * 1e-6 / program_radius; // программное k
  lambda[1] = lambda0 * 1e-6 / ni;                   // реальное lambda [м]
  lambda[1] =
      lambda[1] / (real_radius * 1e-6 / program_radius); // программное lambda

  omega0 = omega0 * real_radius * 1e-6 / program_radius; // программное омега

  // Волновой вектор внешней среды
  k_vec[0] = -std::abs(k[0]) * cos(phi) * sin(teta);
  k_vec[1] = -std::abs(k[0]) * sin(phi) * sin(teta);
  k_vec[2] = -std::abs(k[0]) * cos(teta);
}

//===============================================================================================
//--------------------------------Copy
// constructor-----------------------------------------------
//===============================================================================================
ED_Par::ED_Par(const ED_Par &ed_obj) {
  k_medium = ed_obj.k_medium;
  phi = ed_obj.phi;
  teta = ed_obj.teta;
  omega0 = ed_obj.omega0;
  lambda0 = ed_obj.lambda0;

  epr_phiStart = ed_obj.epr_phiStart;
  epr_phiEnd = ed_obj.epr_phiEnd;
  epr_step = ed_obj.epr_step;

  k_vec[0] = ed_obj.k_vec[0], k_vec[1] = ed_obj.k_vec[1],
  k_vec[2] = ed_obj.k_vec[2];
  e0_H[0] = ed_obj.e0_H[0], e0_H[1] = ed_obj.e0_H[1], e0_H[2] = ed_obj.e0_H[2];
  e0_V[0] = ed_obj.e0_V[0], e0_V[1] = ed_obj.e0_V[1], e0_H[2] = ed_obj.e0_V[2];

  mu_d = new std::complex<double>[k_medium];
  eps_d = new std::complex<double>[k_medium];
  k = new std::complex<double>[k_medium];
  lambda = new std::complex<double>[k_medium];

  for (int i = 0; i < k_medium; i++) {
    mu_d[i] = ed_obj.mu_d[i];
    eps_d[i] = ed_obj.eps_d[i];
    k[i] = ed_obj.k[i];
    lambda[i] = ed_obj.lambda[i];
  }
}

//===============================================================================================
//------------------------------------Destructor-------------------------------------------------
//===============================================================================================
ED_Par::~ED_Par() {
  delete[] eps_d;
  delete[] mu_d;
  delete[] k;
  delete[] lambda;
}
