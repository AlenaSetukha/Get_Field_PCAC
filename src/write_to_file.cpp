#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ED_Par.h"
#include "constants.h"
#include "element_geom.h"

void write_j_toFiles(const std::string &filename,
                     const std::vector<std::complex<double>[3]> &j) {
  // Создание файлов для записи токов
  std::ofstream fout_j_real(filename + "_real.gv");
  std::ofstream fout_j_image(filename + "_image.gv");

  if (!fout_j_real.is_open() || !fout_j_image.is_open()) {
    throw std::runtime_error("Open " + filename +
                             "_image.gv error: write_j_toFiles");
  }

  int num_frm = j.size();
  fout_j_real << 2 << " " << num_frm / 2 << std::endl;
  fout_j_image << 2 << " " << num_frm / 2 << std::endl;

  // Расчет и запись токов
  for (int i = 0; i < num_frm; i++) {
    fout_j_real << j[i][0].real() << " " << j[i][1].real() << " "
                << j[i][2].real() << std::endl;
    fout_j_image << j[i][0].imag() << " " << j[i][1].imag() << " "
                 << j[i][2].imag() << std::endl;
  }

  // Закрытие файлов
  fout_j_real.close();
  fout_j_image.close();
}

void write_j_toFile(const std::string &filename,
                    const std::vector<double[3]> &j) {
  std::ofstream fout_j(filename);
  if (!fout_j.is_open()) {
    throw std::runtime_error("Open " + filename + " error: write_j_toFile");
  }

  int num_frm = j.size();
  fout_j << 2 << " " << num_frm / 2 << std::endl;
  for (int i = 0; i < num_frm; i++) {
    fout_j << j[i][0] << " " << j[i][1] << " " << j[i][2] << std::endl;
  }
  fout_j.close();
}

void write_field_toFiles(const std::string &filename,
                         const std::vector<std::complex<double>[3]> &field_val,
                         const int n1, const int n2) {
  // Создание файлов с результатами
  std::ofstream fout_u_real(filename + "_real.gv");
  std::ofstream fout_u_image(filename + "_image.gv");
  std::ofstream fout_u_abs(filename + "_abs.gdr");

  if (!fout_u_abs.is_open() || !fout_u_image.is_open() ||
      !fout_u_real.is_open()) {
    throw std::runtime_error("Open " + filename +
                             " error: write_field_toFiles");
  }

  fout_u_real << n1 << " " << n2 << std::endl;
  fout_u_image << n1 << " " << n2 << std::endl;
  fout_u_abs << n1 << " " << n2 << std::endl;

  int n = n1 * n2; // количество точек, в которых рассчитано поле

  for (int i = 0; i < n; i++) {
    fout_u_real << field_val[i][0].real() << " " << field_val[i][1].real()
                << " " << field_val[i][2].real() << std::endl;
    fout_u_image << field_val[i][0].imag() << " " << field_val[i][1].imag()
                 << " " << field_val[i][2].imag() << std::endl;
    fout_u_abs << vec_length(field_val[i]) << std::endl;
  }

  // Закрытие файлов
  fout_u_real.close();
  fout_u_image.close();
  fout_u_abs.close();
}

void write_ED_Params_toFile(const std::string &filename_out,
                            const ED_Par &ed_param) {
  std::ofstream fout_params(filename_out);
  if (!fout_params.is_open()) {
    throw std::runtime_error("Open " + filename_out +
                             " error: write_ED_Params_toFile");
  }

  fout_params << "k_medium: " << ed_param.k_medium << std::endl;
  for (int i = 0; i < ed_param.k_medium; i++) {
    fout_params << "medium " << i << ":" << std::endl;
    fout_params << "    eps = " << ed_param.eps_d[i] * Constants::eps0
                << std::endl;
    fout_params << "    mu = " << ed_param.mu_d[i] * Constants::m0 << std::endl;
    fout_params << "    k = " << ed_param.k[i] << std::endl;
    fout_params << "    lambda = " << ed_param.lambda[i] << std::endl;
  }

  fout_params << "k_vec: " << ed_param.k_vec[0] << " " << ed_param.k_vec[1]
              << " " << ed_param.k_vec[2] << std::endl;
  fout_params << "e0_H: " << ed_param.e0_H[0] << " " << ed_param.e0_H[1] << " "
              << ed_param.e0_H[2] << std::endl;
  fout_params << "e0_V: " << ed_param.e0_V[0] << " " << ed_param.e0_V[1] << " "
              << ed_param.e0_V[2] << std::endl;
  fout_params << "phi = " << ed_param.phi / (M_PI / 180.) << std::endl;
  fout_params << "teta = " << ed_param.teta / (M_PI / 180.) << std::endl;
  fout_params << "omega0 = " << ed_param.omega0 << std::endl;
  fout_params << "lambda0 = " << ed_param.lambda0 << std::endl;
  fout_params.close();
}
