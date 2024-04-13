#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <string>
#include <vector>

#include "common_type.h"


//-----------------------Подсчет и запись токов, разложенных по базису---------------------------
void  get_j_basis(const TGrid_DC_Full& a, const std::complex<double>* b, const std::string filename,
    std::vector<std::complex<double>[3]> &j_vec)
{
    //Создание файлов для записи токов
    std::ofstream fout_j_real(filename + "_real.gv");
    if (!fout_j_real.is_open()) {
        std::cout << "Open " << filename << "_real.gv error" << std::endl;
        return;
    }

    std::ofstream fout_j_image(filename + "_image.gv");
    if (!fout_j_image.is_open()) {
        std::cout << "Open " << filename << "_image.gv error" << std::endl;
        return;
    }

    fout_j_real << 2 << " " << a.num_frm / 2 << std::endl;
    fout_j_image << 2 << " " << a.num_frm / 2 << std::endl;

    //Расчет и запись токов
    int i0;
    for (int i = 0; i < a.num_frm; i++)
    {
        i0 = 2 * i;
        j_vec[i][0] = b[i0] * a.cell_list[i].tau[0][0] + b[i0 + 1] * a.cell_list[i].tau[1][0];
        j_vec[i][1] = b[i0] * a.cell_list[i].tau[0][1] + b[i0 + 1] * a.cell_list[i].tau[1][1];
        j_vec[i][2] = b[i0] * a.cell_list[i].tau[0][2] + b[i0 + 1] * a.cell_list[i].tau[1][2];

        fout_j_real << j_vec[i][0].real() << " " << j_vec[i][1].real() << " " << j_vec[i][2].real() << std::endl;
        fout_j_image << j_vec[i][0].imag() << " " << j_vec[i][1].imag() << " " << j_vec[i][2].imag() << std::endl;
    }

    //Закрытие файлов
    fout_j_real.close();
    fout_j_image.close();
}





#include <iostream>
#include <iomanip>
#include <limits>


//============================================================================================================
//-------------------------Дополнительная функция считывания готовых токов------------------------------------
//============================================================================================================

void get_j_from_files(const std::string filename_real,
        const std::string filename_image,
        std::vector<std::complex<double>[3]> &j_vec)
{
    std::ifstream fin_real(filename_real);                // Файл с током, реальная часть
    if (!fin_real.is_open()) {
        std::cout << "Open " << filename_real << " error" << std::endl;
        return;
    }

    std::ifstream fin_image(filename_image);                // Файл с током, мнимая часть
    if (!fin_image.is_open()) {
        std::cout << "Open " << filename_image << " error" << std::endl;
        return;
    }


    const auto digits = std::numeric_limits<double>::digits10;
    fin_real >> std::fixed >> std::setprecision(digits);
    fin_image >> std::fixed >> std::setprecision(digits);

    int n1_real, n2_real, n1_image, n2_image, n;
    fin_real >> n1_real >> n2_real;
    fin_image >> n1_image >> n2_image;
    if (n1_real != n1_image || n2_real != n2_image)
    {
        return;
    }
    n = n1_real * n2_real;
    double x1_real, x2_real, x3_real, x1_image, x2_image, x3_image;
    
    for (int i = 0; i < n; i++)
    {
        fin_real >> x1_real >> x2_real >> x3_real;
        fin_image >> x1_image >> x2_image >> x3_image;
        j_vec[i][0] = std::complex<double>(x1_real, x1_image);
        j_vec[i][1] = std::complex<double>(x2_real, x2_image);
        j_vec[i][2] = std::complex<double>(x3_real, x3_image);
    }

    fin_real.close();
    fin_image.close();
}
