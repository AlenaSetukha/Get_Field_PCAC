#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <fstream>
#include <limits>
#include <iomanip>
#include <limits>
#include <numbers>



//-----------------------Подсчет токов, разложенных по базису---------------------------
void  get_j_basis(const std::vector<double[2][3]> &tau,
        const std::complex<double>* b,
        std::vector<std::complex<double>[3]> &j_vec)
{
    int i0, num_frm = tau.size();
    for (int i = 0; i < num_frm; i++)
    {
        i0 = 2 * i;
        j_vec[i][0] = b[i0] * tau[i][0][0] + b[i0 + 1] * tau[i][1][0];
        j_vec[i][1] = b[i0] * tau[i][0][1] + b[i0 + 1] * tau[i][1][1];
        j_vec[i][2] = b[i0] * tau[i][0][2] + b[i0 + 1] * tau[i][1][2];
    }
}






//-------------------------Дополнительная функция считывания готовых токов------------------------------------
void get_j_from_files(const std::string &filename_real,
        const std::string &filename_image, std::vector<std::complex<double>[3]> &j_vec)
{
    std::ifstream fin_real(filename_real);                // Файл с током, реальная часть
    std::ifstream fin_image(filename_image);                // Файл с током, мнимая часть

    if (!fin_real.is_open() || !fin_image.is_open())
    {
        throw std::runtime_error("Read j.txt error: get_j_from_files");
    }


    const auto digits = std::numeric_limits<double>::digits10;
    fin_real >> std::fixed >> std::setprecision(digits);
    fin_image >> std::fixed >> std::setprecision(digits);

    int n1_real, n2_real, n1_image, n2_image, n;
    fin_real >> n1_real >> n2_real;
    fin_image >> n1_image >> n2_image;

    if (n1_real != n1_image || n2_real != n2_image)
    {
        throw std::runtime_error("Read j.txt error: get_j_from_files");
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
