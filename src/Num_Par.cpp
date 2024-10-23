#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <string>

#include "Num_Par.h"

Num_Par::Num_Par(const std::string &filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Read num_param.txt error" << std::endl;
        exit(1);
    }
    fin >> eps >> n_start >> n_start_seg >> p_max >> p_max_seg >> rs >> rs_seg >> T >> k >> kappa >> M;
    dt = T / k;
    fin.close();
}

Num_Par::Num_Par(const double eps_in, const double rs_in, const double rs_seg_in,
            const int n_start_in, const int n_start_seg_in, const int p_max_in,
            const int p_max_seg_in)
{
    eps = eps_in, rs = rs_in, rs_seg = rs_seg_in;
    n_start = n_start_in, n_start_seg = n_start_seg_in;
    p_max = p_max_in, p_max_seg = p_max_seg_in;
    T = 0., dt = 0., kappa = 0., M = 0.;
    k = 0;
}

Num_Par::Num_Par(const Num_Par& obj)
{
    eps = obj.eps, rs = obj.rs, rs_seg = obj.rs_seg;
    n_start = obj.n_start, n_start_seg = obj.n_start_seg;
    p_max = obj.p_max, p_max_seg = obj.p_max_seg;
    T = obj.T, dt = obj.dt, kappa = obj.kappa, M = obj.M;
    k = obj.k;
}


