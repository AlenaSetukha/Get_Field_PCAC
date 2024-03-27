#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <string>

#include "ED_Par.h"
#include "constants.h"
#include "element_geom.h"

//===============================================================================================
//--------------------------------Filling ED_Par-------------------------------------------------
//===============================================================================================
ED_Par::ED_Par(const std::string filename)
{
    Constants c;
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Read " << filename <<  " error" << std::endl;
        exit(1);
    }
    fin >> alpha >> omega >> alpha_start >> alpha_end;

    fin >> k_medium;

    m_d = new std::complex<double>[k_medium];
    eps_d = new std::complex<double>[k_medium];
    k = new std::complex<double>[k_medium];
    lambda = new std::complex<double>[k_medium];

    for (int i = 0; i < k_medium; i++)
    {
        fin >> eps_d[i];
        fin >> m_d[i];
        k[i] = omega * sqrt(eps_d[i] * m_d[i] * c.eps0 * c.m0);
        lambda[i] = c.pi * 2. / k[i];
    }

    k_vec[0] =  - abs_tmp(k[0]) * cos(alpha);
    k_vec[1] =  - abs_tmp(k[0]) * sin(alpha);
    k_vec[2] = 0.;
    fin.close();
    return;
}


//===============================================================================================
//--------------------------------Copy constructor-----------------------------------------------
//===============================================================================================
ED_Par::ED_Par(const ED_Par& ed_obj)
{
    k_medium = ed_obj.k_medium;
    omega = ed_obj.omega;
    alpha = ed_obj.alpha;
    alpha_start = ed_obj.alpha_start;
    alpha_end = ed_obj.alpha_end;

    k_vec[0] = ed_obj.k_vec[0], k_vec[1] = ed_obj.k_vec[1], k_vec[2] = ed_obj.k_vec[2]; 

    m_d = new std::complex<double>[k_medium];
    eps_d = new std::complex<double>[k_medium];
    k = new std::complex<double>[k_medium];
    lambda = new std::complex<double>[k_medium];

    for (int i = 0; i < k_medium; i++)
    {
        m_d[i] = ed_obj.m_d[i];
        eps_d[i] = ed_obj.eps_d[i];
        k[i] = ed_obj.k[i];
        lambda[i] = ed_obj.lambda[i];
    }

    return;
}


//===============================================================================================
//------------------------------------Destructor-------------------------------------------------
//===============================================================================================
ED_Par::~ED_Par()
{
    delete[] eps_d;
    delete[] m_d;
    delete[] k;
    delete[] lambda;
    return;
}

