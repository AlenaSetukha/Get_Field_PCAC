#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "e0.h"
#include "element_geom.h"

E0_Polar_Type::E0_Polar_Type(const int type_in, const double* k_in)
{
    double ort[3], vv[3];
    ort[0] = 0.;
    ort[1] = 0.;
    ort[2] = 1.;
    if (!type_in)//H
    {
        vec_prod(k_in, ort, vv);
        double k_len = vec_length(vv);
        e0[0] = vv[0] / k_len;
        e0[1] = vv[1] / k_len;
        e0[2] = vv[2] / k_len;
    } else {
        e0[0] = 0.;
        e0[1] = 0.;
        e0[2] = 1.;
    }
}


E0_Polar_Type::E0_Polar_Type(const double* e0_in)
{
    e0[0] = e0_in[0], e0[1] = e0_in[1], e0[2] = e0_in[2]; 
    type = 2;
}

