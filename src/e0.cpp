#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "e0.h"
#include "element_geom.h"

E0_PolarType::E0_PolarType(const int type_in, const double* k_in)
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


E0_PolarType::E0_PolarType(const double* e0_in)
{
    e0[0] = e0_in[0], e0[1] = e0_in[1], e0[2] = e0_in[2]; 
    type = 2;
}

E0_PolarType::E0_PolarType()
{
    type = 0;
    e0[0] = 0., e0[1] = 0., e0[2] = 0.;
}

E0_PolarType::E0_PolarType(const E0_PolarType& obj)
{
    type = obj.type;
    e0[0] = obj.e0[0], e0[1] = obj.e0[1], e0[2] = obj.e0[2];
}



