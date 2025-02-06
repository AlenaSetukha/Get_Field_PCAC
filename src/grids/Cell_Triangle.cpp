#include "Cell_Triangle.h"
#include "element_geom.h"


//===============================================================================================
//--------------------------------Copy constructor-----------------------------------------------
//===============================================================================================
Cell_Triangle::Cell_Triangle(const Cell_Triangle& obj) {
    for (int i = 0; i < 3; i++) {
        cell[i][0] = obj.cell[i][0], cell[i][1] = obj.cell[i][1], cell[i][2] = obj.cell[i][2]; 
    }
    for (int i = 0; i < 3; i++) {
        tau[0][i] = obj.tau[0][i], tau[1][i] =  obj.tau[1][i];
        rkt[i] = obj.rkt[i];
        norm[i] = obj.norm[i];
    }
    s = obj.s;
}




//===============================================================================================
//--------------------------------Функция заполнения---------------------------------------------
//===============================================================================================
void Cell_Triangle::cell_fill(const double (&cell_in)[3][3]) {
    // cell [3][3]
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cell[i][j] = cell_in[i][j];
        }
    }

    // rkt
    for (int i = 0; i < 3; i++) {
        rkt[i] = (cell[0][i] + cell[1][i] + cell[2][i]) / 3.;
    }

    // norm
    double ac[3], ab[3], len;

    for (int k = 0; k < 3; k++) {
        ab[k] = cell[1][k] - cell[0][k];
        ac[k] = cell[2][k] - cell[0][k];
    }

    vec_prod(ab, ac, norm);
    len = vec_length(norm);
    
    for (int k = 0; k < 3; k++) {
        norm[k] /= len;
    }

    // s
    s = len / 2;

    // tau[0] и tau[1]
    len = vec_length(ab);
    for (int k = 0; k < 3; k++) {
        tau[0][k] = ab[k] / len;
    }

    vec_prod(norm, tau[0], tau[1]);
    len = vec_length(tau[1]);
    for (int k = 0; k < 3; k++) {
        tau[1][k] /= len;
    }
}
