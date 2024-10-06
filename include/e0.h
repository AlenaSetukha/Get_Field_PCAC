#ifndef _POLAR_TYPE_H_
#define _POLAR_TYPE_H_


//=====================================================================================
//-----------------------------Поляризация падающей волны------------------------------
//=====================================================================================
/**
 * Поля:
 *      type: {0, 1, 2} - горизонтальная/вертикальная/произвольная
 *      e0 - орт вектора поляризации
 */

class E0_PolarType {
public:
    int type;
    double e0[3];

    E0_PolarType();
    E0_PolarType(const double* e0_in);
    E0_PolarType(const int type_in, const double* k_vec);
    E0_PolarType(const E0_PolarType& obj);
};
#endif

