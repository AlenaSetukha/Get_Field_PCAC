#ifndef _WRITE_TO_FILE_H_
#define _WRITE_TO_FILE_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <vector>
#include <string>



#include "element_geom.h"
#include "ED_Par.h"

//===============================================================================================
//--------------------Запись результатов в файл в формате Aereco---------------------------------
//===============================================================================================
/**
 * Функции:
 *      write_j_toFiles - запись комплексных токов на ячейках в файлы "_real.gv", "_image.gv"
 *      write_j_toFile - запись вещественых токов на ячейках в файл ".gv"
 *      write_rkt_toFile - запись точек коллокаций в файл "body.gr"
 *      write_field_toFiles - запись косплексных значений поля на ячейках в файлы "_real.gv", "_image.gv", "_abs.gdr"
 *      write_ED_Params_toFile - запись э/д параметров задачи в файл "ed_param.res"
 */


void write_j_toFiles(const std::string& filename,
        const std::vector<std::complex<double>[3]> &j);

void write_j_toFile(const std::string& filename,
        const std::vector<double[3]> &j);


template<size_t CellPoints>
void write_rkt_toFile(const std::string& filename,
        const std::vector<double[CellPoints][3]> &cells)
{
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Open " + filename + " error: write_rkt_toFile");
    }

    int num_frm = cells.size();
    double cm[3];
    fout << 2 << " " << num_frm / 2 << std::endl;

    for (int i = 0; i < num_frm; i++)
    {
        get_center_mass(cells[i], cm);
        fout << cm[0] << " " << cm[1] << " " << cm[2] << std::endl;
    }
    fout.close();
}



void write_field_toFiles(const std::string& filename,
        const std::vector<std::complex<double>[3]> &field_val, 
        const int n1, const int n2);


void write_ED_Params_toFile(const std::string& filename_out,
        const ED_Par& ed_param);




#endif