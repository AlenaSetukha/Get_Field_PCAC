cmake_minimum_required(VERSION 3.16)# Проверка версии CMake.
									# Если версия установленной программы
									# старее указаной, произойдёт аварийный выход.

set(CMAKE_CXX_STANDARD 11)          # стандарт компилятора

project(ED_OperatorLib VERSION 1.0)	# Название проекта


# Путь к заголовочным файлам
include_directories(include)

# Создание библиотеки
add_library(ED_OperatorLib STATIC
    src/element_geom.cpp
    src/f_par.cpp
    src/integral_par.cpp
    src/integrals_analytic.cpp
    src/K.cpp
    src/kernel_lib.cpp
    src/Num_Par.cpp
)

# Установите свойства для библиотеки
set_target_properties(ED_OperatorLib PROPERTIES
    VERSION 1.0
    SOVERSION 1
    PUBLIC_HEADER include/ED_OperatorLib/constants.h
    PUBLIC_HEADER include/ED_OperatorLib/element_geom.h
    PUBLIC_HEADER include/ED_OperatorLib/f_par.h
    PUBLIC_HEADER include/ED_OperatorLib/integral_par.h
    PUBLIC_HEADER include/ED_OperatorLib/integral_universal_pnt.h
    PUBLIC_HEADER include/ED_OperatorLib/integral_universal_seg_pnt.h
    PUBLIC_HEADER include/ED_OperatorLib/integrals_analytic.h
    PUBLIC_HEADER include/ED_OperatorLib/K.h
    PUBLIC_HEADER include/ED_OperatorLib/k0.h
    PUBLIC_HEADER include/ED_OperatorLib/kernel_lib.h
    PUBLIC_HEADER include/ED_OperatorLib/Num_Par.h
    PUBLIC_HEADER include/ED_OperatorLib/R.h
)


# Параметры установки
install(TARGETS ED_OperatorLib
    EXPORT ED_OperatorLibTargets
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/ED_OperatorLib
)



install(EXPORT ED_OperatorLibTargets
    FILE ED_OperatorLibTargets.cmake
    NAMESPACE ED_OperatorLib::
    DESTINATION lib/cmake/ED_OperatorLib
)













# Укажите дополнительную информацию (например, для CPack)
include(CPack)