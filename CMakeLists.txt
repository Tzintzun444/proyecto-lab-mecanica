cmake_minimum_required(VERSION 3.15)
project(SimuladorAstrofisico VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Esto asegura que el binario de Python (.so o .pyd) se exporte directo a la carpeta raíz
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# 1. Encontrar dependencias críticas
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

find_package(OpenMP REQUIRED)
find_package(pybind11 REQUIRED)

# 2. Crear el módulo de Python 
pybind11_add_module(simulador_impacto 
    bindings/bindings.cpp 
    src/simulador.cpp
)

# Incluir los headers locales
target_include_directories(simulador_impacto PRIVATE 
    src 
)

# Linkear OpenMP
target_link_libraries(simulador_impacto PRIVATE OpenMP::OpenMP_CXX)

# 3. Optimización Interprocedural (Link-Time Optimization)
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET simulador_impacto PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ==========================================
# 4. Optimizaciones de Alto Rendimiento (Multiplataforma)
# ==========================================
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(simulador_impacto PRIVATE 
        -O3 
        -march=native     # Le dice al compilador que use las instrucciones del CPU actual (AVX/SIMD)
        -ffast-math       
        -funroll-loops    
        -Wall -Wextra     # Advertencias profesionales de código limpio
    )
    
    if(MINGW)
        target_link_options(simulador_impacto PRIVATE -static -static-libgcc -static-libstdc++)
    endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(simulador_impacto PRIVATE 
        /O2 
        /fp:fast 
        /arch:AVX2        # Fuerza instrucciones vectoriales en Windows
        /W3
        /openmp:experimental
    )
endif()