#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "simulador.hpp"

namespace py = pybind11;
using namespace SimuladorAstrofisico;

// ========================================================================
// Módulo PyBind11
// ========================================================================
PYBIND11_MODULE(simulador_impacto, m) {
    m.doc() = "Módulo HPC de simulación de impactos de asteroides mediante Monte Carlo y Verlet."; // Docstring del módulo

    // --------------------------------------------------------------------
    // Exponiendo la Estructura de Datos SoA (Entrada)
    // --------------------------------------------------------------------
    py::class_<EstadoSistemaSoA>(m, "EstadoSistemaSoA", 
        "Estructura orientada a datos (SoA) para las condiciones iniciales del sistema.")
        .def(py::init<>(), "Constructor por defecto.")
        .def_readwrite("num_cuerpos", &EstadoSistemaSoA::num_cuerpos)
        // Gracias a <pybind11/stl.h>, Python interpretará estos std::vector como listas nativas de Python
        .def_readwrite("x", &EstadoSistemaSoA::x)
        .def_readwrite("y", &EstadoSistemaSoA::y)
        .def_readwrite("z", &EstadoSistemaSoA::z)
        .def_readwrite("vx", &EstadoSistemaSoA::vx)
        .def_readwrite("vy", &EstadoSistemaSoA::vy)
        .def_readwrite("vz", &EstadoSistemaSoA::vz)
        .def_readwrite("ax", &EstadoSistemaSoA::ax)
        .def_readwrite("ay", &EstadoSistemaSoA::ay)
        .def_readwrite("az", &EstadoSistemaSoA::az)
        .def_readwrite("masa", &EstadoSistemaSoA::masa)
        .def_readwrite("radio", &EstadoSistemaSoA::radio)
        .def_readwrite("nombres", &EstadoSistemaSoA::nombres);

    // --------------------------------------------------------------------
    // Exponiendo la Estructura de Resultados (Salida)
    // --------------------------------------------------------------------
    py::class_<ResultadosMonteCarlo>(m, "ResultadosMonteCarlo", 
        "Contenedor estadístico de los resultados de la simulación masiva.")
        // Usamos def_readonly para que en Python estos valores no puedan ser mutados accidentalmente
        .def_readonly("total_simulaciones", &ResultadosMonteCarlo::total_simulaciones)
        .def_readonly("impactos_tierra", &ResultadosMonteCarlo::impactos_tierra)
        .def_readonly("probabilidad_impacto", &ResultadosMonteCarlo::probabilidad_impacto)
        .def_readonly("desglose_colisiones", &ResultadosMonteCarlo::desglose_colisiones,
            "Diccionario con el conteo de choques por cuerpo celeste.");

    // --------------------------------------------------------------------
    // Exponiendo el Motor Principal
    // --------------------------------------------------------------------
    py::class_<MotorVerletMonteCarlo>(m, "MotorVerletMonteCarlo",
        "Motor numérico de alta precisión para N-Cuerpos.")
        
        // Constructor con py::arg para permitir named arguments (kwargs) en Python
        .def(py::init<EstadoSistemaSoA, double, double>(),
             py::arg("estado_inicial"), 
             py::arg("paso_tiempo"), 
             py::arg("t_max"),
             "Inicializa el motor. Valida la consistencia de los arreglos SoA.")
        
        // Método de ejecución masiva
        .def("ejecutar_monte_carlo", &MotorVerletMonteCarlo::ejecutarMonteCarlo,
             py::arg("num_corridas"), 
             py::arg("sigma_pos"), 
             py::arg("sigma_vel"),
             
             // ¡LA LÍNEA MÁS IMPORTANTE DEL BINDING!
             // Libera el Global Interpreter Lock (GIL) de Python antes de ejecutar C++.
             // Esto permite que OpenMP dispare todos los núcleos del CPU al 100% 
             // sin que la máquina virtual de Python los frene. Al terminar, recupera el GIL.
             py::call_guard<py::gil_scoped_release>(),
             
             "Ejecuta el bloque de simulaciones asíncronas de Monte Carlo.");
}