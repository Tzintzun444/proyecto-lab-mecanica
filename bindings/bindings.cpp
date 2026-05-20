#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Cabecera esencial: Permite la traducción mágica de std::vector de C++ a Listas de Python.
#include "simulador.hpp"

// Alias de espacio de nombres para mantener el código de envoltura limpio
namespace py = pybind11;
using namespace SimuladorAstrofisico;

// ========================================================================
// DEFINICIÓN DEL MÓDULO PYBIND11 (Exportación a Python)
// ========================================================================
// El nombre "simulador_impacto" debe ser exactamente igual al definido en el CMakeLists.
PYBIND11_MODULE(simulador_impacto, m) {
    m.doc() = "Módulo HPC para dinámica orbital N-Cuerpos. Proyecto de Laboratorio de Mecánica."; 

    // --------------------------------------------------------------------
    // Exportación de la Estructura de Entrada (EstadoSistemaSoA)
    // --------------------------------------------------------------------
    py::class_<EstadoSistemaSoA>(m, "EstadoSistemaSoA", 
        "Carga los parámetros físicos (Masa, Radio, Vectores Cinéticos) hacia la memoria RAM (C++).")
        .def(py::init<>(), "Inicializa una estructura vacía en la memoria.")
        
        // Mapeamos los atributos de C++ para que Python pueda escribir sobre ellos directamente
        .def_readwrite("num_cuerpos", &EstadoSistemaSoA::num_cuerpos)
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
    // Exportación de la Estructura de Salida (ResultadosMonteCarlo)
    // --------------------------------------------------------------------
    py::class_<ResultadosMonteCarlo>(m, "ResultadosMonteCarlo", 
        "Contiene el desglose probabilístico retornado tras concluir el procesamiento.")
        
        // Usamos def_readonly para proteger la integridad estadística. 
        // Evita que un usuario en Python altere los resultados accidentalmente.
        .def_readonly("total_simulaciones", &ResultadosMonteCarlo::total_simulaciones)
        .def_readonly("impactos_tierra", &ResultadosMonteCarlo::impactos_tierra)
        .def_readonly("probabilidad_impacto", &ResultadosMonteCarlo::probabilidad_impacto)
        .def_readonly("desglose_colisiones", &ResultadosMonteCarlo::desglose_colisiones);

    // --------------------------------------------------------------------
    // Exportación de la Clase Principal (El Integrador)
    // --------------------------------------------------------------------
    py::class_<MotorVerletMonteCarlo>(m, "MotorVerletMonteCarlo",
        "Motor de integración Velocity Verlet acoplado a estadística de Monte Carlo.")
        
        // Exponemos el constructor con argumentos etiquetados (kwargs) para mejor legibilidad
        .def(py::init<EstadoSistemaSoA, double, double>(),
             py::arg("estado_inicial"), 
             py::arg("paso_tiempo"), 
             py::arg("t_max"))
        
        // Vinculamos la función pesada del proyecto
        .def("ejecutar_monte_carlo", &MotorVerletMonteCarlo::ejecutarMonteCarlo,
             py::arg("num_corridas"), 
             py::arg("sigma_pos"), 
             py::arg("sigma_vel"),
             
             // ** OPTIMIZACIÓN MULTITHREADING EXTREMA **
             // Liberamos el GIL (Global Interpreter Lock) antes de entrar a las matemáticas.
             // Sin esta directiva, Python frenaría los hilos de OpenMP dejándolos secuenciales.
             py::call_guard<py::gil_scoped_release>(),
             
             "Dispara la evaluación asíncrona de las órbitas usando los núcleos del CPU.");
}