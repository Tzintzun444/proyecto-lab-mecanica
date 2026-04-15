#pragma once

#ifndef SIMULADOR_IMPACTO_HPP
#define SIMULADOR_IMPACTO_HPP

#include <vector>
#include <string>
#include <map>
#include <cstddef>

namespace SimuladorAstrofisico {

    // ========================================================================
    // Constantes Físicas Fundamentales
    // ========================================================================
    // Constante Gravitacional Gaussiana (k). 
    // Usando k^2 obtenemos G = 0.0002959122082855911 
    // Unidades: Masas Solares, Unidades Astronómicas (UA), y Días.
    constexpr double G_GAUSSIANA = 0.0002959122082855911;

    // ========================================================================
    // Estructuras de Datos Orientadas a Rendimiento (Data-Oriented Design)
    // ========================================================================

    struct alignas(64) EstadoSistemaSoA {
        std::size_t num_cuerpos;

        // Cinemática
        std::vector<double> x, y, z;
        std::vector<double> vx, vy, vz;
        
        // Dinámica
        std::vector<double> ax, ay, az;

        // Propiedades Físicas
        std::vector<double> masa;  // Masas Solares
        std::vector<double> radio; // Radios físicos en UA (para detección de colisión)

        // Metadatos
        std::vector<std::string> nombres; // Ej. "Sol", "Tierra", "Apophis"
    };

    struct ResultadosMonteCarlo {
        std::size_t total_simulaciones;
        std::size_t impactos_tierra;
        double probabilidad_impacto;
        std::map<std::string, std::size_t> desglose_colisiones; // Mapea nombre -> conteo de choques
    };

    // ========================================================================
    // Motor de Simulación Principal
    // ========================================================================

    class MotorVerletMonteCarlo {
    private:
        EstadoSistemaSoA estado_base_; // Estado "limpio" del cual parten todos los hilos
        double dt_;                    // Paso de tiempo de integración (en días)
        double tiempo_maximo_;         // Tiempo máximo a simular (en días)
        std::size_t indice_tierra_;    // Caché del índice de la Tierra para búsquedas rápidas
        std::size_t indice_asteroide_; // Caché del índice del asteroide

        // --------------------------------------------------------------------
        // Funciones Internas de Alto Rendimiento (Inlines)
        // --------------------------------------------------------------------
        
        void calcularAceleraciones(EstadoSistemaSoA& estado) const noexcept;

        [[nodiscard]] std::string ejecutarSimulacionAislada(EstadoSistemaSoA& estado_local) const noexcept;

    public:
        // --------------------------------------------------------------------
        // Constructor y Configuración
        // --------------------------------------------------------------------
        MotorVerletMonteCarlo(EstadoSistemaSoA estado_inicial, double paso_tiempo, double t_max);

        // --------------------------------------------------------------------
        // Interfaz de Ejecución Asíncrona
        // --------------------------------------------------------------------

        [[nodiscard]] ResultadosMonteCarlo ejecutarMonteCarlo(
            std::size_t num_corridas, 
            double sigma_pos, 
            double sigma_vel
        ) const;
    };

}

#endif // SIMULADOR_IMPACTO_HPP