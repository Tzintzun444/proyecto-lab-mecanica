#pragma once

#ifndef SIMULADOR_IMPACTO_HPP
#define SIMULADOR_IMPACTO_HPP

#include <vector>
#include <string>
#include <map>
#include <cstddef> 

namespace SimuladorAstrofisico {

    // ========================================================================
    // CONSTANTES FÍSICAS
    // ========================================================================
    // Se utiliza la Constante Gravitacional Gaussiana (k) al cuadrado para
    // mantener la precisión numérica. Si usáramos la G normal (6.67e-11), 
    // los exponentes tan pequeños causarían errores de truncamiento en C++.
    // Unidades del sistema: Masas Solares, Unidades Astronómicas (UA) y Días.
    constexpr double G_GAUSSIANA = 0.0002959122082855911;

    // ========================================================================
    // ESTRUCTURAS DE DATOS 
    // ========================================================================

    /**
     * @brief Estructura de Arreglos (SoA - Structure of Arrays).
     * Decidimos usar SoA en lugar de objetos individuales (AoS) para que los 
     * datos de posición y velocidad estén contiguos en la memoria RAM.
     * La directiva alignas(64) asegura que el bloque de memoria cuadre con las 
     * líneas de caché del procesador, permitiendo vectorización SIMD eficiente.
     */
    struct alignas(64) EstadoSistemaSoA {
        std::size_t num_cuerpos; // Cantidad total de cuerpos celestes en la simulación

        // Vectores cinemáticos (Posición en 3D)
        std::vector<double> x, y, z;
        
        // Vectores cinemáticos (Velocidad en 3D)
        std::vector<double> vx, vy, vz;
        
        // Vectores dinámicos (Aceleración en 3D calculada en cada iteración)
        std::vector<double> ax, ay, az;

        // Propiedades físicas de los cuerpos
        std::vector<double> masa;  // Masas relativas al Sol
        std::vector<double> radio; // Radios físicos (importante para detectar si chocan)

        // Identificadores de los cuerpos para el registro de colisiones
        std::vector<std::string> nombres; 
    };

    /**
     * @brief Estructura para almacenar y retornar los resultados del Monte Carlo.
     * Esta estructura será empaquetada por PyBind11 para ser leída desde Python.
     */
    struct ResultadosMonteCarlo {
        std::size_t total_simulaciones;       // Total de universos generados
        std::size_t impactos_tierra;          // Conteo específico de choques en la Tierra
        double probabilidad_impacto;          // Relación de impactos / total_simulaciones
        std::map<std::string, std::size_t> desglose_colisiones; // Diccionario con todos los choques registrados
    };

    // ========================================================================
    // CLASE PRINCIPAL: MOTOR DEL SIMULADOR
    // ========================================================================

    /**
     * @class MotorVerletMonteCarlo
     * @brief Implementa la integración numérica (Velocity Verlet) y paralelización.
     */
    class MotorVerletMonteCarlo {
    private:
        EstadoSistemaSoA estado_base_; // Estado original (sin ruido) extraído de las efemérides
        double dt_;                    // Delta de tiempo (paso de integración en días)
        double tiempo_maximo_;         // Límite temporal de la simulación (días a futuro)
        
        // Caché de los índices para no hacer búsquedas de strings en cada iteración
        std::size_t indice_tierra_;    
        std::size_t indice_asteroide_; 

        /**
         * @brief Calcula la aceleración gravitacional de todos los cuerpos.
         * Se marca como 'noexcept' para informar al compilador que no lanzará 
         * excepciones, permitiendo una optimización más agresiva del código ensamblador.
         */
        void calcularAceleraciones(EstadoSistemaSoA& estado) const noexcept;

        /**
         * @brief Evoluciona un universo individual a través del tiempo usando Verlet.
         * @param estado_local El universo específico con ruido estadístico aplicado.
         * @return El nombre del planeta impactado, o cadena vacía si no hubo impacto.
         */
        [[nodiscard]] std::string ejecutarSimulacionAislada(EstadoSistemaSoA& estado_local) const noexcept;

    public:
        /**
         * @brief Constructor del motor físico. Inicializa y valida la memoria.
         */
        MotorVerletMonteCarlo(EstadoSistemaSoA estado_inicial, double paso_tiempo, double t_max);

        /**
         * @brief Orquesta la simulación paralela inyectando ruido gaussiano (Monte Carlo).
         * @param num_corridas La cantidad de universos paralelos a simular.
         * @param sigma_pos Margen de error/incertidumbre en la posición inicial.
         * @param sigma_vel Margen de error/incertidumbre en la velocidad inicial.
         * @return El consolidado estadístico de las colisiones.
         */
        [[nodiscard]] ResultadosMonteCarlo ejecutarMonteCarlo(
            std::size_t num_corridas, 
            double sigma_pos, 
            double sigma_vel
        ) const;
    };

} // namespace SimuladorAstrofisico

#endif // SIMULADOR_IMPACTO_HPP