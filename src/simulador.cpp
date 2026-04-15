#include "simulador.hpp"
#include <cmath>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace SimuladorAstrofisico {

// ========================================================================
// Constructor: Inicialización y Validación
// ========================================================================
MotorVerletMonteCarlo::MotorVerletMonteCarlo(EstadoSistemaSoA estado_inicial, double paso_tiempo, double t_max)
    : estado_base_(std::move(estado_inicial)), dt_(paso_tiempo), tiempo_maximo_(t_max) {
    
    std::size_t N = estado_base_.num_cuerpos;

    // Validación de integridad del SoA (Previene SegFaults por mala memoria)
    if (N < 2) {
        throw std::invalid_argument("El sistema requiere al menos 2 cuerpos.");
    }
    if (estado_base_.x.size() != N || estado_base_.masa.size() != N || estado_base_.nombres.size() != N) {
        throw std::length_error("Inconsistencia en el tamaño de los arreglos SoA.");
    }

    // Cachear los índices críticos para evitar búsquedas de strings dentro de los ciclos de física
    auto it_tierra = std::find(estado_base_.nombres.begin(), estado_base_.nombres.end(), "Tierra");
    auto it_asteroide = std::find(estado_base_.nombres.begin(), estado_base_.nombres.end(), "Asteroide"); // O "Apophis"

    if (it_tierra == estado_base_.nombres.end() || it_asteroide == estado_base_.nombres.end()) {
        throw std::runtime_error("No se encontró 'Tierra' o 'Asteroide' en los nombres provistos.");
    }

    indice_tierra_ = std::distance(estado_base_.nombres.begin(), it_tierra);
    indice_asteroide_ = std::distance(estado_base_.nombres.begin(), it_asteroide);

    // Calcular aceleraciones iniciales en el estado base para arrancar Verlet correctamente
    calcularAceleraciones(estado_base_);
}

// ========================================================================
// Motor Físico: Ley de Gravitación Universal Vectorizada
// ========================================================================
void MotorVerletMonteCarlo::calcularAceleraciones(EstadoSistemaSoA& estado) const noexcept {
    const std::size_t N = estado.num_cuerpos;

    for (std::size_t i = 0; i < N; ++i) {
        double ax = 0.0, ay = 0.0, az = 0.0;
        
        // Copiamos la posición del cuerpo actual a registros rápidos
        const double xi = estado.x[i];
        const double yi = estado.y[i];
        const double zi = estado.z[i];

        // Directiva OpenMP para forzar vectorización SIMD
        #pragma omp simd reduction(+:ax, ay, az)
        for (std::size_t j = 0; j < N; ++j) {
            if (i == j) continue; // Un cuerpo no se atrae a sí mismo

            double dx = estado.x[j] - xi;
            double dy = estado.y[j] - yi;
            double dz = estado.z[j] - zi;

            double dist_sq = dx*dx + dy*dy + dz*dz;
            
            // 1.0 / sqrt() es frecuentemente optimizado por el compilador a 
            // la instrucción de hardware rsqrt (Reciprocal Square Root).
            double inv_dist = 1.0 / std::sqrt(dist_sq);
            double inv_dist3 = inv_dist * inv_dist * inv_dist;

            // f contiene a G y la masa. Omitimos la masa_i porque se cancela con F=ma
            double f = G_GAUSSIANA * estado.masa[j] * inv_dist3;

            ax += f * dx;
            ay += f * dy;
            az += f * dz;
        }

        estado.ax[i] = ax;
        estado.ay[i] = ay;
        estado.az[i] = az;
    }
}

// ========================================================================
// Simulación de un Universo (Iteración de Verlet y Early Exit)
// ========================================================================
std::string MotorVerletMonteCarlo::ejecutarSimulacionAislada(EstadoSistemaSoA& estado_local) const noexcept {
    const std::size_t N = estado_local.num_cuerpos;
    double t = 0.0;

    // Constantes precalculadas
    const double dt_medio = 0.5 * dt_;
    const std::size_t idx_ast = indice_asteroide_;

    while (t < tiempo_maximo_) {
        // 1. Verlet: Avanzar posiciones y medio paso de velocidad
        #pragma omp simd
        for (std::size_t i = 0; i < N; ++i) {
            estado_local.vx[i] += estado_local.ax[i] * dt_medio;
            estado_local.vy[i] += estado_local.ay[i] * dt_medio;
            estado_local.vz[i] += estado_local.az[i] * dt_medio;

            estado_local.x[i] += estado_local.vx[i] * dt_;
            estado_local.y[i] += estado_local.vy[i] * dt_;
            estado_local.z[i] += estado_local.vz[i] * dt_;
        }

        // 2. Early Exit: Detección de Colisión
        const double x_ast = estado_local.x[idx_ast];
        const double y_ast = estado_local.y[idx_ast];
        const double z_ast = estado_local.z[idx_ast];
        const double r_ast = estado_local.radio[idx_ast];

        for (std::size_t i = 0; i < N; ++i) {
            if (i == idx_ast) continue;

            double dx = estado_local.x[i] - x_ast;
            double dy = estado_local.y[i] - y_ast;
            double dz = estado_local.z[i] - z_ast;
            double dist_sq = dx*dx + dy*dy + dz*dz;

            double limite_colision = r_ast + estado_local.radio[i];
            
            // Comparar cuadrados evita un std::sqrt() costoso
            if (dist_sq <= limite_colision * limite_colision) {
                return estado_local.nombres[i]; // Choque detectado, retornar el nombre del planeta
            }
        }

        // 3. Verlet: Actualizar aceleraciones con las nuevas posiciones
        calcularAceleraciones(estado_local);

        // 4. Verlet: Completar el medio paso restante de velocidad
        #pragma omp simd
        for (std::size_t i = 0; i < N; ++i) {
            estado_local.vx[i] += estado_local.ax[i] * dt_medio;
            estado_local.vy[i] += estado_local.ay[i] * dt_medio;
            estado_local.vz[i] += estado_local.az[i] * dt_medio;
        }

        t += dt_;
    }

    return ""; // Sobrevivió el tiempo máximo sin colisiones
}

// ========================================================================
// Orquestador Monte Carlo (Multithreading y Estadística)
// ========================================================================
ResultadosMonteCarlo MotorVerletMonteCarlo::ejecutarMonteCarlo(
    std::size_t num_corridas, double sigma_pos, double sigma_vel) const {
    
    ResultadosMonteCarlo resultados = {num_corridas, 0, 0.0, {}};

    // Se inicializa un motor de base real para dispersar las semillas de los hilos
    std::random_device rd;
    unsigned int base_seed = rd();

    #pragma omp parallel
    {
        // Thread-Local Storage (TLS): Aislamiento absoluto de variables por hilo
        // Evita el data race y el falso uso compartido (false sharing) de caché
        std::map<std::string, std::size_t> mapa_impactos_local;
        std::size_t impactos_tierra_local = 0;

        // Generador estadístico aislado por hilo. 
        // La semilla usa el ID del hilo para garantizar divergencia estocástica.
        std::mt19937_64 rng(base_seed + omp_get_thread_num());
        std::normal_distribution<double> dist_pos(0.0, sigma_pos);
        std::normal_distribution<double> dist_vel(0.0, sigma_vel);

        long long total_corridas = static_cast<long long>(num_corridas);
        
        // Repartición dinámica del trabajo entre los núcleos del CPU
        #pragma omp for schedule(dynamic)
        for (long long corrida = 0; corrida < total_corridas; ++corrida) {
            
            // Generar el universo paralelo local copiando el estado limpio
            EstadoSistemaSoA estado_universo = estado_base_;

            // Aplicar el ruido gaussiano a las condiciones del asteroide
            estado_universo.x[indice_asteroide_] += dist_pos(rng);
            estado_universo.y[indice_asteroide_] += dist_pos(rng);
            estado_universo.z[indice_asteroide_] += dist_pos(rng);

            estado_universo.vx[indice_asteroide_] += dist_vel(rng);
            estado_universo.vy[indice_asteroide_] += dist_vel(rng);
            estado_universo.vz[indice_asteroide_] += dist_vel(rng);

            // Correr la simulación y capturar colisión
            std::string colision_resultado = ejecutarSimulacionAislada(estado_universo);

            if (!colision_resultado.empty()) {
                mapa_impactos_local[colision_resultado]++;
                if (colision_resultado == "Tierra") {
                    impactos_tierra_local++;
                }
            }
        }

        // Sección Crítica: Agrupar resultados al mapa global y contador global
        // Solo un hilo puede ejecutar este bloque a la vez. Es rápido porque 
        // solo ocurre una vez por hilo al finalizar su lote de trabajo.
        #pragma omp critical
        {
            resultados.impactos_tierra += impactos_tierra_local;
            for (const auto& [planeta, conteo] : mapa_impactos_local) {
                resultados.desglose_colisiones[planeta] += conteo;
            }
        }
    }

    // Post-procesamiento estadístico
    resultados.probabilidad_impacto = 
        static_cast<double>(resultados.impactos_tierra) / static_cast<double>(num_corridas);

    return resultados;
}

} // namespace SimuladorAstrofisico