#include "simulador.hpp"
#include <cmath>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace SimuladorAstrofisico {

// ========================================================================
// INICIALIZACIÓN Y VALIDACIÓN DEL SISTEMA
// ========================================================================
MotorVerletMonteCarlo::MotorVerletMonteCarlo(EstadoSistemaSoA estado_inicial, double paso_tiempo, double t_max)
    : estado_base_(std::move(estado_inicial)), dt_(paso_tiempo), tiempo_maximo_(t_max) {
    
    std::size_t N = estado_base_.num_cuerpos;

    // Validación de seguridad: Evitar que el programa falle por falta de datos
    if (N < 2) {
        throw std::invalid_argument("Error: Se necesitan al menos 2 cuerpos (ej. Sol y Asteroide).");
    }
    // Verificar que todos los vectores tengan el mismo tamaño (consistencia de SoA)
    if (estado_base_.x.size() != N || estado_base_.masa.size() != N || estado_base_.nombres.size() != N) {
        throw std::length_error("Error: Tamaño inconsistente en los arreglos de datos del sistema.");
    }

    // Buscamos los índices de la Tierra y el Asteroide una sola vez al inicio.
    // Hacer esto dentro del ciclo de la física arruinaría el rendimiento computacional.
    auto it_tierra = std::find(estado_base_.nombres.begin(), estado_base_.nombres.end(), "Tierra");
    auto it_asteroide = std::find(estado_base_.nombres.begin(), estado_base_.nombres.end(), "Asteroide");

    if (it_tierra == estado_base_.nombres.end() || it_asteroide == estado_base_.nombres.end()) {
        throw std::runtime_error("Error: Faltan los identificadores 'Tierra' o 'Asteroide' en los datos.");
    }

    // Guardamos la posición en memoria de los cuerpos críticos
    indice_tierra_ = std::distance(estado_base_.nombres.begin(), it_tierra);
    indice_asteroide_ = std::distance(estado_base_.nombres.begin(), it_asteroide);

    // Arrancamos el sistema calculando la aceleración inicial (t=0) para que 
    // el algoritmo de Verlet tenga los datos base para empezar a iterar.
    calcularAceleraciones(estado_base_);
}

// ========================================================================
// NÚCLEO DINÁMICO: CÁLCULO DE FUERZAS GRAVITACIONALES
// ========================================================================
void MotorVerletMonteCarlo::calcularAceleraciones(EstadoSistemaSoA& estado) const noexcept {
    const std::size_t N = estado.num_cuerpos;

    // Evaluamos la atracción gravitacional mutua para cada cuerpo en el sistema
    for (std::size_t i = 0; i < N; ++i) {
        double ax = 0.0, ay = 0.0, az = 0.0;
        
        // Guardamos las coordenadas locales en variables temporales para acceso rápido
        const double xi = estado.x[i];
        const double yi = estado.y[i];
        const double zi = estado.z[i];

        // Instruimos al procesador a vectorizar este ciclo forzando operaciones SIMD.
        // Reduction suma los resultados paralelos en las variables de aceleración.
        #pragma omp simd reduction(+:ax, ay, az)
        for (std::size_t j = 0; j < N; ++j) {
            // Físicamente, un cuerpo no puede ejercer fuerza de gravedad sobre sí mismo
            if (i == j) continue; 

            // Distancias relativas en los tres ejes espaciales
            double dx = estado.x[j] - xi;
            double dy = estado.y[j] - yi;
            double dz = estado.z[j] - zi;

            // Distancia al cuadrado (Teorema de Pitágoras en 3D)
            double dist_sq = dx*dx + dy*dy + dz*dz;
            
            // Calculamos el inverso de la raíz cuadrada. Los procesadores modernos 
            // tienen instrucciones dedicadas para esto, haciéndolo extremadamente rápido.
            double inv_dist = 1.0 / std::sqrt(dist_sq);
            double inv_dist3 = inv_dist * inv_dist * inv_dist; // Inverso al cubo

            // Aplicamos la Ley de Gravitación de Newton: a = G * M / r^2
            // Notar que la masa propia (i) se cancela, quedando solo la masa perturbadora (j)
            double f = G_GAUSSIANA * estado.masa[j] * inv_dist3;

            // Proyectamos la magnitud de la aceleración sobre los componentes vectoriales
            ax += f * dx;
            ay += f * dy;
            az += f * dz;
        }

        // Actualizamos el estado de aceleración en la memoria principal
        estado.ax[i] = ax;
        estado.ay[i] = ay;
        estado.az[i] = az;
    }
}

// ========================================================================
// INTEGRACIÓN NUMÉRICA: EL ALGORITMO VELOCITY VERLET
// ========================================================================
std::string MotorVerletMonteCarlo::ejecutarSimulacionAislada(EstadoSistemaSoA& estado_local) const noexcept {
    const std::size_t N = estado_local.num_cuerpos;
    double t = 0.0;

    // Se calculan las constantes fuera del ciclo 'while' para no repetir 
    // multiplicaciones innecesarias durante la integración.
    const double dt_medio = 0.5 * dt_;
    const std::size_t idx_ast = indice_asteroide_;
    const double r_ast = estado_local.radio[idx_ast];

    const double x_ast_init = estado_local.x[idx_ast];
    const double y_ast_init = estado_local.y[idx_ast];
    const double z_ast_init = estado_local.z[idx_ast];

    for (std::size_t i = 0; i < N; ++i) {
        if (i == idx_ast) continue;
        double dx = estado_local.x[i] - x_ast_init;
        double dy = estado_local.y[i] - y_ast_init;
        double dz = estado_local.z[i] - z_ast_init;
        double dist_sq = dx*dx + dy*dy + dz*dz;
        
        double limite_colision = r_ast + estado_local.radio[i];
        if (dist_sq <= limite_colision * limite_colision) {
            return estado_local.nombres[i]; // Nació chocando
        }
    }

    // Ciclo de evolución temporal del universo actual
    while (t < tiempo_maximo_) {
        
        // PASO 1 de Verlet: Calculamos la velocidad a medio paso (t + dt/2) 
        // y actualizamos la posición al paso completo (t + dt).
        #pragma omp simd
        for (std::size_t i = 0; i < N; ++i) {
            estado_local.vx[i] += estado_local.ax[i] * dt_medio;
            estado_local.vy[i] += estado_local.ay[i] * dt_medio;
            estado_local.vz[i] += estado_local.az[i] * dt_medio;

            estado_local.x[i] += estado_local.vx[i] * dt_;
            estado_local.y[i] += estado_local.vy[i] * dt_;
            estado_local.z[i] += estado_local.vz[i] * dt_;
        }

        // CONDICIÓN DE PARADA TEMPRANA (Early Exit): Detección de Colisión Geométrica
        // Obtenemos las coordenadas y el radio físico del asteroide en el tiempo actual
        const double x_ast = estado_local.x[idx_ast];
        const double y_ast = estado_local.y[idx_ast];
        const double z_ast = estado_local.z[idx_ast];
        const double r_ast = estado_local.radio[idx_ast];

        // Verificamos si el asteroide chocó contra algún planeta masivo
        for (std::size_t i = 0; i < N; ++i) {
            if (i == idx_ast) continue;

            double dx = estado_local.x[i] - x_ast;
            double dy = estado_local.y[i] - y_ast;
            double dz = estado_local.z[i] - z_ast;
            double dist_sq = dx*dx + dy*dy + dz*dz;

            double limite_colision = r_ast + estado_local.radio[i];
            
            // Comparamos el cuadrado de la distancia contra el cuadrado del límite.
            // Esto se hace para evitar calcular std::sqrt(), ahorrando ciclos de CPU valiosos.
            if (dist_sq <= limite_colision * limite_colision) {
                return estado_local.nombres[i]; // Impacto detectado, abortamos la simulación
            }
        }

        // PASO 2 de Verlet: Recalculamos el campo de fuerzas (aceleración) 
        // usando las nuevas posiciones generadas en el Paso 1.
        calcularAceleraciones(estado_local);

        // PASO 3 de Verlet: Completamos el ciclo de la velocidad con el dt/2 restante
        // usando las nuevas aceleraciones calculadas en el Paso 2.
        #pragma omp simd
        for (std::size_t i = 0; i < N; ++i) {
            estado_local.vx[i] += estado_local.ax[i] * dt_medio;
            estado_local.vy[i] += estado_local.ay[i] * dt_medio;
            estado_local.vz[i] += estado_local.az[i] * dt_medio;
        }

        // Avanzamos el reloj de la simulación
        t += dt_;
    }

    return ""; // Si termina el bucle while sin retornar, no hubo colisiones.
}

// ========================================================================
// ORQUESTADOR ESTADÍSTICO (MONTE CARLO PARALELO)
// ========================================================================
ResultadosMonteCarlo MotorVerletMonteCarlo::ejecutarMonteCarlo(
    std::size_t num_corridas, double sigma_pos, double sigma_vel) const {
    
    // Contenedor global de resultados para devolver a Python
    ResultadosMonteCarlo resultados = {num_corridas, 0, 0.0, {}};

    // Se inicializa un motor aleatorio hardware/físico para servir de semilla base
    std::random_device rd;
    unsigned int base_seed = rd();

    // Iniciamos la región asíncrona. A partir de este bloque, OpenMP levanta todos los hilos.
    #pragma omp parallel
    {
        // Thread-Local Storage (Memoria Privada del Hilo):
        // Creamos contadores locales para que los hilos no choquen al intentar escribir
        // sobre las mismas variables simultáneamente (evitando 'data races').
        std::map<std::string, std::size_t> mapa_impactos_local;
        std::size_t impactos_tierra_local = 0;

        // Generador de ruido estocástico privado para este hilo. 
        // Usamos el ID del hilo para asegurar que las secuencias pseudoaleatorias no se repitan.
        std::mt19937_64 rng(base_seed + omp_get_thread_num());
        
        // Creamos las distribuciones gaussianas para perturbar las condiciones iniciales
        std::normal_distribution<double> dist_pos(0.0, sigma_pos);
        std::normal_distribution<double> dist_vel(0.0, sigma_vel);

        // Casteamos a signed int (long long) para que el compilador MSVC en Windows lo acepte
        long long total_corridas = static_cast<long long>(num_corridas);

        // Repartimos dinámicamente las iteraciones (universos) entre los procesadores
        #pragma omp for schedule(dynamic)
        for (long long corrida = 0; corrida < total_corridas; ++corrida) {
            
            // Creamos una copia fresca del sistema solar para esta simulación particular
            EstadoSistemaSoA estado_universo = estado_base_;

            // INYECCIÓN DE CAOS (Monte Carlo): Alteramos la posición inicial del asteroide
            estado_universo.x[indice_asteroide_] += dist_pos(rng);
            estado_universo.y[indice_asteroide_] += dist_pos(rng);
            estado_universo.z[indice_asteroide_] += dist_pos(rng);

            // Inyectamos caos en los vectores de velocidad inicial del asteroide
            estado_universo.vx[indice_asteroide_] += dist_vel(rng);
            estado_universo.vy[indice_asteroide_] += dist_vel(rng);
            estado_universo.vz[indice_asteroide_] += dist_vel(rng);

            // Desplegamos la física y esperamos a ver qué sucede
            std::string colision_resultado = ejecutarSimulacionAislada(estado_universo);

            // Si el motor físico devuelve un nombre, registramos el evento trágico
            if (!colision_resultado.empty()) {
                mapa_impactos_local[colision_resultado]++;
                if (colision_resultado == "Tierra") {
                    impactos_tierra_local++;
                }
            }
        }

        // REDUCCIÓN CRÍTICA: Los hilos terminan su trabajo y consolidan la información.
        // El pragma 'critical' actúa como un embudo: solo deja pasar a un hilo a la vez
        // para que sumen sus resultados locales al contenedor global de manera segura.
        #pragma omp critical
        {
            resultados.impactos_tierra += impactos_tierra_local;
            for (const auto& [planeta, conteo] : mapa_impactos_local) {
                resultados.desglose_colisiones[planeta] += conteo;
            }
        }
    }

    // Cálculo estadístico final: Probabilidad P(A) = Casos Favorables / Casos Posibles
    resultados.probabilidad_impacto = 
        static_cast<double>(resultados.impactos_tierra) / static_cast<double>(num_corridas);

    return resultados;
}

} // namespace SimuladorAstrofisico