import time
import sys

# Importamos nuestro motor C++ de alto rendimiento
try:
    import simulador_impacto
except ImportError as e:
    print(f"Error al importar el módulo C++. Asegúrate de que el archivo .pyd/.so esté en esta carpeta.\nDetalle: {e}")
    sys.exit(1)

def crear_condiciones_iniciales():
    """
    Carga los vectores de estado (Posición y Velocidad) aproximados
    para la época J2000. 
    Unidades: Masas Solares, Unidades Astronómicas (UA) y Días.
    """
    estado = simulador_impacto.EstadoSistemaSoA()
    
    # 4 Cuerpos: Sol, Tierra, Júpiter, Asteroide
    estado.num_cuerpos = 4
    estado.nombres = ["Sol", "Tierra", "Jupiter", "Asteroide"]
    
    # Masas en Masas Solares
    estado.masa = [
        1.0,                    # Sol
        3.0027e-6,              # Tierra
        9.543e-4,               # Júpiter
        0.0                     # Asteroide (masa despreciable para su propia cinemática)
    ]
    
    # Radios físicos en UA (Cruciales para la detección de colisión)
    estado.radio = [
        0.00465,                # Sol (~696,000 km)
        4.2588e-5,              # Tierra (~6,371 km)
        4.67e-4,                # Júpiter (~69,911 km)
        2.0e-7                  # Asteroide grande (~30 km)
    ]

    # Posiciones iniciales (X, Y, Z) en UA
    estado.x = [0.0,  0.983,   5.20,  0.995]
    estado.y = [0.0,  0.174,   0.00,  0.050]
    estado.z = [0.0, -0.0001, -0.12,  0.001]
    
    # Velocidades iniciales (Vx, Vy, Vz) en UA/día
    # La Tierra viaja a ~0.0172 UA/día (unos 30 km/s)
    estado.vx = [0.0, -0.003,   0.000, -0.002]
    estado.vy = [0.0,  0.0172,  0.0075, 0.0165]
    estado.vz = [0.0,  0.000,   0.000,  0.001]
    
    # Inicializamos las aceleraciones en 0.0
    estado.ax = [0.0] * estado.num_cuerpos
    estado.ay = [0.0] * estado.num_cuerpos
    estado.az = [0.0] * estado.num_cuerpos
    
    return estado

def main():
    print("="*60)
    print(" MOTOR DE SIMULACIÓN MONTE CARLO - N-CUERPOS (C++ / OpenMP)")
    print("="*60)
    
    # 1. Preparar datos
    print("[1] Cargando efemérides base...")
    estado_base = crear_condiciones_iniciales()
    
    # 2. Configurar el integrador físico
    # dt = 0.05 días (~1.2 horas por iteración)
    # t_max = 365.25 * 10: Simular 10 años en el futuro
    dt_dias = 0.05
    tiempo_total_dias = 365.25 * 10 
    
    print(f"[2] Inicializando integrador Velocity Verlet...")
    print(f"    - Paso de tiempo (dt): {dt_dias} días")
    print(f"    - Tiempo a simular:    {tiempo_total_dias} días (10 años)")
    
    motor = simulador_impacto.MotorVerletMonteCarlo(estado_base, dt_dias, tiempo_total_dias)
    
    # 3. Configurar la estadística de Monte Carlo
    num_simulaciones = 100_000  # Numero de universos a simular
    
    # Ruido estadístico
    # sigma_pos = 1e-4 UA (~15,000 km de error de medición en la posición inicial)
    # sigma_vel = 1e-5 UA/día de error en la velocidad detectada
    sigma_pos = 1e-4 
    sigma_vel = 1e-5 
    
    print(f"[3] Disparando Monte Carlo en paralelo (C++ liberando el GIL)...")
    print(f"    - Simulaciones: {num_simulaciones:,}")
    print(f"    - Ruido Posición: {sigma_pos} UA | Ruido Velocidad: {sigma_vel} UA/d")
    
    # 4. EJECUCIÓN DEL MOTOR
    inicio_cómputo = time.perf_counter()
    
    resultados = motor.ejecutar_monte_carlo(
        num_corridas=num_simulaciones, 
        sigma_pos=sigma_pos, 
        sigma_vel=sigma_vel
    )
    
    fin_cómputo = time.perf_counter()
    tiempo_ejecucion = fin_cómputo - inicio_cómputo
    
    # 5. Análisis de Resultados
    print("\n" + "="*60)
    print(" RESULTADOS DEL ANÁLISIS ESTADÍSTICO")
    print("="*60)
    
    print(f"Tiempo de Cómputo CPU: {tiempo_ejecucion:.4f} segundos")
    print(f"Velocidad de Simulación: {(num_simulaciones / tiempo_ejecucion):,.0f} universos/segundo\n")
    
    print("Desglose de Colisiones (Early Exits):")
    if not resultados.desglose_colisiones:
        print("  - Ninguna colisión registrada en ninguna simulación.")
    else:
        for planeta, conteo in resultados.desglose_colisiones.items():
            porcentaje = (conteo / num_simulaciones) * 100
            print(f"  - Impactos en {planeta}: {conteo:,} ({porcentaje:.4f}%)")
            
    print("-" * 60)
    # Mostramos el resultado vital
    prob = resultados.probabilidad_impacto * 100
    print(f"🚀 PROBABILIDAD DE IMPACTO EN LA TIERRA: {prob:.15f} %")
    print("="*60)

if __name__ == "__main__":
    main()