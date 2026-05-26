import time
import sys

# Importamos la librería compilada dinámicamente en C++ (.pyd en Windows / .so en Linux)
try:
    import simulador_impacto
except ImportError as e:
    print(f"Error crítico: No se encontró el módulo C++ compilado. Revisa que el archivo exista.\nDetalles: {e}")
    sys.exit(1)

def crear_condiciones_iniciales():
    """
    Construye las condiciones límite basadas en efemérides astronómicas.
    Unidades físicas obligatorias:
    - Masas: Masas Solares (M_sun = 1.0)
    - Distancias: Unidades Astronómicas (UA)
    - Tiempos: Días terrestres
    Modelo Idealizado (Elementos Medios J2000).
    Datos referenciados de NASA Planetary Fact Sheet.
    """
    estado = simulador_impacto.EstadoSistemaSoA()
    
    estado.num_cuerpos = 4
    estado.nombres = ["Sol", "Tierra", "Jupiter", "Asteroide"]
    
    # Masas (Relativas al Sol)
    estado.masa = [
        1.0,                    # Sol
        3.0027e-6,              # Tierra (Masa real)
        9.543e-4,               # Júpiter (Masa real)
        0.0                     # Asteroide (Partícula de prueba)
    ]
    
    # Radios físicos en UA
    estado.radio = [
        0.00465,                # Sol
        4.2588e-5,              # Tierra
        4.67e-4,                # Júpiter
        2.0e-7                  # Asteroide genérico
    ]

    if True:
        # Posiciones J2000 Idealizadas (Plano Eclíptico X-Y)
        # Colocamos a los planetas en el eje X
        estado.x = [
            0.0,                    # Sol en el origen
            1.0,                    # Tierra exactamente a 1 UA
            5.2,                    # Júpiter a 5.2 UA
            1.05                    # Asteroide iniciando un poco más allá de la Tierra
        ]
        estado.y = [0.0, 0.0, 0.0, 0.05] # Asteroide ligeramente desfasado en Y
        estado.z = [0.0, 0.0, 0.0, 0.0]  # Sistema coplanar ideal
        
        # Velocidades orbitales ideales para órbitas circulares (Eje Y)
        # v = sqrt(G * M_sol / r) -> Convertido a UA/día
        estado.vx = [
            0.0, 
            0.0, 
            0.0, 
            -0.015  # Asteroide moviéndose hacia adentro (hacia el Sol)
        ]
        estado.vy = [
            0.0, 
            0.0172, # Tierra (~30 km/s)
            0.0075, # Júpiter (~13 km/s)
            0.0100  # Asteroide con velocidad transversal
        ]
        estado.vz = [0.0, 0.0, 0.0, 0.0]
        
        # Inicialización de aceleraciones
        estado.ax = [0.0] * 4
        estado.ay = [0.0] * 4
        estado.az = [0.0] * 4
        
    else:
        # === ESCENARIO DE CHOQUE SEGURO ===
        # Tomamos las coordenadas exactas de la Tierra
        x_tierra = 0.983
        y_tierra = 0.174
        z_tierra = -0.0001
        
        # Colocamos al asteroide a una distancia de 0.00001 UA (1,500 km) de la Tierra.
        # El radio de la Tierra es 0.00004 UA (~6,300 km). 
        # Por lo tanto, el asteroide literalmente está naciendo en el manto terrestre.
        estado.x = [0.0, x_tierra,  5.20, x_tierra + 0.00001]
        estado.y = [0.0, y_tierra,  0.00, y_tierra]
        estado.z = [0.0, z_tierra, -0.12, z_tierra]
        
        # Velocidades (las igualamos para que no haya escape)
        estado.vx = [0.0, -0.003,  0.000, -0.003]
        estado.vy = [0.0,  0.0172, 0.0075, 0.0172]
        estado.vz = [0.0,  0.000,  0.000,  0.000]
        
        estado.ax = [0.0] * 4
        estado.ay = [0.0] * 4
        estado.az = [0.0] * 4
            
    return estado

def main():
    print("="*60)
    print(" MOTOR HPC: SIMULACIÓN DE IMPACTOS MEDIANTE MONTE CARLO")
    print(" Proyecto Integrador - Laboratorio de Mecánica Clásica")
    print("="*60)
    
    # [FASE 1] Carga de memoria RAM (C++) desde variables en el intérprete (Python)
    print("\n[1] Cargando condiciones iniciales en el formato Structure of Arrays...")
    estado_base = crear_condiciones_iniciales()
    
    # [FASE 2] Configuración temporal de la Ecuación Diferencial
    dt_dias = 0.05                 # Resolución del integrador: Una iteración cada 1.2 horas físicas.
    tiempo_total_dias = 365.25 * 10 # Tiempo futuro de proyección: 10 años (considerando bisiestos).
    
    print(f"[2] Configurando integrador simpléctico (Velocity Verlet)...")
    print(f"    - dt: {dt_dias} días terrestres")
    print(f"    - Límite temporal: {tiempo_total_dias} días (10 años)")
    
    # Instanciamos el objeto que contiene la lógica astrofísica y el control de OpenMP
    motor = simulador_impacto.MotorVerletMonteCarlo(estado_base, dt_dias, tiempo_total_dias)
    
    # [FASE 3] Parametrización del Ruido Estocástico (Campana de Gauss)
    num_simulaciones = 100000  # Universo estadístico
    
    # Desviación estándar impuesta a las observaciones iniciales del asteroide.
    # Un sigma_pos mayor crea una 'nube de probabilidad' más expansiva.
    if True:
        sigma_pos = 1e-4  # Incertidumbre en posición (~15,000 km)
        sigma_vel = 1e-5  # Incertidumbre en la velocidad telescópica 
    else:
        # Incertidumbre nula
        sigma_pos = 0.0
        sigma_vel = 0.0
    
    print(f"[3] Liberando GIL e iniciando procesamiento asíncrono con OpenMP...")
    print(f"    - Ciclos Monte Carlo a ejecutar: {num_simulaciones:,}")
    
    # Control del cronómetro de hardware para evaluar rendimiento HPC
    inicio_cómputo = time.perf_counter()
    
    # [FASE 4] EJECUCIÓN PURA (Se delega el 100% de los núcleos del procesador al backend de C++)
    resultados = motor.ejecutar_monte_carlo(
        num_corridas=num_simulaciones, 
        sigma_pos=sigma_pos, 
        sigma_vel=sigma_vel
    )
    
    fin_cómputo = time.perf_counter()
    tiempo_ejecucion = fin_cómputo - inicio_cómputo
    
    # [FASE 5] Reporte y Formateo Analítico
    print("\n" + "="*60)
    print(" RESULTADOS ESTADÍSTICOS FINALES")
    print("="*60)
    
    print(f"-> Tiempo en CPU: {tiempo_ejecucion:.4f} segundos")
    print(f"-> Tasa de rendimiento físico:   {(num_simulaciones / tiempo_ejecucion):,.0f} universos/segundo\n")
    
    print("Reporte de Early Exits (Choques):")
    if not resultados.desglose_colisiones:
        print("  - Ningún asteroide intersectó barreras físicas durante los 10 años simulados.")
    else:
        for planeta, conteo in resultados.desglose_colisiones.items():
            porcentaje_local = (conteo / num_simulaciones) * 100
            print(f"  - [{planeta}] colisiones absolutas: {conteo:,} ({porcentaje_local:.4f}%)")
            
    print("-" * 60)
    # Extracción de la probabilidad porcentual 
    prob = resultados.probabilidad_impacto * 100
    print(f" PROBABILIDAD REAL DE IMPACTO EN LA TIERRA: {prob:.4f} %")
    print("="*60)

if __name__ == "__main__":
    main()