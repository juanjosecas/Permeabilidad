#!/usr/bin/env python
# coding: utf-8

import os
import psutil
import re
import time
from datetime import datetime

import MDAnalysis as mda
import numpy as np

# donde están los archivos
# si vienen de Gromacs, .gro es la topología y .xtc es
# la trayectoria
# si viene de NAMD don .psf y .dcd
GRO = 'Mamifero.gro'
XTC = 'Mamifero.xtc'
SUF = 'Mamifero'  # el sufijo de los archivos, o sea que todos se van a llamar permeation_[sufijo].csv
# u otra cosa de ser necesario, pero con ese sufijo

MAX = 20  # vamos a ver si sirve para corregir las PBC,
# tiene que ser muy cercano al valor de la longitud del poro que queremos analizar
# esto nos deja, de nuevo, con el problema de definir qué es el poro
# dónde empieza y donde termina. Para la AQP1 la longitud del poro, medida de diferentes formas,
# arroja alrededor de 20 angstrom

OUTPUT_FILE = "permeation_" + SUF + ".csv"
LOGFILE = "permeation_" + SUF + ".log"  # log file
MEMORY = False  # vamos a ver si evitamos cargar la cosa en la ram.


# Sirve si es muy grande la simulación y tenemos pocos recursos

# ordenar en forma humana
def atoi(text):
    return int(text) if text.isdigit( ) else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]


# función para truncar los números. Es cruda, simple y sirve para no poblar la pantalla de
# tantos números. Dudo mucho que la vaya a mejorar, porque funciona para algo rápido
# no redondea, sino que corta los decimales.
def truncate(n, dec=3):
    if n >= 1.0:
        multiplier = 10 ** dec
        return int(n * multiplier) / multiplier
    else:
        multiplier = 100 ** dec
        return int(n * multiplier) / multiplier


# printear información básica
print("Cargando los archivos de la trayectoria en la memoria...")
print('Estadísticas de los archivos')
gro_stat = os.stat(GRO)
xtc_stat = os.stat(XTC)
print(GRO + ":", truncate(gro_stat.st_size / (1024 * 1024)), "MB")
print(XTC + ":", truncate(xtc_stat.st_size / (1024 * 1024)), "MB")

# vamos a registrar el tiempo de ejecución
print('Hoy es', datetime.today( ).strftime('%Y-%m-%d %H:%M:%S'))  # la fecha y hora de hoy
start = time.time( )

# cargar los archivos en un Universo de MDAnalysis
u1 = mda.Universe(GRO, XTC, in_memory=MEMORY)

# informar qué hace
process = psutil.Process(os.getpid( ))  # recupera el ID del proceso para medir el impacto en el sistema
print("\n", truncate(psutil.virtual_memory( )[0] / 1000000), "MB RAM total")  # in megabytes
print(truncate(process.memory_info( ).rss / 1000000), "MB RAM ocupada por este proceso")  # in megabytes
print(psutil.virtual_memory( )[2], "% RAM ocupada")  # in megabytes
time.sleep(1)

# guardar las posiciones de las moléculas que nos interesan
print("Obteniendo coordenadas de las partes del sistema involucradas...")

# Donde empieza y donde termina la membrana, Se usa el átomo de fósforo.
# por qué no con la proteína? cambiaría? esto no es más genérico?
limit_memb = u1.select_atoms("name P")

# A ver si la membrana/sistema límite está centrado
xbox, ybox, zbox = u1.coord.dimensions[:3]

# el grupo superior e inferior
limit_sup = limit_memb.select_atoms("prop z > %f" % (zbox / 2), updating=True)  # seleccionar lo que está por arriba
limit_inf = limit_memb.select_atoms("prop z < %f" % (zbox / 2), updating=True)  # seleccionar lo que está por abajo

assert limit_inf.atoms.n_atoms != 0, "Eh, parece que no hay nada"

# Seleccionar las aguas. Acordarse de verificar los tipos de átomos.
# Si usamos TIP3 tiene como oxígeno a OH2 en CHARMM

solve_select = "resname TIP3 and name OH2"

# seleccionamos los oxígenos
waterox = u1.select_atoms(solve_select)

# los arrays con las posiciones de las aguas y los fósforos de los lípidos
nframes = u1.trajectory.n_frames - 1  # cuantos frames-1 porque el primero no cuenta ya que es la topología!
nb_waters = len(waterox.atoms.resids)  # cuantas aguas hay?
print("La selección de las moléculas se hizo con el comando de selección:", solve_select)
print('\nFrames para analizar:', nframes, 'n\Aguas a analizar:', nb_waters)

# arrays con valores iniciales, todos con 0 porque no hay información previa
z = np.zeros((nb_waters, nframes))  # un array de numpy lleno con ceros, porque al principio no pasa nada
old_z = np.zeros((nb_waters,
                  nframes))  # un array de numpy lleno con ceros, porque al principio no pasa nada pero ahora con la condición "vieja"
z_inf_traj = np.zeros(nframes)
z_sup_traj = np.zeros(nframes)

# unidades de tiempo, complicado si la trayectoria los tiene borrados,
# Las unidades por default de Gromacs son los picosegundos, los pasamos a nanosegundos (ns = ps*0.001).
# Si modificamos la simulación en un postprocesamiento, entonces tal vez esto arroje un error
dt_frame = u1.trajectory.ts.dt * 0.001

# captar las coordenadas "viejas"
old_z[:, 0] = waterox.positions[:, 2]

# de nuevo, enfatizo el excluir el primer frame
for iframe, ts in enumerate(u1.trajectory[1:]):
    # La media de la coordenada z de los fósforos
    z_inf_traj[iframe] = limit_inf.positions[:, 2].mean( )  # inferior
    z_sup_traj[iframe] = limit_sup.positions[:, 2].mean( )  # superior

    # recupero la coordenada z para el átomo de oxígeno x y z tienen el orden 0 1 2
    z[:, iframe] = waterox.positions[:, 2]

old_z[:, 1:nframes] = z[:, :nframes - 1]  # el viejo frame es el nuevo menos uno

# Chequeo la RAM a ver cómo viene
process = psutil.Process(os.getpid( ))
print("\n", truncate(psutil.virtual_memory( )[0] / 1000000), "MB RAM total\n")  # in megabytes
print(truncate(process.memory_info( ).rss / 1000000), "MB RAM ocupada por este proceso\n")  # in megabytes
print(psutil.virtual_memory( )[2], "% RAM ocupada\n")  # in megabytes
time.sleep(1)

# Tratemos de eliminar el problema de las PBC usando los datos de dónde está la membrana

test_disp = abs(z - old_z) > MAX  # actual menos el anterior es mayor a la longitud del poro
# algo interesante en esta idea porque establece que está lejos de la membrana y se asegura
# de corregir las PBC

# de nuevo, inicializamos con arrays de ceros y llenamos con True False para saber si están arriba o abajo
# es mejor que empezar a poner los z y luego transformarlos
# Numpy ocupa menos recursos que pandas, aunque no estoy seguro como para jurar por dios y los santos evangelios
is_inf = np.zeros((nb_waters, nframes), dtype=bool)  # está abajo?
is_sup = np.zeros((nb_waters, nframes), dtype=bool)  # está arriba?
is_rel = np.zeros((nb_waters, nframes), dtype=bool)  # es relativo a?

# la lógica del siguiente loop, en forma matemática es:
# inferior?: z < zinf, z de los P de abajo
# superior?: z > zsup, z de los P de arriba
# relativo? z < (zinf+zsup)/2
for i in range(z.shape[0]):
    is_inf[i] = z[i, :nframes] < z_inf_traj[:nframes]
    is_sup[i] = z[i, :nframes] > z_sup_traj[:nframes]
    is_rel[i] = (z[i, :nframes] < (z_inf_traj + z_sup_traj)[:nframes] / 2)

# con respecto a la membrana
# la función invert está relacionada con la conversión de números binarios
# demasiado críptica para explicar mejor, buscar en algún lugar tipo tutorial
# para que nos quede claro y no copiar y pegar ciegamente
is_memb = np.invert(is_inf + is_sup)
is_memb_inf = is_rel * np.invert(is_inf)
is_memb_sup = np.invert(is_rel) * np.invert(is_sup)

# Ahora la parte de llamar a las cosas por su nombre
# evidentemente no podemos escapar a la idea de Zhu
# así que voy a hacer lo mismo, pero en vez de 2,1,0,-1,-2
# que parece algo poco intuitivo, vamos a hacer:
'''
   -------------------- = límite superior del sistema
           O = átomo de oxígeno
        Arriba = 1
            O
    pppppp|   |pppppppp = fosfatos
          | O |
          |   |  Poro = 0
          | O |
    pppppp|   |pppppppp = fosfatos
            O    O
        Abajo = -1  O
   -------------------- = límite inferior del sistema

'''
# simplificación de la trayectoria
traj_simp = (is_inf[:, :nframes - 1] * -1  # si está debajo de la membrana
             + is_sup[:, :nframes - 1] * 1  # si está arriba de la membrana
             + is_memb[:, :nframes - 1] * 0  # si está en la membrana
             + test_disp[:, :nframes - 1] * (
                     is_memb_sup[:, :nframes - 1] * 1 + is_memb_inf[:, :nframes - 1] * -1)  # ok, dónde cambió
             + test_disp[:, 1:nframes] * (is_memb_sup[:, :nframes - 1] * 1 + is_memb_inf[:, :nframes - 1] * -1))

# vamos a procesar la trayectoria
print("Procesando trayectoria...")

stdout = open(OUTPUT_FILE, "w")
stdout.write("water_index resid time_per_frame(ns) time(frame) time(ns) duration(frames) direction\n")

nmol, nmoves = traj_simp.shape

time_perm = []  # lista de permeaciones, se inicia vacía
count = 0  # conteos iniciales
count_rel = 0

for index, value in enumerate(traj_simp[:, :]):
    history = np.zeros(nmoves, dtype=int)
    history[0] = value[0]
    label_previous = value[0]  # posición inicial

    for i, position in enumerate(value[:]):
        # historia de las posiciones ; -1 abajo, +1 arriba
        history[i] = 0  # empezar el registro de historia
        if (position == 1):
            history[i] = 1
        elif (position == -1):
            history[i] = -1
        else:
            if (label_previous > 0):
                history[i] = label_previous + 1
            elif (label_previous < 0):
                history[i] = label_previous - 1
        label_previous = history[i]

    # poblar el array de cálculos con el tipo de fenómeno
    diff = value[1:nmoves] - value[:nmoves - 1]

    # Cómo llamamos a cada evento
    # condiciones:
    # - la permeación ocurre cuando parte de la membrana (value-diff=0)
    # - el salto (diff) es de signo opuesto al anterior estado (history) 
    events_sparse = (value[1:] - diff == 0) * (history[:nmoves - 1] * diff < 0) * diff * (np.arange(nmoves - 1) + 2)
    events_bool = (events_sparse != 0)
    events = events_sparse[events_bool]
    duration_frames = abs((events_bool * history[:nmoves - 1] * diff)[events_bool])
    direction_perm = diff[events_bool]

    if events.size:
        for i, time_frame in enumerate(events):
            count += 1  # contar
            count_rel += direction_perm[i]
            resid, time_ns = waterox.resids[index], time_frame * dt_frame  # en unidades de tiempo
            # poblar el archivo de salida con los datos
            stdout.write("%5d %5d %5d %5d %6.3f %6.3f %5d\n"
                         % (
                             index, resid, dt_frame, abs(time_frame), abs(time_ns), duration_frames[i],
                             direction_perm[i]))

# cerrar el archivo!!!!
stdout.close( )

# registrar el final del conteo
end = time.time( )

print('\nHa finalizado el conteo de eventos de permeación')
print('Los cálculos tardaron', end - start, 'segundos en completarse')

# vamos a graficar!!!
import matplotlib.pyplot as plt
import pandas as pd

TIME_OFFSET = 0  # a partir de cuándo quiero ignorar los conteos

# cargo el sufijo en nueva variable.
# innecesario, pero las cosas son así porque no tengo ganas de cambiarlo
# además recuerdo tan abajo del script para qué sirve
FLAG_SIM = SUF

X_data = pd.read_csv(OUTPUT_FILE, sep=' ', header=0, skipinitialspace=True)

X = X_data[X_data["time(ns)"] > TIME_OFFSET]

# Suma total
X = X.sort_values(by='time(ns)')  # ordenados en función
X['permeations_tot'] = [i + 1 for i in range(len(X))]

print("Guardando en PNG, PDF y un dataset en formato CSV puesto lindo...\n")
print('Los archivos de salida son', 'permeation_' + FLAG_SIM + '_tot.png', 'y', 'permeation_' + FLAG_SIM + '_tot.pdf')

plt.figure( )
plt.step(X['time(ns)'], X['permeations_tot'], '-', where='post')
plt.title('Cumulated Total number of Permeations \n in %s' % FLAG_SIM, fontsize=18)
plt.xlabel('time(ns)', fontsize=16)
plt.ylabel('# permeations', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(TIME_OFFSET, )
plt.ylim(0, )

plt.savefig('permeation_' + FLAG_SIM + '_tot.png')
plt.savefig('permeation_' + FLAG_SIM + '_tot.pdf')
