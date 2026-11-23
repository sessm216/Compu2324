# ================================================================================
# ANIMACION ISING  (versión que guarda la animación como vídeo .mp4)
#
# Genera una animación a partir de un fichero de datos con la configuración
# del retículo en cada instante de tiempo.
#
# El fichero debe estructurarse de la siguiente forma:
#
#   s(1,1)_1, s(1,2)_1, ..., s(1,M)_1
#   s(2,1)_1, s(2,2)_1, ..., s(2,M)_1
#   (...)
#   s(N,1)_1, s(N,2)_1, ..., s(N,M)_1
#
#   s(1,1)_2, s(1,2)_2, ..., s(1,M)_2
#   (...)
#
# donde s(i,j)_k es el valor del spin en la fila i-ésima y la columna
# j-ésima en el instante k. Los valores del spin deben ser +1 ó -1.
# Los bloques (instantes de tiempo) se separan con una línea en blanco.
#
# Si solo hay un instante de tiempo, se genera una imagen en pdf en lugar de
# una animación.
#
# Se puede configurar cambiando las variables de la sección "Parámetros".
# ================================================================================

import io
import os
import re

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# Si matplotlib no encuentra ffmpeg, usa el que trae el paquete imageio-ffmpeg
# (pip install imageio-ffmpeg). Si no está instalado, se usa el ffmpeg del sistema.
try:
    import imageio_ffmpeg
    matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

# Parámetros
# ========================================
file_in = "spins_v2.dat"  # Nombre del fichero de datos
file_out = "isingT23"               # Nombre del fichero de salida (sin extensión)
interval = 100                      # Tiempo entre fotogramas en milisegundos
save_to_file = True                 # True: guarda el vídeo mp4, False: lo muestra por pantalla
dpi = 150                           # Calidad del vídeo de salida (dots per inch)
every = 1                           # Usa 1 de cada "every" fotogramas (p. ej. 5 acorta el vídeo)

# Las rutas se resuelven respecto a la carpeta de este script,
# no respecto al directorio desde el que se ejecuta
base_dir = os.path.dirname(os.path.abspath(__file__))
path_in = os.path.join(base_dir, file_in)

# Lectura del fichero de datos
# ========================================
with open(path_in, "r") as f:
    data_str = f.read()

# Lista con los datos de cada fotograma: frames_data[j] es el fotograma j-ésimo
frames_data = list()

# Los bloques se separan por una o más líneas en blanco (también las que solo
# tienen espacios). Se ignoran los bloques vacíos.
for frame_data_str in re.split(r"\n\s*\n", data_str.strip()):
    # Quita espacios y la coma final de cada fila, y descarta filas vacías
    cleaned = "\n".join(
        line.strip().rstrip(",")
        for line in frame_data_str.splitlines()
        if line.strip()
    )
    if cleaned:
        frames_data.append(np.loadtxt(io.StringIO(cleaned), delimiter=","))

frames_data = frames_data[::every]
nframes = len(frames_data)
print(f"Leídos {nframes} fotogramas de {frames_data[0].shape}")

# Creación de la animación/gráfico
# ========================================
fig, ax = plt.subplots()
ax.axis("off")  # No muestra los ejes

# Representa el primer fotograma
im = ax.imshow(frames_data[0], cmap="binary", vmin=-1, vmax=+1)


# Función que actualiza la configuración del sistema en la animación
def update(j_frame, frames_data, im):
    im.set_data(frames_data[j_frame])
    return im,


# Si hay más de un instante de tiempo, genera la animación
if nframes > 1:
    animation = FuncAnimation(
        fig, update,
        fargs=(frames_data, im), frames=nframes, blit=True, interval=interval)

    if save_to_file:
        path_out = os.path.join(base_dir, "{}.mp4".format(file_out))
        # yuv420p: formato de píxel compatible con el Reproductor de Windows,
        # "Películas y TV", navegadores, etc.
        writer = FFMpegWriter(
            fps=1000 / interval,
            codec="libx264",
            extra_args=["-pix_fmt", "yuv420p"],
        )
        animation.save(path_out, writer=writer, dpi=dpi)
        print(f"Vídeo guardado en: {path_out}")
    else:
        plt.show()
# En caso contrario, muestra o guarda una imagen
else:
    if save_to_file:
        path_out = os.path.join(base_dir, "{}.pdf".format(file_out))
        fig.savefig(path_out)
        print(f"Imagen guardada en: {path_out}")
    else:
        plt.show()
