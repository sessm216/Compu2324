# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:20:20 2024

@author: lucia
"""


#************HACER LAS GRAFICAS*************

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

nc=50
# Valores de x
l_1 = np.linspace(0, 0.999999, 100)  # 100 puntos entre 0 y 10
y_teoricos_1 = (4*(1-l_1))/((4*(1-l_1))+(l_1**2*(np.sin((2*np.pi/5)*nc*(1-l_1)**0.5))**2))


l_2 = np.linspace(1.00000001, 10.5, 100)  # 100 puntos entre 0 y 10


# Valores teóricos (por ejemplo, y = sin(x))
y_teoricos_2 = (4*(l_2-1))/((4*(l_2-1))+(l_2**2*(np.sinh((2*np.pi/5)*nc*(l_2-1)**0.5))**2))

print(y_teoricos_1)

# Valores experimentales (supongamos que tienes 10 puntos experimentales)
x_experimentales = [0.1,0.3,0.5,1,5,10]
y_experimentales = [0.981,0.947,0.9879,0.248,0.003,0]

plt.figure(figsize=(10, 6))

# Graficar los valores teóricos
plt.plot(l_1, y_teoricos_1, color='green', linestyle='-', linewidth=2)
plt.plot(l_2, y_teoricos_2, label='Valores Teóricos', color='green', linestyle='-', linewidth=2)
# Graficar los valores experimentales
plt.scatter(x_experimentales, y_experimentales, label='Valores Experimentales', color='blue', marker='o')

# Añadir títulos y etiquetas
#plt.title('Valores Teóricos y Experimentales')
plt.xlabel('$\lambda$')
plt.ylabel('K')
plt.legend()

# Mostrar la gráfica
plt.show()
