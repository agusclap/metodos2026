import csv
import matplotlib.pyplot as plt

# Grafica cualquier CSV generado por edo.cpp (o por trapecio.cpp/simpson.cpp/
# lagrangeGrilla.cpp, que usan el mismo formato: primera columna "x", el
# resto son curvas a dibujar). Si hay una columna "error", no se grafica
# junto con las demas (otra escala) -- se imprime el maximo nomas.

nombre = input("Nombre del archivo CSV a graficar: ").strip()

x = []
columnas = {}

with open(nombre, newline="") as f:
    lector = csv.reader(f)
    encabezados = next(lector)
    for col in encabezados[1:]:
        columnas[col] = []
    for fila in lector:
        x.append(float(fila[0]))
        for i, col in enumerate(encabezados[1:], start=1):
            columnas[col].append(float(fila[i]))

plt.figure()
for nombre_col, valores in columnas.items():
    if nombre_col == "error":
        print(f"Error maximo ({nombre_col}): {max(valores):.6g}")
        continue
    estilo = "--" if nombre_col == "exacta" else "o-"
    plt.plot(x, valores, estilo, label=nombre_col, markersize=4)

plt.xlabel("x")
plt.ylabel("y")
plt.title(nombre)
plt.legend()
plt.grid(True)

salida = nombre.rsplit(".", 1)[0] + ".png"
plt.savefig(salida)
print(f"Grafico guardado en {salida}")
plt.show()
