# Checklist de preparación antes del parcial (con internet disponible)

Hacer esto durante el tiempo previo en la compu del laboratorio, mientras todavía hay
internet — así en el parcial (sin internet) ya está todo resuelto.

## 1. Compilador de C++ (`g++`)

```bash
g++ --version
```

Si no está instalado:

```bash
sudo apt update && sudo apt install -y g++
```

## 2. Python 3 + matplotlib (para graficar.py)

```bash
python3 -c "import matplotlib; print(matplotlib.__version__)"
```

Si tira error:

```bash
sudo apt install -y python3-matplotlib
```

Si no hay permisos de `sudo` en esa máquina (alternativa sin administrador):

```bash
pip3 install --user matplotlib
```

## 3. GeoGebra portable

Ya está armado y probado (instrucciones completas en el apunte, sección "GeoGebra portable (Linux)").
Llevarlo en el pendrive — pesa ~111MB, puede tardar en copiarse, hacerlo con tiempo.

```bash
unzip GeoGebra-Linux64-Portable-6-0-804-0.zip
cd GeoGebra-linux-x64
chmod +x GeoGebra
./GeoGebra
# si tira error de "sandbox":
./GeoGebra --no-sandbox
```

## 4. Copiar todos los archivos

Al pendrive (o directo al escritorio de esa máquina, si el tiempo previo lo permite):

- Todas las carpetas con los `.cpp` (LocalizacionRaices, Matrices, Interpolacion,
  AjusteDeCurvas, Integracion, EcuacionesDiferenciales)
- `apunte-metodos-numericos.html`
- `EcuacionesDiferenciales/graficar.py`
- El `.zip` de GeoGebra portable

## 5. Probar todo una vez, ahí mismo, en esa máquina real

No confiar en que "debería andar" — probar de verdad, mientras todavía hay internet
para solucionar cualquier cosa rara:

```bash
g++ Integracion/trapecio.cpp -o trapecio && ./trapecio    # compila y corre algo
python3 EcuacionesDiferenciales/graficar.py                # con un CSV de prueba cualquiera
```

Esto descubre cualquier diferencia de esa instalación puntual (versión de Ubuntu,
falta algún paquete, etc.) mientras se puede arreglar — no en medio del parcial.

**Nota:** `graficar.py` además de abrir la ventana interactiva, siempre guarda un PNG
del gráfico — así que aunque esa máquina no tenga bien configurada la parte gráfica
interactiva (pasa en instalaciones mínimas de Ubuntu), el archivo con el gráfico va a
estar igual, abrible con cualquier visor de imágenes.
