#include <iostream>
#include <fstream>
using namespace std;

// Genera el sistema de banda del Problema 4 (pentadiagonal: diagonal=12,
// +-1 offset=-2, +-2 offset=1, b=5) y lo escribe como matriz aumentada
// [A|b] en un archivo de texto, listo para leer con leerSistema().
int main() {
    int n;
    cout << "Ingrese n (tamano del sistema, el del enunciado es 50): ";
    cin >> n;

    string nombreArchivo;
    cout << "Nombre del archivo a generar: ";
    cin >> nombreArchivo;

    ofstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo crear el archivo." << endl;
        return 1;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            double valor = 0;
            if (j == i)          valor = 12;
            else if (j == i - 1) valor = -2;
            else if (j == i - 2) valor = 1;
            else if (j == i + 1) valor = -2;
            else if (j == i + 2) valor = 1;
            archivo << valor << " ";
        }
        archivo << 5 << endl; // b_i = 5
    }

    archivo.close();
    cout << "Listo: " << nombreArchivo << " generado (n=" << n << ")." << endl;
    return 0;
}
