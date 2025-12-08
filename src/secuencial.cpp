#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>  
#include <chrono>
#include <set>
#include <cmath>

using namespace std;

// Para guardar un valor > umbral
struct Hit {
    string gene;   // nombre del gen
    double value;  // valor de expresión
    int row;       // número de fila de datos
};

int main() {
    // Inicialización de tiempo
    auto start_time = std::chrono::high_resolution_clock::now();

    string filename = "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"; //archivo

    int topkMode = 2; // Modo de uso del algoritmo de topk, 1 para filas y 2 para columnas
    // Tamaño del Top-K
    const int K = 5;

    const double THRESHOLD = 100.0;  // umbral
    int maxLines = 200;              // 200 filas para probar

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: No pude abrir el archivo." << endl;
        return 1;
    }
    
    string header;
    if (!getline(file, header)) {
        cerr << "Error: El archivo está vacío o no tiene encabezado." << endl;
        return 1;
    }

    stringstream ss(header);
    string col;
    vector<string> headerNames;

    while (getline(ss, col, '\t')) {
        headerNames.push_back(col);
    }

    int numColumns = (int)headerNames.size();
    cout << "Columnas detectadas: " << numColumns << endl;

    if (numColumns < 2) {
        cerr << "Error: Se esperaban al menos 2 columnas (gene_id + expresión)." << endl;
        return 1;
    }

    // Crear un Top-K por columna usando multiset (orden ascendente)
    vector<multiset<double>> topK(numColumns);
    // Crear un Top-k para filas
    multiset<pair<double,int>> topKRows;
    //limites de filas para el procesamiento de topK por filas
    long int limMin = 5;
    long int limMax = 50;

    // vectores para calculo de media y desviación estándar
    vector<double> mean(numColumns, 0.0);
    vector<double> sd(numColumns, 0.0);

    // Columna 0: gene_id
    // Columnas +1 ya son numericos
    int numExprCols = numColumns - 1;

    // Conteo de valores > umbral por columna
    vector<long long> countAbove(numExprCols, 0);
    vector<vector<Hit>> hits(numExprCols);

    // vector para almacenar todas las filas numéricas
    vector<vector<double>> data;

    string line;
    int counter = 0;  // cuenta las filas de datos

    while (getline(file, line) && counter < maxLines) {
    //while (getline(file, line) ) {

        if (line.empty())
            continue;  // saltar líneas vacías

        stringstream row(line);
        string cell;
        vector<double> rowData;

        // Leer la primera columna: el nombre del gen
        string geneName;
        if (!getline(row, geneName, '\t')) {
            continue;  
        }

        // Leer las columnas de expresión
        int exprColIndex = 0; 

        while (getline(row, cell, '\t') && exprColIndex < numExprCols) {

            // Omitir NA o celdas vacías
            if (cell == "NA" || cell.empty()) {
                exprColIndex++;
                // Se toman las celdas vacias o con NA como 0
                rowData.push_back(0);
                continue;
            }

            // Convertir a double
            double val = strtod(cell.c_str(), nullptr);
            // Agregar val a vector de datos
            rowData.push_back(val);
            // Verificar contra el umbral
            if (val > THRESHOLD) {
                countAbove[exprColIndex]++;

                Hit h;
                h.gene = geneName;
                h.value = val;
                h.row = counter;  // número de fila de datos

                hits[exprColIndex].push_back(h);
            }

            // Implementación de topk por filas
            if (topkMode==1){
                // Referencia al Top-K de la columna
                auto& tk = topK[exprColIndex];

                // Si el Top-K tiene menos de K valores → insertar
                if (tk.size() < K) {
                    tk.insert(val);
                }
                // Si está lleno, verificamos si reemplazamos
                else {
                    auto smallest = tk.begin();  // menor valor actual
                    if (val > *smallest) {
                        tk.erase(smallest);      // eliminar menor
                        tk.insert(val);          // agregar el nuevo
                    }
                }
            }
            // Implementación de topk por rango de filas
            else{
                if ( limMin < exprColIndex && limMax < exprColIndex ){
                    auto& tk = topKRows;
                    if ((int)tk.size() < K) {
                        tk.insert({val, exprColIndex});
                    } else {
                        auto smallest = tk.begin();
                        if (val > smallest->first) {
                            tk.erase(smallest);
                            tk.insert({val, exprColIndex});
                        }
                    }
                }
            }

            mean[exprColIndex] += val;

            exprColIndex++;
            // Se llega al final de la linea
            if(numExprCols <= exprColIndex){
                data.push_back(rowData);
            }
        }

        counter++;
    }

    file.close();
    // counter != data.size()
    // filas data.size();
    // columnas data[1].size() es una matriz cuadrada
    
    // cálculo de medias
    for (auto i=0; i < data[0].size(); i++){
        mean[i] = mean[i]/data.size();
    }

    // Calcular desviaciones estándar
    for (int c = 0; c < data[1].size(); c++) {
        double sumSq = 0.0;
        for (int r = 0; r < data.size()-1; r++) {
            sumSq += pow(data[r][c] - mean[c], 2);
        }
        sd[c] = sqrt(sumSq / data[1].size());
        if (sd[c] == 0) sd[c] = 1;   // evitar división por cero
    }

    // Generar matriz de Z-SCORE
    vector<vector<double>> zdata(data.size(), vector<double>(numExprCols));

    for (int r = 0; r < data.size(); r++) {
        for (int c = 0; c < data[1].size(); c++) {
            zdata[r][c] = (data[r][c] - mean[c]) / sd[c];
        }
    }

    ofstream out("results.csv");

    if (!out.is_open()) {
        cerr << "Error: No pude crear el archivo CSV." << endl;
        return 1;
    }

    out << "Filas de datos procesadas: " << counter << endl;
    out << "Umbral usado: " << THRESHOLD << endl;

    out << "\n=== RESUMEN: VALORES > UMBRAL POR COLUMNA ===\n\n";

    for (int i = 0; i < numExprCols; i++) {
        int colReal = i + 1;  // porque la 0 es gene_id
        out << "Columna " << colReal
             << " (" << headerNames[colReal] << "): "
             << countAbove[i] << " valores > " << THRESHOLD << "\n";
    }

    if (topkMode == 1){
        out << "\n=== TOP K por columna ===\n\n";

        for (int i = 0; i < numColumns; i++) {
            out << "Columna " << i << " (" << headerNames[i] << "): ";
            for (double v : topK[i])
            out << v << " ";
            out << "\n";
        }
    }
    else{
        out << "\n=== TOP K agrupados por filas ===\n\n";
        for (int i = 0; i < counter; i++) {
            out << "Fila " << i << ": ";
            for (auto &p : topKRows)
                out << p.first << ":" << headerNames[p.second] << "  ";
            out << endl;
        }
    }
    // escribir datos estandarizados
    for (int r = 0; r < data.size(); r++) {
        for (int c = 0; c < data[1].size(); c++) {
            out << zdata[r][c];
            if (c < data[1].size()) 
                out << "\t";
        }
        out << "\n";
    }
    out.close();

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> serial_duration = end_time - start_time;
    cout << "Serial time: " << chrono::duration<double>(serial_duration).count()  << std::endl;
    return 0;
}
