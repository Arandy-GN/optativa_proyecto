#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <omp.h> 

using namespace std;

//estructuras
struct Hit {
    string gene;
    string patient;
    double value;
};

struct TopKItem {
    double value;
    int geneIndex;
};

struct Dataset {
    int rows = 0;
    int cols = 0;
    vector<string> colNames;
    vector<string> geneNames;
    // optimizacion de memoria
    // Usamos un vector unidimensional en lugar de vector<vector>.
    vector<double> data; 
};

// funcion para cargar el archivo 
double safe_stod(string s) {
    if (s.empty()) return 0.0;
    if (s.back() == '\r') s.pop_back();
    if (s == "NA" || s == "na" || s == "null" || s == "NaN") return 0.0;
    try { return stod(s); } catch (...) { return 0.0; }
}

Dataset loadCSV(const string& filename, int maxLines) {
    Dataset db;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error FATAL: No se pudo abrir " << filename << endl;
        exit(1);
    }
    string line;
    if (getline(file, line)) { // Headers
        if (!line.empty() && line.back() == '\r') line.pop_back();
        stringstream ss(line);
        string cell;
        getline(ss, cell, '\t'); 
        while (getline(ss, cell, '\t')) db.colNames.push_back(cell);
    }
    db.cols = db.colNames.size();

    cout << "Cargando datos a RAM..." << endl;
    while (getline(file, line) && (maxLines == -1 || db.rows < maxLines)) {
        if (line.empty()) continue;
        size_t tabPos = line.find('\t');
        if (tabPos == string::npos) continue;
        db.geneNames.push_back(line.substr(0, tabPos));
        size_t start = tabPos + 1;
        size_t end;
        while ((end = line.find('\t', start)) != string::npos) {
            db.data.push_back(safe_stod(line.substr(start, end - start)));
            start = end + 1;
        }
        db.data.push_back(safe_stod(line.substr(start)));
        db.rows++;
    }
    file.close();
    return db;
}

//procesamiento paralelo
void process_parallel(Dataset& db, double threshold, vector<Hit>& global_hits, vector<vector<TopKItem>>& all_topk, int num_threads) {
    
    omp_set_num_threads(num_threads);

    // Memoria auxiliar compartida
    vector<double> means(db.cols, 0.0);
    vector<double> std_devs(db.cols, 0.0);
    const int K = 10;

    
    #pragma omp parallel
    {
        // CALCULAR MEDIAS Y DESVIACIONES
        // Usamos 'schedule(static)' porque la carga de trabajo es idéntica
        // para cada columna (mismo número de sumas), así que el reparto fijo es óptimo.
        #pragma omp for schedule(static)
        for (int j = 0; j < db.cols; j++) {
            double sum = 0.0;
            for (int i = 0; i < db.rows; i++) sum += db.data[i * db.cols + j];
            means[j] = sum / db.rows;

            double sum_sq = 0.0;
            for (int i = 0; i < db.rows; i++) {
                double val = db.data[i * db.cols + j];
                sum_sq += pow(val - means[j], 2);
            }
            double variance = (db.rows > 1) ? (sum_sq / (db.rows - 1)) : 0.0;
            std_devs[j] = sqrt(variance);
            if (std_devs[j] == 0.0) std_devs[j] = 1.0; 
        }

        // Z-SCORE y BÚSQUEDA
        // Variable PRIVADA por hilo. Evita 'Race Conditions' al escribir hallazgos.
        vector<Hit> local_hits; 
        
        // Usamos 'schedule(dynamic)' porque el 'if (val > threshold)' introduce irregularidad
        
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < db.rows; i++) {
            for (int j = 0; j < db.cols; j++) {
                int idx = i * db.cols + j;
                double val_original = db.data[idx];

                if (val_original > threshold) {
                    local_hits.push_back({db.geneNames[i], db.colNames[j], val_original});
                }
                // Sobrescribimos en RAM
                db.data[idx] = (val_original - means[j]) / std_devs[j];
            }
        }
        
        // Un solo hilo a la vez vacía su vector local en el global.
        #pragma omp critical
        global_hits.insert(global_hits.end(), local_hits.begin(), local_hits.end());

        // BARRERA DE SINCRONIZACIÓN:
        // Nadie puede empezar el Top-K hasta que TODOS hayan terminado de calcular
        // los Z-Scores en la matriz 'db.data'.
        #pragma omp barrier 
        // TOP-K POR PACIENTE
        #pragma omp for schedule(dynamic)
        for (int j = 0; j < db.cols; j++) {
            
            // Recolectar datos de la columna 'j' en un vector temporal
            vector<TopKItem> col_vals;
            col_vals.reserve(db.rows);
            for (int i = 0; i < db.rows; i++) {
                col_vals.push_back({db.data[i * db.cols + j], i});
            }

            // ALGORITMO: PARTIAL SORT
            // En lugar de ordenar todo el vector (O(N log N)), solo ordenamos
            // los primeros K elementos (O(N log K)). Esto es mucho más rápido para HPC.
            std::partial_sort(col_vals.begin(), 
                              col_vals.begin() + K, 
                              col_vals.end(), 
                              [](const TopKItem& a, const TopKItem& b) {
                                  return a.value > b.value; // Orden Descendente
                              });
            
            
            col_vals.resize(K);
            all_topk[j] = col_vals;
        }
    } 
}

// MAIN 
int main(int argc, char* argv[]) {
    int n_threads = 4; 
    if (argc > 1) n_threads = atoi(argv[1]);

    string filename = "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"; 
    
    // CONFIGURACIÓN
    int maxLines = -1; // -1 = Leer archivo completo
    double umbral = 1000.0;

    // CARGA (Secuencial)
    Dataset db = loadCSV(filename, maxLines);
    cout << "Datos cargados: " << db.rows << " genes x " << db.cols << " pacientes." << endl;

    vector<Hit> hits;
    vector<vector<TopKItem>> all_topk(db.cols);

    // PROCESAMIENTO (Medición de Speedup)
    cout << "Procesando con " << n_threads << " hilos..." << endl;
    
    // Cronómetro solo para la parte paralela 
    auto start = chrono::high_resolution_clock::now();

    process_parallel(db, umbral, hits, all_topk, n_threads);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << ">>> Tiempo de Cómputo: " << duration.count() << " s <<<" << endl;

    // ESCRITURA DE ARCHIVOS (Fuera del cronómetro)
    cout << "Escribiendo search_hits.csv (" << hits.size() << " registros)..." << endl;
    ofstream f_hits("search_hits.csv");
    f_hits << "Gene,Patient,Value\n";
    for(const auto& h : hits) {
        f_hits << h.gene << "," << h.patient << "," << h.value << "\n";
    }
    f_hits.close();

    cout << "Escribiendo topk_results.csv..." << endl;
    ofstream f_topk("topk_results.csv");
    f_topk << "Patient,Rank,Gene,Z_Score\n";
    
    for(int j=0; j < db.cols; j++) {
        for(int k=0; k < all_topk[j].size(); k++) {
            TopKItem item = all_topk[j][k];
            f_topk << db.colNames[j] << "," 
                   << (k+1) << "," 
                   << db.geneNames[item.geneIndex] << "," 
                   << item.value << "\n";
        }
    }
    f_topk.close();

    cout << "Archivos generados exitosamente." << endl;
    return 0;
}