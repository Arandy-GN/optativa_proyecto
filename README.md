<<<<<<< HEAD
# Proyecto OpenMP Gen�mica 
=======
# Optimización Paralela de Procesamiento Genómico a Gran Escala con OpenMP
## Integrantes
- Garcia Navarro Arandy 
- Hernández Pérez Fernando Emmanuel
- Martínez Cruz Aldo
- Tenorio Camacho Javier
## Resumen
El proyecto implementa técnicas de cómputo paralelo en C/C++ utilizando OpenMP para acelerar el análisis de datos genómicos masivos del dataset EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv. Se implementa y compara versiones secuenciales y paralelas de tareas esenciales: cálculo de z-scores, búsquedas por umbral, conteos y métricas por ventana configurables. El objetivo es reducir tiempos de procesamiento y mejorar la escalabilidad aprovechando múltiples núcleos. Se mide el rendimiento variando hilos y distintas estrategias de scheduling, evaluando speedup, eficiencia y escalabilidad en la base de datos mencionada anteriormente de más de 10⁸ elementos. El programa permite configurar parámetros, generar resultados y tiempos en archivos CSV.
## Instrucciones para ejecutar los códigos
Dentro del repositorio se encuentran 2 códigos, ambos siendo implementaciones de los mismos algoritmos, donde uno corresponde a la versión paralela y el otro a la versión secuencial.
1. Primero se debe de descargar o clonar este repositorio
2. Para ejecutar cada uno de los archivos cpp primero se debe de descargar la base de datos dentro de la misma carpeta, la base de datos se llama EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp y se descarga del siguiente sitio web https://gdc.cancer.gov/about-data/publications/pancanatlas
3. Posteriormente se debe de crear el ejecutable de cada uno de los archivos cpp. Para cada uno de los archivos el ejecutable se obtiene usando el siguiente código
```
g++ secuencial.cpp -o secuencial.exe
```
```
g++ -fopenmp paralelo.cpp -o paralelo.exe
```
4. Después se debe de ejecutar el siguiente código para correr el ejecutable de cada uno. Para la versión paralelo.exe, el # debe de ser sustituido por un número entre 1 a el mayor número de hilos soportado por el equipo en donde se esté trabajando. Estos se encuentran en la carpeta src.
```
./secuencial.exe
```
```
./paralelo.exe #
```
>>>>>>> 57ae23c25f8760dcd470f912951480b12517c6eb
