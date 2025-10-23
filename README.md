# Aplicación de heurísticas y técnicas de HPC para la segmentación del universo en el Minimum Set Cover Problem

## Instalación
1. Clona el repositorio y dirígete a la rama **grasp**:
    ```bash
    git clone https://github.com/isiIH/mscp.git
    cd mscp
    git switch grasp

## Pasos previos
- Descomprime `test.zip`, que contiene ejemplos extraídos de la OR-Library y datasets generados mediante el archivo `generate_SC_file.py`.
- Crea un archivo llamado `results.txt` en la raíz del proyecto para almacenar los resultados de los experimentos.

## OR-Library
Los datasets ubicados en la carpeta `test/` corresponden a instancias de la **OR-Library**, específicamente de los problemas **Set covering** y **Set partitioning**.

Referencia: 
- [Set covering](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)
- [Set Partitioning](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/sppinfo.html)

## Generar datasets aleatorios

Puedes generar datasets personalizados utilizando el archivo `generate_SC_file.py`. 

Dentro del archivo se encuentra la siguiente sección que puedes modificar:

```python
generate_set_cover_file("nombre_archivo", n=100, m=200, num_groups=2)
```

donde:
- nombre_archivo: Nombre base del archivo que se almacenará en la carpeta `test/`. Destacar que al nombre final se le incluirá automáticamente información sobre el número de columnas y grupos (Ej: nombre_archivo_200_2.txt).
- n: Número de elementos del universo.
- m: Número de subconjuntos.
- num_groups: Controla el número de grupos generados, útil para realizar pruebas de segmentación.

## Archivo de configuración

Todas las variables de configuración se encuentran en el archivo `config.h`.

- PRINT (0 ó 1): Imprime información principal durante la ejecución.
- CHECK (0 ó 1): Imprime información adicional de los pasos intermedios realizados durante la ejecución.
- TEST (0 ó 1): Guarda información sobre el archivo ejecutado en `result.txt`.
- MAX_RM: 0.5 por default. Indica en porcentaje el máximo número de elementos a eliminar dentro de la solución actual.
- MAX_ITER: 300 por default. Indica el máximo número de iteraciones que realiza GRASP en la fase de mejora.
- GROUP_SEG: Ejecutar sin (0) o con (1) segmentación del universo.
- SEG_TYPE: Si GROUP_SEG = 1, selecciona el método de segmentación, ya sea con UNION-FIND (0) o el método con MST (1).

## Ejecución
Después de haber seteado las variables en `config.h`, compila y ejecuta el programa:

```bash
make
./grasp <filename> <nt>
```

donde:
- filename: Nombre base del archivo almacenado en la carpeta `test/`.
- nt: Número de threads a utilizar.

### Ejemplo

```bash
./grasp rail582 32
    