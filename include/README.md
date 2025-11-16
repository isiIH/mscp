# Explicación del código

La arquitectura del código está modularizada para separar la representación de datos, las estructuras de soporte, el núcleo del algoritmo GRASP y las rutinas de segmentación.

## Representación de datos de la instancia

1. Set:
   - Encargada de la representación y manipulación de un subconjunto individual así como también del universo de elementos, utilizando de base el archivo BasicCDS para la gestión de la secuencia de bits.
   - Cuando se inicializa un subconjunto, automáticamente se setean todos los bits a 0.
   - Soporta la adición, substracción y checkeo de elementos individuales, y permite obtener el número total de elementos activos dentro del bitstring. Además, permite realizar operaciones de conjuntos como la unión, substracción e intersección entre diferentes objetos de la clase Set.

2. SCP:
   - Contiene información de la instancia original del problema.
   - Se encarga de la lectura y el análisis del problema para la representación sucinta del universo y subconjuntos.
   - Para la lectura, se utiliza la función `readFileScp` o `readFilePartition` dependiendo del tipo de instancia que se utilice, ya que se leen de distinta forma.

| Método | Descripción |
| :--- | :--- |
| `readFileScp` | Lee archivos pertenecientes a los **[Set Covering Problems (SCP)](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)**. <br>**Formato de la instancia:** <br>1. Número de filas (m), número de columnas (n). <br>2. El costo de cada columna c(j), j=1,...,n. <br>3. Por cada fila i (i=1,...,m): el número de columnas que cubren la fila i, seguido de la lista de las columnas que la cubren. |
| `readFilePartition` | Lee archivos pertenecientes a los **[Set Partitioning Problems (SPP)](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/sppinfo.html)**. <br>**Formato de la instancia:** <br>1. Número de filas, número de columnas (n). <br>2. Para cada columna j (j=1,...,n): <br>&nbsp;&nbsp;&nbsp;&nbsp;a) Costo de la columna. <br>&nbsp;&nbsp;&nbsp;&nbsp;b) Número de filas cubiertas por j. <br>&nbsp;&nbsp;&nbsp;&nbsp;c) Lista de las filas cubiertas por j. |

3. RowCovering:
   - Implementa la estructura RowMap descrito en el artículo de tesis.
   - Agrega el método `countIntersection` utilizado para crear la lista de relaciones o `Edges` utilizado en la parte de segmentación. Con esto se puede saber el número de subconjuntos que contienen a dos elementos.

4. SetCover:
   - Representa a la solución de una instancia específica. Entre sus atributos se encuentran el universo del problema o subproblema específico `X`, el universo residual `U`, el rowMap asociado al universo, una lista con los subconjuntos considerados únicos y un objeto Set con los subconjuntos excluidos de la solución.
   - Cuando se inicializa un objeto, automáticamente se inicia el preprocesamiento. En este método se realiza los siguientes pasos:
     1. Ejecución de la dominación de subconjuntos.
     2. Creación del RowMap del universo.
     3. Ejecución de la inclusión de subconjuntos.
   - Otras funciones desarrolladas son la eliminación de subconjuntos de la solución, la obtención de un Set con la unión de todos los subconjuntos dentro de la solución, el chequeo de si se forma un set cover con los subconjuntos actuales, la obtención del número de subconjuntos de la solución, y la función para remover subconjuntos redundantes.

## Estructuras para el algoritmo principal y segmentación
