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

1. Grasp:
   - La ejecución del algoritmo propuesto se centraliza en la función `search`. Dependiendo de la configuración seleccionada, esta función ejecutará `searchPerGroup` (para el GRASP con segmentación, Algoritmo 6 de la tesis) o el GRASP tradicional sin segmentación (Algoritmo 2 de la tesis).
   - Al finalizar cualquiera de los dos métodos de búsqueda, los subconjuntos únicos (uniqueSets) que fueron identificados y reservados durante la fase de preprocesamiento se añaden a la solución final.
   - El método `searchPerGroup` funciona de la siguiente manera:
     1. Primero, se realiza la segmentación del universo llamando a la función `findGroups`, utilizando la heurística de Union-Find o la de MST según se haya configurado.
     2. Una vez definidos los subuniversos (grupos), la función `distributeSubsets` se encarga de asignar cada subconjunto de la instancia original al grupo con el que tenga mayor relación.
     3. Cada subproblema se resuelve en forma paralela. Para ello, la información del grupo es copiada a su propia instancia de solución o SetCover.
     4. Cada solución parcial encontrada es unida a la solución global. Si se utilizó la segmentación basada en MST, se ejecuta una verificación final para remover subconjuntos redundantes. Esto es necesario porque la división del universo por MST puede generar solapamiento entre subproblemas, a diferencia de Union-Find.
   - El método `updateSolution` representa el núcleo de la fase de búsqueda local iterativa del algoritmo GRASP, correspondiente a las líneas 10-18 del Algoritmo 2.
   - El método `randSuccinctSC` implementa el Algoritmo 1 (RANDSUCCINCTSC). Recordar que si no se encuentran mejoras en la solución (IMPROVE == False), en lugar de escoger al subconjunto que cubra la mayor cantidad de elementos restantes, el subconjunto se escoge por probabilidad en base a su cobertura para favorecer la exploración y diversidad.

2. Group:
   - Implementa las heurísticas clave para la segmentación del universo, ya sea a través de Union-Find o Maximum Spanning Tree (MST), dependiendo de la configuración utilizada.
   - La heurística Union-Find tiene su clase propia dentro del archivo. Para esta clase se crean los métodos find y unite descritas en el Algoritmo 3 y 4 de la tesis.
   - Por otro lado, la clase Group ejecuta la segmentación de grupos para ambas heurísticas dependiendo de la configuración utilizada. La función principal corresponde a `findGroups`, el cual ejecuta los siguientes pasos:
     1. Primero, se crea el grafo de elementos a partir del método `createGraph`, el cual identifica las relaciones entre elementos a partir del método `countIntersection` descrito en la clase RowCovering.
     2. Los siguientes pasos dependen de la heurística utilizada. Para el caso del MST, se realizan los pasos descritos en el Algoritmo 5 de la tesis. En el caso de Union-Find, por cada relación de cobertura encontrada entre elementos, se realiza la función `unite`. Posteriormente, para definir los universos de cada grupo (componentes disjuntas), se mapea a cada elemento del universo con su representante utilizando la función `find`.
     3. Por último, se procede con el método `distributeSubsets`. Este método también depende de la heurística. Para el MST, la distribución es más compleja debido al posible solapamiento. Primero, se asegura que el universo de cada grupo sea cubierto, asignando a los subconjuntos que cubran más elementos del grupo. Una vez cubiertos todos los elementos, se terminan de agregar los subconjuntos sobrantes en una segunda ronda. Por otro lado, la distribución de Union-Find es más sencilla porque no hay solapamiento entre grupos. Cada subconjunto se agrega al grupo al que pertenece uno de los elementos que cubre.
