#include <iostream>
#include <random>
#include <ctime>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "include/BasicCDS.h"

using namespace std;
using namespace cds;

#define PRINT 1
#define TEST 1

uint bM; 		// bits for MAX

// Structure with all globals parameters program
typedef struct {
	ulong *A;
	ulong *X;
	ulong n;
	ulong MAX;

	ulong sizeA, sizeX;
	ulong nWX;		// number of Words for X[]
} ParProg;

void genArrays(ParProg *par);

void quickSort(ulong *A, int l, int r);
void quickSortX(ulong *X, int l, int r);

/*
Aveces el entero se lee una sola vez y luego se trabaja con una variavble int o long int, por tanto el tiempo que se pierde es despreciable respecto a lo que se gana en espacio cuando los enteros son relativamente pequeños.
- probar utilizar 5 bytes para numeros grandes
*/

int main(int argc, char** argv){
	/*int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	
	int nCells = n/64;

	ulong **X = new ulong*[m];
	for(int i=0; i<m; i++){
		x[i] = new ulong[nCells];
		for (int j=0; j<nCells; j++)
			X[i][j] = 0;
	}
*/



	ParProg *par = new ParProg();
	int i,j;
	ulong k;

	if(argc != 3){
		cout << "Execution Error! call: ./readWriteInt <n> <MAX>" << endl;
		exit(EXIT_FAILURE);
	}
	par->n = atoi(argv[1]);
	par->MAX = atoi(argv[2]);

	cout << "Parameters..." << endl;
	cout << "n = " << par->n << endl;
	cout << "MAX = " << par->MAX << endl;

	genArrays(par);
	
	quickSort(par->A, 0, par->n-1);
	cout << "A[] = ";
	for (i=0; i<par->n; i++)
		cout << par->A[i] << " ";
	cout << endl;
	
	quickSortX(par->X, 0, par->n-1);
	cout << "X[] = ";
	//read values from X[]...
	for (i=j=0; i<par->n; i++, j+=bM){
		cout << getNum64(par->X, j, bM) << " ";
	}
	cout << endl;
	
	if(TEST){
		// check all the position A[i] == X[i]...
		for (i=j=0; i<par->n; i++, j+=bM){
			k = getNum64(par->X, j, bM);
			if (k!=par->A[i]){
				cout << "ERROR. A[" <<i<< "] = " << par->A[i] << " != X[i] = " << k << endl;
				exit(-1);
			}
		}

		cout << "Test OK !!" << endl;
	}

	cout << "################## " << endl;
	return 0;
}


void genArrays(ParProg *par){
	ulong i, j, k;

	par->A = new ulong[par->n];
	par->sizeA = par->n*sizeof(ulong);
	for (i=0; i<par->n; i++)
		par->A[i] = rand()%par->MAX;

	// **************************************************************************
	// create X[] array...
	bM = 1+log2(par->MAX);
	par->nWX = (par->n*bM)/(sizeof(ulong)*8); // number of words for X
	if ((par->n*bM)%(sizeof(ulong)*8)>0)
		par->nWX++;

	par->X = new ulong[par->nWX];
	par->sizeX = par->nWX*sizeof(ulong);

	cout << "bM = " << bM << endl;
	cout << "nWX = " << par->nWX << endl;
	cout << " size for A[] = " << par->sizeA/(1024.0*1024.0) << " MiB (using ulong per cell)" << endl;
	cout << " size for X[] = " << par->sizeX/(1024.0*1024.0) << " MiB" << endl;

	// store values from A[] into X[] (calling the method setNum64())...
	for (i=j=0; i<par->n; i++, j+=bM)
		setNum64(par->X, j, bM, par->A[i]);

	if (PRINT){
		cout << "A[] = ";
		for (i=0; i<par->n; i++)
			cout << par->A[i] << " ";
		cout << endl;

		cout << "X[] = ";
		//read values from X[]...
		for (i=j=0; i<par->n; i++, j+=bM){
			cout << getNum64(par->X, j, bM) << " ";
		}
		cout << endl;

		// print bits using printBitsUlong()...
		for (i=0; i<par->nWX; i++){
			printBitsUlong(par->X[i]);
			cout << " - ";
		}
		cout << endl;
	}

	if(TEST){
		// check all the position A[i] == X[i]...
		for (i=j=0; i<par->n; i++, j+=bM){
			k = getNum64(par->X, j, bM);
			if (k!=par->A[i]){
				cout << "ERROR. A[" <<i<< "] = " << par->A[i] << " != X[i] = " << k << endl;
				exit(-1);
			}
		}

		cout << "Test OK !!" << endl;
	}
}


// realiza la particion de A[l..r] retornando la posision de la particion t,
// dejando A tal que todo elemento en A[l..t-1] < A[t] <= A[t+1...r]
int partition(ulong *A, int l, int r){
	int i, p;
	ulong pv;

	p = l;
	pv = A[p];
	for(i=l; i<=r; i++){
		if (A[i] < pv){
			swap(A[i], A[p+1]);
			p++;
		}
	}
	swap(A[l], A[p]);

	return p;
}

// ordena los elementos de A con el algoritmo quickSort
void quickSort(ulong *A, int l, int r){
	if (l<r){
		int p = partition(A, l, r);
		quickSort(A, l, p-1);
		quickSort(A, p+1, r);
	}
}



int partitionX(ulong *X, int l, int r){
	int i, p;
	ulong pv, u, v;

	p = l;
	pv = getNum64(X, p*bM, bM); // A[p];
	for(i=l; i<=r; i++){
		u = getNum64(X, i*bM, bM); // A[i];
		if (u < pv){
			v = getNum64(X, (p+1)*bM, bM); // A[p+1];
				
			//swap(A[i], A[p+1]);
			setNum64(X, i*bM, bM, v);
			setNum64(X, (p+1)*bM, bM, u);
			
			p++;
		}
	}
	
	u = getNum64(X, l*bM, bM);
	v = getNum64(X, p*bM, bM);
	//swap(A[l], A[p]);
	setNum64(X, l*bM, bM, v);
	setNum64(X, p*bM, bM, u);

	return p;
}

void quickSortX(ulong *X, int l, int r){
	if (l<r){
		int p = partitionX(X, l, r);
		quickSortX(X, l, p-1);
		quickSortX(X, p+1, r);
	}
}




