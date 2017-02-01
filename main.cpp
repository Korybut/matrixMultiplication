#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "header.h"

clock_t start, stop;
double czas;

int main(void) {
    srand(time(NULL));
    int value;
    cout << "Podaj wartosc: ";
    cin >> value;
    
    // tu powinna zacz¹æ siê funkcja strassen
    int size = value;
	
	float** MA = allocMatrixEX(size, size);
    float** MB = allocMatrixEX(size, size);
    float** MC = allocMatrixEX(size, size);
    
    generateMX(MA, size);
    generateMX(MB, size);
    //showM(MA, size);
	//showM(MB, size);
    //system("pause");
    
    start = clock();
	MC = strassen(MC, MA, MB, size);
	stop = clock();
	czas = (double)(stop-start)/CLOCKS_PER_SEC;
	cout << "Najszybsze mnozenie macierzy EVER: " << czas << "!";
	cout << "\nSrednia wartosc: " << midvalue(MC, size);
	cout << endl;
//	showM(MC, size);
	
	deleteMatrixEX(MA, size);
	deleteMatrixEX(MB, size);
    deleteMatrixEX(MC, size);
    
    return 0;
}
