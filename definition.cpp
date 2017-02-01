#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "header.h"

// ALOKACJA PAMIÊCI DLA MACIERZY
float** allocMatrixEX(int height, int width){
	float** matrix = new float*[height];
	for (int i = 0; i < height; ++i)
    	matrix[i] = new float[width];
    return matrix;
}

// USUWANIE PAMIÊCI PO MACIERZY
void deleteMatrixEX(float**& matrix, int height){
    for (int i = 0; i < height; ++i)
    		delete matrix[i];
	delete matrix;
}

// GENEROWANIE MACIERZY
void generateMX(float**& matrix, int size){
	float number;
	cout << "Podaj wartosci macierzy: ";
	cin >> number;
	cout << endl;
	for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            matrix[i][j] = number;
}

//LICZENIE ŒREDNIEJ WARTOSCI MACIERZY
float midvalue(float**& MC, int size){
	float mid_v = 0;
	for (int i=0; i<size; i++){
        for (int j=0; j<size; j++)
            mid_v+=MC[i][j];
       // cout << "\nsuma komorek: " << mid_v;
	}
    //cout << "\nkomorka 100,50: " << MC[99][49] << endl;
   // cout << "komorka 50,158: " << MC[49][158] << endl;
	return mid_v/(size*size);
}

// KOPIOWANIE MACIERZY
void copyMX(float**& matrix, float**& MX0, int size){
	for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            matrix[i][j] = MX0[i][j];
}

// MNO¯ENIE MACIERZY
float** multiply(float**& MX0, float**& MX1, float**& MX2, int size){
	/* for(int y=0; y<size; y++){
		for(int j=0; j<size; j++){
			for(int i=0; i<size; i++){
				MX0[y][j]+=MX1[y][i]*MX2[i][j];
			}
		}
		//cout << "mnozenie kolumny " << y << " macierzy " << MX0 << endl;
	}*/
	for(int i=0; i<size; i++){ // przesuniêcie obliczania o jeden wiersz w dó³ macierzyA
		for(int j=0; j<size; j++){ // przesuniêcie obliczania o jedn¹ kolumnê w prawo macierzyB
			MX0[i][j]=0;
			for(int k=0; k<size; k++){ // pêtla mno¿enia wierszaA i kolumnyB
				MX0[i][j]+=MX1[i][k]*MX2[k][j]; // mno¿enie
			}
		}
	}
	
	return MX0;
}

// DODAWANIE MACIERZY
float** add(float**& MX0, float**& MX1, float**& MX2, int size){
	for (int i=0; i<size/2; i++){
        for (int j=0; j<size/2; j++)
            MX0[i][j] = MX1[i][j]+MX2[i][j];
    	//cout << "dodawanie kolumny " << i << " macierzy " << MX0 << endl;
    }
	return MX0;
}

// ODEJMOWANIE MACIERZY
float** subtract(float**& MX0, float**& MX1, float**& MX2, int size){
	for (int i=0; i<size/2; i++){
        for (int j=0; j<size/2; j++)
            MX0[i][j] = MX1[i][j]-MX2[i][j];
    	//cout << "odejmowanie kolumny " << i << " macierzy " << MX0 << endl;
    }
	return MX0;
}

// PODZIA£ G£OWNYCH MACIERZY NA ÆWIARTKI
void subdivide(float**MX0, float**& matrix, int i, int j, int size){
	int x=i;
	for (int a=0; a<size; a++)
	{
		int y=j;
		for (int b=0; b<size; b++)
		{
        	MX0[a][b] = matrix[x][y];
        	y++;
		}
		x++;
	}
}

// KOÑCOWE SUMOWANIE MACIERZY WYNIKOWEJ
void summary(float**& MOC, float**& c11, float**& c12, float**& c21, float**& c22, int size){
	int x=0;
	for (int a=0; a<size/2; a++)
	{
		int y=0;
		for (int b=0; b<size/2; b++)
		{
        	MOC[a][b] = c11[x][y];
        	y++;
		}
		x++;
	}
	x=0;
	for (int a=0; a<size/2; a++)
	{
		int y=0;
		for (int b=size/2; b<size; b++)
		{
        	MOC[a][b] = c12[x][y];
        	y++;
		}
		x++;
	}
	x=0;
	for (int a=size/2; a<size; a++)
	{
		int y=0;
		for (int b=0; b<size/2; b++)
		{
        	MOC[a][b] = c21[x][y];
        	y++;
		}
		x++;
	}
	x=0;
	for (int a=size/2; a<size; a++)
	{
		int y=0;
		for (int b=size/2; b<size; b++)
		{
        	MOC[a][b] = c22[x][y];
        	y++;
		}
		x++;
	}
}

// WYŒWIETLANIE MACIERZY
void showM(float**& MX0, int size){
	for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
        	cout << MX0[i][j] << " | ";
		}
    	cout << endl;
	}
	cout << "**********************************************\n";
}

// FUNKCJA STRASSEN

float** strassen(float**& MXC, float**& MXA, float**& MXB, int value){
int size = value;
    
    /* ALOKACJA PAMIÊCI NA POTRZEBY MACIERZY OPERACYJNYCH */
    float** matrixA = allocMatrixEX(size, size);
    float** matrixB = allocMatrixEX(size, size);
    float** matrixC = allocMatrixEX(size, size);
    
    float** A11 = allocMatrixEX(size/2, size/2);
    float** A12 = allocMatrixEX(size/2, size/2);
    float** A21 = allocMatrixEX(size/2, size/2);
    float** A22 = allocMatrixEX(size/2, size/2);
    float** B11 = allocMatrixEX(size/2, size/2);
    float** B12 = allocMatrixEX(size/2, size/2);
    float** B21 = allocMatrixEX(size/2, size/2);
    float** B22 = allocMatrixEX(size/2, size/2);
    
    float** M1 = allocMatrixEX(size/2,size/2);
    float** M2 = allocMatrixEX(size/2,size/2);
    float** M3 = allocMatrixEX(size/2,size/2);
    float** M4 = allocMatrixEX(size/2,size/2);
    float** M5 = allocMatrixEX(size/2,size/2);
    float** M6 = allocMatrixEX(size/2,size/2);
    float** M7 = allocMatrixEX(size/2,size/2);
    /* KONIEC ALOKACJI */
    
    /* PRZYPISANIE WARTOŒCI DO MACIERZY I ICH PODZIA£ NA ÆWIARTKI */
    
    copyMX(matrixA, MXA, size);
    	subdivide(A11, matrixA, 0, 0, size/2);
		subdivide(A12, matrixA, 0, size/2, size/2);
		subdivide(A21, matrixA, size/2, 0, size/2);
		subdivide(A22, matrixA, size/2, size/2, size/2);
		
		//showM(A11, size/2);
    	//showM(A12, size/2);
    	//showM(A21, size/2);
    	//showM(A22, size/2);
    //	system("pause");
	deleteMatrixEX(matrixA, size);
    
    copyMX(matrixB, MXB, size);
    	subdivide(B11, matrixB, 0, 0, size/2);
		subdivide(B12, matrixB, 0, size/2, size/2);
		subdivide(B21, matrixB, size/2, 0, size/2);
		subdivide(B22, matrixB, size/2, size/2, size/2);
		
		//showM(B11, size/2);
    	//showM(B12, size/2);
    	//showM(B21, size/2);
    	//showM(B22, size/2);
	//	system("pause");
    deleteMatrixEX(matrixB, size);		
	/* KONIEC PRZYPISYWANIA I PODZIA£U MACIERZY WEJŒCIOWYCH */
    
    /* OBLICZANIE MACIERZY OPERACYJNYCH */
    // jeœli ma dzieliæ dalej to mno¿enie odbywa siê za pomoc¹ funkcji strassen(), w innym wypadku za pomoc¹ funkcji multiply() //
    
    float** temp0 = allocMatrixEX(size/2,size/2);
    float** temp1 = allocMatrixEX(size/2,size/2);
    
    if(size>512){ /// HOW DEEP!!!!!!!!!!!!
    	//krok 1
		temp0 = add(temp0, A11, A22, size);
		temp1 = add(temp1, B11, B22, size);
			M1=strassen(M1, temp0, temp1, size/2);
		//krok 2
		temp0 = add(temp0, A21, A22, size);
			M2=strassen(M2, temp0, B11, size/2);
		//krok 3
		temp1 = subtract(temp1, B12, B22, size);
			M3=strassen(M3, A11, temp1, size/2);
		//krok 4
		temp1 = subtract(temp1, B21, B11, size);
			M4=strassen(M4, A22, temp1, size/2);
		//krok 5
		temp1 = add(temp1, A11, A12, size);
			M5=strassen(M5, temp1, B22, size/2);
		//krok 6
		temp0 = subtract(temp0, A21, A11, size);
		temp1 = add(temp1, B11, B12, size);
			M6=strassen(M6, temp0, temp1, size/2);
		//krok 7
		temp0 = subtract(temp0, A12, A22, size);
		temp1 = add(temp1, B21, B22, size);
			M7=strassen(M7, temp0, temp1, size/2);
	} else{
    	//krok 1
		temp0 = add(temp0, A11, A22, size);
		temp1 = add(temp1, B11, B22, size);
			M1=multiply(M1, temp0, temp1, size/2);
		//krok 2
		temp0 = add(temp0, A21, A22, size);
			M2=multiply(M2, temp0, B11, size/2);
		//krok 3
		temp1 = subtract(temp1, B12, B22, size);
			M3=multiply(M3, A11, temp1, size/2);
		//krok 4
		temp1 = subtract(temp1, B21, B11, size);
			M4=multiply(M4, A22, temp1, size/2);
		//krok 5
		temp1 = add(temp1, A11, A12, size);
			M5=multiply(M5, temp1, B22, size/2);
		//krok 6
		temp0 = subtract(temp0, A21, A11, size);
		temp1 = add(temp1, B11, B12, size);
			M6=multiply(M6, temp0, temp1, size/2);
		//krok 7
		temp0 = subtract(temp0, A12, A22, size);
		temp1 = add(temp1, B21, B22, size);
			M7=multiply(M7, temp0, temp1, size/2);
	}
	
	/* KONIEC OBLICZEÑ OPERACYJNYCH */
    
    deleteMatrixEX(A11, size/2);
    deleteMatrixEX(A12, size/2);
    deleteMatrixEX(A21, size/2);
    deleteMatrixEX(A22, size/2);
    deleteMatrixEX(B11, size/2);
    deleteMatrixEX(B12, size/2);
    deleteMatrixEX(B21, size/2);
    deleteMatrixEX(B22, size/2);
    deleteMatrixEX(temp0, size/2);
    deleteMatrixEX(temp1, size/2);
    
    float** C11 = allocMatrixEX(size/2,size/2);
    float** C12 = allocMatrixEX(size/2,size/2);
    float** C21 = allocMatrixEX(size/2,size/2);
    float** C22 = allocMatrixEX(size/2,size/2);
    float** temp2 = allocMatrixEX(size/2,size/2);
    float** temp3 = allocMatrixEX(size/2,size/2);
    
	/* OSTATECZNE SUMOWANIE WARTOŒCI WYJŒCIOWYCH */
	temp2 = add(temp2, M1, M4, size);
	temp3 = subtract(temp3, temp2, M5, size);
	C11=add(C11, temp3, M7, size);
	C12= add(C12, M3, M5, size);
	C21= add(C21, M2, M4, size);
	temp2 = subtract(temp2, M1, M2, size);
	temp3 = add(temp3, temp2, M3, size);
	C22=add(C22, temp3, M6, size);
	/* KONIEC SUMOWANIA */
	
	deleteMatrixEX(M1, size/2);
    deleteMatrixEX(M2, size/2);
    deleteMatrixEX(M3, size/2);
    deleteMatrixEX(M4, size/2);
    deleteMatrixEX(M5, size/2);
    deleteMatrixEX(M6, size/2);
    deleteMatrixEX(M7, size/2);
    deleteMatrixEX(temp2, size/2);
    deleteMatrixEX(temp3, size/2);
    
	/* SK£ADANIE MACIERZY WYJŒCIOWEJ */
	summary(MXC, C11, C12, C21, C22, size);
    
    deleteMatrixEX(C11, size/2);
    deleteMatrixEX(C12, size/2);
    deleteMatrixEX(C21, size/2);
    deleteMatrixEX(C22, size/2);
    
    //cout << ".";
	//showM(MXC, size);
	copyMX(matrixC, MXC, size);
	//deleteMatrixEX(MXC, size);
    return matrixC;
}
