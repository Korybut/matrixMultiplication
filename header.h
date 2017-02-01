#ifndef HEADER_H_
#define HEADER_H_
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;

float** allocMatrixEX(int, int);
void deleteMatrixEX(float**&, int size);

void generateMX(float**& matrix, int size);
void copyMX(float**& matrix, float**& MX0, int size);

void subdivide(float**MX0, float**& matrix, int i, int j, int size);

float** add(float**& MX0, float**& MX1, float**& MX2, int size);
float** subtract(float**& MX0, float**& MX1, float**& MX2, int size);
float** multiply(float**& MX0, float**& MX1, float**& MX2, int size);
float** strassen(float**& MXC, float**& MXA, float**& MXB, int value);
void summary(float**& MOC, float**& C11, float**& C12, float**& C21, float**& C22, int size);
float midvalue(float**& MC, int size);

void showM(float**& MX0, int size);

#endif
