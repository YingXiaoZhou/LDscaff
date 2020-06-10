#ifndef CAL_H
#define CAL_H

#include <stdio.h>
#include <stdlib.h>

using namespace std;

float mean(int arr[], int size);
int filt(int arr1[], int arr2[], int* farr1, int* farr2, int size);
int balance(int arr1[], int arr2[], int* barr1, int* barr2, int size);
float corr(int arr1[], int arr2[], int size);
float fcorr(int arr1[], int arr2[], int size);
float bcorr(int arr1[], int arr2[], int size);
// float LD(int arr1[], int arr2[], int size);

#endif

