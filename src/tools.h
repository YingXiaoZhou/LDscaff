/* Author: Zachary_ZHAO */

#ifndef TOOLS_H
#define TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "math.h"

using namespace std;

float factorial(int n);
float combine(int n, int m);

float cal_pe(float ref_freq, float alt_freq);
int get_index(string name, string array[], int size);

string& rTrim(string &ss);
string& lTrim(string &ss);
string& trim(string &ss);
void parseLine(string strLine, string deLimiter, vector<string> &vecString);

#endif
