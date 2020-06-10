/* Author: Zachary_ZHAO */
#include <algorithm>
#include <functional>
#include <cmath>

#include "tools.h"

using namespace std;


// log trick
float factorial(int n) {
    if (n == 1 || n == 0)
        return 0;
    // return factorial(n-1) * (long long int)n;
    return factorial(n-1) + log(n);
}


// log trick
float combine(int n, int m) {
    // return factorial(n) / (factorial(m) * factorial(n-m));
    return factorial(n) - (factorial(m) + factorial(n-m));
}


string& lTrim(string &ss){
    string::iterator iter = find_if(ss.begin(), ss.end(), not1(ptr_fun<int, int>(isspace)));
    ss.erase(ss.begin(), iter);
    return ss;
}


string& rTrim(string &ss){
    string::reverse_iterator iter = find_if(ss.rbegin(), ss.rend(), not1(ptr_fun<int, int>(isspace)));
    ss.erase(iter.base(), ss.end());
    return ss;
}


string& trim(string &ss){
    lTrim(rTrim(ss));
    return ss;
}


void parseLine(string strLine, string deLimiter, vector<string> &vecString){
    size_t pos = 0;
    string tmp;
    vecString.clear();
    // cout << "pos = " << strLine.find("\n") << endl;
    while ((pos = strLine.find(deLimiter)) != string::npos){
        tmp = strLine.substr(0, pos);
        vecString.push_back(trim(tmp));
        strLine.erase(0, pos + deLimiter.length());
    }
    if (strLine != "") vecString.push_back(trim(strLine));
}


float cal_pe(float ref_freq, float alt_freq) {
    float pe = 0;
    pe += ref_freq * pow((1 - ref_freq), 2);
    pe += alt_freq * pow((1 - alt_freq), 2);
    pe -= pow(ref_freq, 2) * pow(alt_freq, 2) * (4 - 3 * ref_freq - 3 * alt_freq);
    return pe;
}


int get_index(string name, string array[], int size) {
    for (int i = 0; i != size; i ++) {
        if (array[i] == name)
            return i;
    }
    return -1;
}

