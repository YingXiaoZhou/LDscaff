#include "cal.h"
// #include "ortools/linear_solver/linear_solver.h"
#include <math.h>
#include <iostream>

using namespace std;

float mean(int arr[], int size) {
    if (size == 0) 
        return 0;
    float sum = 0;
    for (int i = 0; i != size; i ++)
        sum += (float)arr[i];
    return sum / size;
}

int filt(int arr1[], int arr2[], int* farr1, int* farr2, int size) { // filter out heterozygous markers
    int total = 0;
    for (int i = 0; i != size; i ++) {
        // if (arr1[i] == 1 || arr2[i] == 1) continue;
        if (arr1[i] < 0 || arr2[i] < 0) continue;
        else {
            farr1[total] = arr1[i];
            farr2[total] = arr2[i];
            total ++;
        }
    }
    return total;
}

int balance(int arr1[], int arr2[], int* barr1, int* barr2, int size) { // balance major allele rate
    int total = 0;

    // cout << "arr1: " << endl;
    // for (int i = 0; i != size; i ++) 
    //     cout << arr1[i] << ",";
    // cout << endl;
    // cout << "arr2: " << endl;
    // for (int i = 0; i != size; i ++) 
    //     cout << arr2[i] << ",";
    // cout << endl;

    // major allele count
    float pA = 0, pB = 0;
    for (int i = 0; i != size; i ++) {
        if (arr1[i] == 0) pA += 2;
        else if (arr1[i] == 1) pA += 1;
        if (arr2[i] == 0) pB += 2;
        else if (arr2[i] == 1) pB += 1;
    }
    pA /= (2 * size);
    pB /= (2 * size);

    cout << "pA: " << pA << endl;
    cout << "pB: " << pB << endl;

    int m1 = pA >= 0.5 ? 0 : 1;
    int m2 = pB >= 0.5 ? 0 : 1;

    pA = pA >= 0.5 ? pA : (1 - pA);
    pB = pB >= 0.5 ? pB : (1 - pB);
        
    int* t1 = new int[size];
    int* t2 = new int[size];

    if (pA == pB) {
        int x1 = (pA - 0.5) * size;
        for (int i = 0; i != size && x1 > 0; i ++) {
            if (arr1[i] == 2 * m1 && arr2[i] == 2 * m2) x1 --;
            else {
                barr1[total] = arr1[i];
                barr2[total] = arr2[i];
                total ++;
            }
        }
    }
    if (pA > pB) {
        int xx = (pA - pB) * size;
        int x1, i2, i6, tt = 0;
        for (int i = 0, i2 = 0, i6 = 0; i != size && xx > 0; i ++) {
            if (arr1[i] == 2 * m1 && arr2[i] == 1) { i2 ++; xx --;}
            else if (arr1[i] == 1 && arr2[i] == 2 * (1 - m2)) {i6 ++; xx --;}
            else {
                t1[tt] = arr1[i];
                t2[tt] = arr2[i];
                tt ++;
            }
        }
        x1 = (pB - 0.5) * size - 0.5 * i2;
        for (int i = 0; i != tt && x1 > 0; i ++) {
            if (arr1[i] == 2 * m1 && arr2[i] == 2 * m2) x1 --;
            else {
                barr1[total] = arr1[i];
                barr2[total] = arr2[i];
                total ++;
            }
        }
    }
    if (pA < pB) {
        int xx = (pB - pA) * size;
        int x1, i4, i8, tt = 0;
        for (int i = 0, i4 = 0, i8 = 0; i != size && xx > 0; i ++) {
            if (arr1[i] == 1 && arr2[i] == 2 * m2) { i4 ++; xx --;}
            else if (arr1[i] == 2 * (1 - m1) && arr2[i] == 1) {i8 ++; xx --;}
            else {
                t1[tt] = arr1[i];
                t2[tt] = arr2[i];
                tt ++;
            }
        }
        x1 = (pA - 0.5) * size - 0.5 * i4;
        for (int i = 0; i != tt && x1 > 0; i ++) {
            if (arr1[i] == 2 * m1 && arr2[i] == 2 * m2) x1 --;
            else {
                barr1[total] = arr1[i];
                barr2[total] = arr2[i];
                total ++;
            }
        }
    }
    delete[] t1;
    delete[] t2;

    return total;
}

float fcorr(int arr1[], int arr2[], int size) {
    int* farr1, *farr2; int total;
    farr1 = new int[size];
    farr2 = new int[size];
    total = filt(arr1, arr2, farr1, farr2, size);
    float mean1 = mean(farr1, total);
    float mean2 = mean(farr2, total);
    float numerator = 0, denominator1 = 0, denominator2 = 0;

    // cout << "original length: " << size << endl;
    // cout << "filtered length: " << total << endl;
    // cout << "filted array 1: " << endl;
    // for (int i = 0; i != total; i ++) {
    //     cout << farr1[i] << ",";
    // }
    // cout << endl;
    // cout << "filted array 2: " << endl;
    // for (int i = 0; i != total; i ++) {
    //     cout << farr2[i] << ",";
    // }
    // cout << endl;

    for (int i = 0; i != total; i ++) {
        numerator += ((float)farr1[i] - mean1) * ((float)farr2[i] - mean2);
        denominator1 += pow(((float)farr1[i] - mean1), 2);
        denominator2 += pow(((float)farr2[i] - mean2), 2);
    }
    if (denominator1 == 0 || denominator2 == 0)
        return 0;
    float r = numerator / sqrt(denominator1 * denominator2);

    delete[] farr1;
    delete[] farr2;

    return r;
}

float bcorr(int arr1[], int arr2[], int size) {
    int* barr1, * barr2; int total;
    barr1 = new int[size];
    barr2 = new int[size];
    total = balance(arr1, arr2, barr1, barr2, size);
    float mean1 = mean(barr1, total);
    float mean2 = mean(barr2, total);
    float numerator = 0, denominator1 = 0, denominator2 = 0;

    for (int i = 0; i != total; i ++) {
        numerator += ((float)barr1[i] - mean1) * ((float)barr2[i] - mean2);
        denominator1 += pow(((float)barr1[i] - mean1), 2);
        denominator2 += pow(((float)barr2[i] - mean2), 2);
    }
    if (denominator1 == 0 || denominator2 == 0)
        return 0;
    float r = numerator / sqrt(denominator1 * denominator2);

    return r;
}

float corr(int arr1[], int arr2[], int size) {
    float mean1 = mean(arr1, size);
    float mean2 = mean(arr2, size);
    // float mean1 = 0, mean2 = 0;
    float numerator = 0, denominator1 = 0, denominator2 = 0;

    for (int i = 0; i != size; i ++) {
    //     if (arr1[i] != -1 && arr2[i] != -1) {
            numerator += ((float)arr1[i] - mean1) * ((float)arr2[i] - mean2);
            denominator1 += pow(((float)arr1[i] - mean1), 2);
            denominator2 += pow(((float)arr2[i] - mean2), 2);
    //     }
    }
    if (denominator1 == 0 || denominator2 == 0)
        return 0;
    float r = numerator / sqrt(denominator1 * denominator2);
    // if (r > 1 || r < -1) 
    //     cout << "r: " << r << ", fix me" << endl;
    
    // normalize
    
    // float pA = 0, pB = 0;
    // for (int i = 0; i != size; i ++) {
    //     if (arr1[i] == 0) pA += 2;
    //     else if (arr1[i] == 1) pA += 1;
    //     if (arr2[i] == 0) pB += 2;
    //     else if (arr2[i] == 1) pB += 1;
    // }
    // pA /= (2 * size);
    // pB /= (2 * size);

    // pA = pA >= 0.5 ? pA : (1 - pA);
    // pB = pB >= 0.5 ? pB : (1 - pB);
    
    // r /= (pA * pB);
    // 
    // // r /= temp;
    // if (r > 1) {
    //     for (int i = 0; i != size; i ++) 
    //         cout << arr1[i] << " ";
    //     cout << endl;
    //     for (int i = 0; i != size; i ++)
    //         cout << arr2[i] << " ";
    //     cout << endl;
    //     cout << "numerator: " << numerator << endl;
    //     cout << "denominator1: " << denominator1 << endl;
    //     cout << "denominator2: " << denominator2 << endl;
    //     cout << "pA: " << pA << endl;
    //     cout << "pB: " << pB << endl;
    //     cout << r << endl;
    //     exit(1);
    // }
    return r;
}

// namespace operations_research {
// float ibcorr(int arr1[], int arr2[], int size, float sigma) {
//     int obs[9] = {0}; // 0: AABB, 1: AABb, 3: AAbb, 4: AaBB, 5: AaBb, 6: Aabb, 7: aaBB, 8: aaBb, 9: aabb
//     for (int i = 0; i != size; i ++) {
//         if (arr1[i] == 0) { // AA
//             if (arr2[i] == 0) obs[0] += 1;
//             if (arr2[i] == 1) obs[1] += 1;
//             if (arr2[i] == 2) obs[2] += 1;
//         }       
//         if (arr1[i] == 1) { // Aa 
//             if (arr2[i] == 0) obs[3] += 1;
//             if (arr2[i] == 1) obs[4] += 1;
//             if (arr2[i] == 2) obs[5] += 1;
//         }
//         if (arr1[i] == 2) {
//             if (arr2[i] == 0) obs[6] += 1;
//             if (arr2[i] == 1) obs[7] += 1;
//             if (arr2[i] == 2) obs[8] += 1;
//         }
//     }
//     // create the linear solver with the GLOP backend.
//     MPSolver solver("corr_withdraw", MPSolver::GLOP_LINEAR_PROGRAMMING);
//     // Create the variables d0 to d8
//     MPVariable* const d0 = solver.MakeNumVar(0.0, obs[0], "d0");
//     MPVariable* const d1 = solver.MakeNumVar(0.0, obs[1], "d1");
//     MPVariable* const d2 = solver.MakeNumVar(0.0, obs[2], "d2");
//     MPVariable* const d3 = solver.MakeNumVar(0.0, obs[3], "d3");
//     MPVariable* const d4 = solver.MakeNumVar(0.0, obs[4], "d4");
//     MPVariable* const d5 = solver.MakeNumVar(0.0, obs[5], "d5");
//     MPVariable* const d6 = solver.MakeNumVar(0.0, obs[6], "d6");
//     MPVariable* const d7 = solver.MakeNumVar(0.0, obs[7], "d7");
//     MPVariable* const d8 = solver.MakeNumVar(0.0, obs[8], "d8");
//     // Create linear constraint: -sigma < 2fA - 1 < sigma; -sigma < 2fB - 1 < sigma
//     //                           2[] - (1 + sigma)total < 0 & 2[] - (1-sigma)total > 0
//     // for fA
//     int neg_obs_fA = 2 * obs[0] + 2 * obs[1] + 2 * obs[2] + obs[3] + obs[4] + obs[5];
//     MPConstraint* const ct1 = solver.MakeRowConstraint(-size, (1 + 2 * sigma) * size - neg_obs_fA, "ct1");
//     ct1->SetCoefficient(d0, 2 * sigma - 1);
//     ct1->SetCoefficient(d1, 2 * sigma - 1);
//     ct1->SetCoefficient(d2, 2 * sigma - 1);
//     ct1->SetCoefficient(d3, 2 * sigma);
//     ct1->SetCoefficient(d4, 2 * sigma);
//     ct1->SetCoefficient(d5, 2 * sigma);
//     ct1->SetCoefficient(d6, 1 + 2 * sigma);
//     ct1->SetCoefficient(d7, 1 + 2 * sigma);
//     ct1->SetCoefficient(d8, 1 + 2 * sigma);
//     MPConstraint* const ct2 = solver.MakeRowConstraint((1 - 2 * sigma) * size - neg_obs_fA, size, "ct2");
//     ct2->SetCoefficient(d0, -2 * sigma - 1);
//     ct2->SetCoefficient(d1, -2 * sigma - 1);
//     ct2->SetCoefficient(d2, -2 * sigma - 1);
//     ct2->SetCoefficient(d3, -2 * sigma);
//     ct2->SetCoefficient(d4, -2 * sigma);
//     ct2->SetCoefficient(d5, -2 * sigma);
//     ct2->SetCoefficient(d6, 1 - 2 * sigma);
//     ct2->SetCoefficient(d7, 1 - 2 * sigma);
//     ct2->SetCoefficient(d8, 1 - 2 * sigma);
//     // for fB
//     int neg_obs_fB = 2 * obs[0] + 2 * obs[3] + 2 * obs[6] + obs[1] + obs[4] + obs[7];
//     MPConstraint* const ct3 = solver.MakeRowConstraint(-size, (1 + 2 * sigma) * size - neg_obs_fB, "ct3");
//     ct3->SetCoefficient(d0, 2 * sigma - 1);
//     ct3->SetCoefficient(d1, 2 * sigma);
//     ct3->SetCoefficient(d2, 2 * sigma + 1);
//     ct3->SetCoefficient(d3, 2 * sigma - 1);
//     ct3->SetCoefficient(d4, 2 * sigma);
//     ct3->SetCoefficient(d5, 2 * sigma + 1);
//     ct3->SetCoefficient(d6, 2 * sigma - 1);
//     ct3->SetCoefficient(d7, 2 * sigma);
//     ct3->SetCoefficient(d8, 2 * sigma + 1);
//     MPConstraint* const ct4 = solver.MakeRowConstraint((1 - 2 * sigma) * size - neg_obs_fB, size, "ct4");
//     ct2->SetCoefficient(d0, -2 * sigma - 1);
//     ct2->SetCoefficient(d1, -2 * sigma);
//     ct2->SetCoefficient(d2, -2 * sigma + 1);
//     ct2->SetCoefficient(d3, -2 * sigma - 1);
//     ct2->SetCoefficient(d4, -2 * sigma);
//     ct2->SetCoefficient(d5, -2 * sigma + 1);
//     ct2->SetCoefficient(d6, -2 * sigma - 1);
//     ct2->SetCoefficient(d7, -2 * sigma);
//     ct2->SetCoefficient(d8, -2 * sigma + 1);
//     // Create objective function
//     MPObjective* const objective = solver.MutableObjective();
//     objective->SetCoefficient(d0, 1);
//     objective->SetCoefficient(d1, 1);
//     objective->SetCoefficient(d2, 1);
//     objective->SetCoefficient(d3, 1);
//     objective->SetCoefficient(d4, 1);
//     objective->SetCoefficient(d5, 1);
//     objective->SetCoefficient(d6, 1);
//     objective->SetCoefficient(d7, 1);
//     objective->SetCoefficient(d8, 1);
//     objective->SetMinimization();
// 
//     solver.Solve();
//     
//     // with draw values
//     int withdraw[9] = { \
//         d0->solution_value(), \
//         d1->solution_value(), \
//         d2->solution_value(), \
//         d3->solution_value(), \
//         d4->solution_value(), \
//         d5->solution_value(), \
//         d6->solution_value(), \
//         d7->solution_value(), \
//         d8->solution_value() \
//     };
//     int total_withdraw = 0;
//     for (int i = 0; i != 0; i ++) total_withdraw += withdraw[i];
//     
//     // drop values
//     int *warr1 = new int[size - total_withdraw];
//     int *warr2 = new int[size - total_withdraw];
//     int wcurr[9] = {0};
//     int curr = 0;
//     for (int i = 0; i != size; i ++) {
//         if (arr1[i] == 0) { // AA
//             if (arr2[i] == 0) {
//                 if (wcurr[0] < withdraw[0]) { wcurr[0] ++; continue; }
//             }
//             if (arr2[i] == 1) {
//                 if (wcurr[1] < withdraw[1]) { wcurr[1] ++; continue; }
//             }
//             if (arr2[i] == 2) {
//                 if (wcurr[2] < withdraw[2]) { wcurr[2] ++; continue; }
//             }
//         }       
//         if (arr1[i] == 1) { // Aa 
//             if (arr2[i] == 0) {  
//                 if (wcurr[3] < withdraw[3]) { wcurr[3] ++; continue; }
//             } 
//             if (arr2[i] == 1) { 
//                 if (wcurr[4] < withdraw[4]) { wcurr[4] ++; continue; }
//             } 
//             if (arr2[i] == 2) {
//                 if (wcurr[5] < withdraw[5]) { wcurr[5] ++; continue; }
//             } 
//         }
//         if (arr1[i] == 2) {
//             if (arr2[i] == 0) {
//                 if (wcurr[6] < withdraw[6]) { wcurr[6] ++; continue; }
//             }
//             if (arr2[i] == 1) {
//                 if (wcurr[7] < withdraw[7]) { wcurr[7] ++; continue; }
//             }
//             if (arr2[i] == 2) {
//                 if (wcurr[8] < withdraw[8]) { wcurr[8] ++; continue; }
//             }
//         }
//         warr1[curr] = arr1[i]; 
//         warr2[curr] = arr2[i]; 
//         curr++;
//     }
//     // ---- debug
//     cout << "obs0: " << obs[0] << ", d0: " << withdraw[0] << ", AABB" << endl;
//     cout << "obs1: " << obs[1] << ", d1: " << withdraw[1] << ", AABb" << endl;
//     cout << "obs2: " << obs[2] << ", d2: " << withdraw[2] << ", AAbb" << endl;
//     cout << "obs3: " << obs[3] << ", d3: " << withdraw[3] << ", AaBB" << endl;
//     cout << "obs4: " << obs[4] << ", d4: " << withdraw[4] << ", AaBb" << endl;
//     cout << "obs5: " << obs[5] << ", d5: " << withdraw[5] << ", Aabb" << endl;
//     cout << "obs6: " << obs[6] << ", d6: " << withdraw[6] << ", aaBB" << endl;
//     cout << "obs7: " << obs[7] << ", d7: " << withdraw[7] << ", aaBb" << endl;
//     cout << "obs8: " << obs[8] << ", d8: " << withdraw[8] << ", aabb" << endl;
//     cout << "total withdraw number: " << total_withdraw << endl;
//     // ----
// 
//     float ret = corr(warr1, warr2, size - total_withdraw);
//     delete[] warr1; 
//     delete[] warr2;
//     exit(1);
//     return ret;
// }
// } // namespace operations_research
