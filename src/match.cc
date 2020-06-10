/* Author: zicheng*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <string.h>

#include "graph.h"

#define INF 65535

using namespace std;

const static char* myopts = "v:m:t:r:h";

static char* vcf_file   = NULL;
static char* tmp_folder = "/tmp";

static int nmarker = 100;

static float minor = 0.2;

static void Usage() {
    fprintf(stderr, 
            "Program: matching (Link scaffolds by LD)\n"
            "Version: 1.0\n"
            "Contact: Zicheng ZHAO <zachazhao2-c@my.cityu.edu.hk>\n"
            "\n"
            "Usage: mathcing -v <in.vcf> -m <nmarker> [options]\n"
            "\n"
            "Options: -t STR    temp directory [default: /tmp]\n"
            "         -r FLOAT  minor allele cutoff [default: 0.2]\n");
    return;
}

static void parse_opt(int argc, char** argv) {
    int ch;
    optarg = NULL;
    int  = 0;
    char* p;

    while ((ch = getopt(argc, argv, myopts)) != -1) {
        switch(ch) {
            case 'v':
                vcf_file = strdup(optarg);
                break;
            case 't':
                temp_folder = strdup(optarg);
                break;
            case 'm':
                nmarker = atoi(optarg);
                break;
            case 'r':
                minor = atod(optarg);
                break;
            case '?':
                fprintf(stderr, "Unrecognized option - %c\n", optopt);
                break;
            default:
                return 1;
        }
    }
    if (vcf_file == NULL) {
        fprintf(stderr, "Need vcf file");
        Usage();
        return 1;
    }
    
    return 0;
}

int main(int argc, char** argv) {
    parse_opt(argc, argv);

    clock_t start, finish;
    double totalTime;
    start = clock();

    // Vcf* vcf = new vcf(vcf_file, nimor, nmarker);
    // preprocessing, parse vcf
    

    finish = clock();
    totalTime = (double)(finish - start) / 1000;
    cout << "[main] Finish, total time: " << totalTime << endl;

    return 0;
}
