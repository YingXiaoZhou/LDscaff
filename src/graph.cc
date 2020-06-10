#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <assert.h>
#include "tools.h"
#include "cal.h"

#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/math.h>

#define INF 65535
#define NSAMPLE 1000

using namespace std;
using namespace lemon;

typedef MaxWeightedMatching<SmartGraph, SmartGraph::EdgeMap<float> > MWM;

static int nmarker = 100;

struct scafGraph{
    string *vexs;
    float **edge;
    int nvex, nedge;
    
    scafGraph(int n) : nvex(n) {
        nedge = nvex * nvex - nvex;
        vexs = new string[nvex];
        edge = new float*[nvex];
        for (int i = 0; i != nvex; i ++)
            edge[i] = new float[nvex];
    }

    ~scafGraph() {
        if (vexs)
            delete[] vexs;
        for (int i = 0; i != nedge; i ++)
            delete[] edge[i];
        delete[] edge;
    }

    void printGraph() {
        for (int i = 0; i != nvex; i ++) {
            for (int j = 0; j != nvex; j ++) 
                cout << edge[i][j] << "\t";
        cout << endl;
        }
    }

    string to_lgf() {
        ostringstream ss;
        // node part
        ss << "@nodes\n";
        ss << "label\n";
        for (int i = 0; i != nvex; i ++)
            ss << i << "\n";

        // edge part
        ss << "@edges\n";
        ss << " \t \tlabel\tweight\n";
        int cc = 0;
        for (int i = 0; i != nvex; i ++) {
            for (int j = 0; j != nvex; j ++) {
                ss << i << "\t" << j << "\t" << cc << "\t" << edge[i][j] << "\n";
                cc += 1;
            }
        }
        return ss.str();
    }
};

scafGraph *createGraph(char** infList, int nsample, int nscaf, ofstream &logf) {
    scafGraph *sg = new scafGraph(2 * nscaf);
    int ***elems = new int**[2 * nscaf]; // each scaffold 2 ends
    for (int i = 0; i != 2 * nscaf; i ++) {
        elems[i] = new int*[nmarker]; // each end nmarker markers
        for (int j = 0; j != nmarker; j ++) {
            elems[i][j] = new int[nsample]; // sample number, array for corr cal
            // init elems with -9 (neg value while some scaffolds may have zero markers or less markers than nmarkers)
            for (int k = 0; k != nsample; k ++) {
                elems[i][j][k] = -9;
            }
        }
    }
    
    for (int i = 0; i != nsample; i ++) { // each file a sample
        // ifstream ifs(infList[i]);
        ifstream ifs(infList[i]);
        
        if (!ifs.is_open()) {
            cerr << "error opening file " << infList[i] << endl;
            exit(1);
        }
        cout << "processing file " << infList[i] << "... ";
        
        string data; 
        int ncount = 0;
        while (getline(ifs, data)) { // each line one scaffold
            if (data == "")
                continue;
            vector<string> tokens;
            parseLine(data, " ", tokens); // should be nmarker tokens
            // assert(tokens.size() >= nmarker);  // scaffold may have less markers than nmarkers
            // sample no.: i, scaf no.: ncount, marker no.: iterator
            for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it ++) {
                if ((int)(it-tokens.begin()) == nmarker)
                    break;
                // for phased, 
                int tmp = atoi((*it).c_str());
                if (tmp == 4) tmp = 3;
                elems[ncount][(int)(it-tokens.begin())][i] = tmp;
            }
            ncount ++;
        }

        cout << "done" << endl;
    }

    cout << "finish processing all files in file list" << endl;
    // cout << "calculation edges, loged in edge.log" << endl;
    // ofstream ledge("edge.log");
    for (int i = 0; i != 2 * nscaf; i ++) {
        for (int j = i; j != 2 * nscaf; j ++) {
            if (i == j) sg->edge[i][j] = 0;
            else if (i % 2 == 0 && j - i == 1) sg->edge[i][j] = 0;
            // if (i < j) directed? undirected?
            else {
                float scorr = 0;
                int comb_count = 0;
                for (int n1 = 0; n1 != nmarker; n1 ++) {
                    if (elems[i][n1][0] == -9) // scaffold i do not have n1 th marker
                        continue;
                    for (int n2 = 0; n2 != nmarker; n2 ++) {
                        if (elems[j][n2][0] == -9)
                            continue;
                        comb_count += 1;
                        // if (i == 0 && j == 40) {
                        //     cout << "scaffold " << i << ", marker " << n1 << endl;
                        //     for (int k = 0; k != nsample; k ++) cout << elems[i][n1][k] << " ";
                        //     cout << endl;
                        //     cout << "scaffold " << j << ", marker " << n2 << endl;
                        //     for (int k = 0; k != nsample; k ++) cout << elems[j][n2][k] << " ";
                        //     cout << endl;
                        // } 
                        float cor = fcorr(elems[i][n1], elems[j][n2], nsample);
                        scorr += abs(cor);
                    }
                }
                if (comb_count < 1 * nmarker) // consider too less loci, not efficient evidence
                    sg->edge[i][j] = 0;
                else {
                    sg->edge[i][j] = scorr / comb_count;
                    logf << i << "\t" << j << "\t" << sg->edge[i][j] << endl;
                }
            }
        }
    }

    for (int i = 0; i != 2 * nscaf; i ++) {
        for (int j = 0; j != nmarker; j ++) 
            delete[] elems[i][j];
        delete[] elems[i];
    }
    delete[] elems;

    cout << "finish construct graph" << endl;

    return sg;
}

// MWM *calMaxMatching(SmartGraph graph, SmartGraph::EdgeMap<float> weight) {
//     // SmartGraph graph;
//     // SmartGraph::EdgeMap<float> weight(graph);
//     // 
//     // istringstream lgfs(sg->to_lgf());
//     // graphReader(graph, lgfs).
//     //   edgeMap("weight", weight).run();
//     
//     // MaxWeightedMatching< SmartGraph, SmartGraph::EdgeMap<float> > *mwm = new MaxWeightedMatching< SmartGraph, SmartGraph::EdgeMap<float> >(graph, weight);
//     MWM *mwm = new MWM(graph, weight);
//     mwm->run();
//     
//     return mwm;
// }

int main(int argc, char* argv[]) {
    ifstream ifs(argv[1]);
    if (!ifs.is_open()) {
        cerr << "Cannot open category file " << argv[1] << endl;
        return 1;
    }

    string tmp;
    char *flist[NSAMPLE];
    for (int i = 0; i != NSAMPLE; i ++)
        flist[i] = new char[100];
    int scount = 0;
    while (ifs >> tmp) {
        // flist[scount] = tmp.c_str();
        strcpy(flist[scount], tmp.c_str());
        scount ++;
    }

    cout << "phasing file list done" << endl;

    int nsample = atoi(argv[2]);
    int nscaf = atoi(argv[3]);

    nmarker = atoi(argv[4]);
    ofstream logf(argv[5]);

    scafGraph *sg = createGraph(flist, nsample, nscaf, logf);
    SmartGraph graph;
    SmartGraph::EdgeMap<float> weight(graph);

    istringstream lgfs(sg->to_lgf());
    graphReader(graph, lgfs).
        edgeMap("weight", weight).run();
    
    MWM *mwm = new MWM(graph, weight);
    mwm->run();
    // MWM *mwm = calMaxMatching(graph, weight);

    for (SmartGraph::EdgeIt e(graph); e != INVALID; ++ e) {
        if (graph.id(graph.u(e)) >= graph.id(graph.v(e))) continue;
        if (mwm->matching(e))
            cout << graph.id(graph.u(e)) << "\t" << graph.id(graph.v(e)) << "\t" << weight[e] << endl;
    }
    
    return 0;
}
