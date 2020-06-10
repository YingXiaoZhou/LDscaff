#include <iostream>
#include <lemon/list_graph.h>

using namespace std;
using namespace lemon;

int main() {
    ListDigraph g;

    ListDigraph::Node u = g.addNode();
    ListDigraph::Node v = g.addNode();
    ListDigraph::Arc  a = g.addArc(u, v);

    cout << "directed graph with " << countNodes(g) << " nodes and " << countArcs(g) << " arc." << endl;

    return 0;
}

