//
// Created by Spencer Halberg on 9/23/21.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include "Graph.h"
#include <string>
#include <algorithm>
#include <queue>



#ifndef FIND_SHORTEST_PATH_FRAMEWORK_H
#define FIND_SHORTEST_PATH_FRAMEWORK_H
using namespace std;

class Framework {
public:
    Framework();

    ~Framework();

    int init(string &, string &);

    int construct_tree();

    int printSteinerTree(string &, Graph *);

    int printOrginalGraph(string &, Graph *);

    int findLargestConnectedComp(Graph &);


    Graph *getGraph();

private:
    int find_shortest_paths(const string &, set<string> &);
    int find_shortest_paths(const string &, set<string> &, set<string> *);
    Graph f_network;
    Graph f_steinerTree;
    int f_numNodes{};
    set<string> f_geneList;
    queue<string> f_startNodes;
    map<pair<string, string>, int> f_distance;
    map<pair<string, string>, list<string> *> f_paths;

    int add_path(const string &, const string &);
};

#endif //FIND_SHORTEST_PATH_FRAMEWORK_H
