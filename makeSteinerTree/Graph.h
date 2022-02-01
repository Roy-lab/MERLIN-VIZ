//
// Created by Spencer Halberg on 10/10/20.
//

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include <vector>
#include <map>


#ifndef PROBLEM4_GRAPH_H
#define PROBLEM4_GRAPH_H



using namespace std;

class Graph {
public:
    Graph();
    ~Graph();
    int init(const list<string>&, const set<string>&);
    int addNode(string&);
    int addEdge(string&, string&);
    int addEdges(set<pair<string, string>> directions);
    int getIndex(const string &Label);
    string getLabel(int Index);
    set<int>* getNeighbors(int Index);
    int getNumNodes();
    set<string>* getNodeSet();
    bool getEdgeDirection(string&, string&);
    set<pair<string, string>> getEdgeDirections();
    bool checkNode(const string &node);
    bool checkNode(int node);
    bool noNeighbors(const string &node);
    int printGraph(const string&, Graph*);
private:
    int m_numNodes;
    set<string> m_Nodes;
    set<pair<string, string>> m_edgeDirection;
    map<string, int> m_Label2Index;
    map<int, string> m_Index2Label;
    vector<set<int>*> m_adjacencyList;
};


#endif //PROBLEM4_GRAPH_H
