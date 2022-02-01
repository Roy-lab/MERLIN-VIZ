//
// Created by Spencer Halberg on 10/10/20.
//


#include "Graph.h"

Graph::Graph(){
    m_numNodes=0;
};

Graph::~Graph(){
    for(set<int>* neighbors:m_adjacencyList){
        delete neighbors;
    }
}

int Graph::init(const list<string>& edgeList, const set<string>& nodeList) {
    m_numNodes=0;
    for(string node:nodeList){
       addNode(node);
    }
    string Label1;
    string Label2;
    for(string edge:edgeList){
        istringstream iss(edge);
        iss >> edge;
        Label1=edge;
        iss >> edge;
        Label2=edge;
        addEdge(Label1, Label2);
    }
    return 0;
}

int Graph::addNode(string &node) {
    m_Label2Index[node]= m_numNodes;
    m_Index2Label[m_numNodes]=node;
    m_Nodes.insert(node);
    m_adjacencyList.push_back(new set<int>);
    m_numNodes+=1;
    return 0;
}

int Graph::addEdge(string &node1, string &node2){
    int index1=getIndex(node1);
    int index2=getIndex(node2);
    if (m_adjacencyList[index1]->find(index2) == m_adjacencyList[index1]->end()) {
        m_adjacencyList[index1]->insert(index2);
    }
    m_edgeDirection.insert(pair<string, string>(node1, node2) );
    if (m_adjacencyList[index2]->find(index1) == m_adjacencyList[index2]->end()){
        m_adjacencyList[index2]->insert(index1);
    }
    return 0;
}

int Graph::addEdges(set<pair<string, string>> directions) {
    for(pair<string, string> d:directions){
        string node1=d.first;
        string node2=d.second;
    if (m_Nodes.find(node1)!=m_Nodes.end() && m_Nodes.find(node2)!=m_Nodes.end()) {
            m_edgeDirection.insert(d);
            int index1 = getIndex(node1);
            int index2 = getIndex(node2);
            if (m_adjacencyList[index1]->find(index2) == m_adjacencyList[index1]->end()) {
                m_adjacencyList[index1]->insert(index2);
            }
            if (m_adjacencyList[index2]->find(index1) == m_adjacencyList[index2]->end()) {
                m_adjacencyList[index2]->insert(index1);
            }
        }
    }
    return 0;
}


int Graph::getIndex(const string &Label){return m_Label2Index[Label];}

string Graph::getLabel(const int Index){return m_Index2Label[Index];}

set<int>* Graph::getNeighbors(const int Index){return m_adjacencyList[Index];}

bool Graph::checkNode(const string &node) {return m_Label2Index.count(node)>0;}

bool Graph::checkNode(int node) {return m_Index2Label.count(node)>0;}

bool Graph::noNeighbors(const string &node) {return m_adjacencyList[m_Label2Index[node]]->empty();}

set<string> *Graph::getNodeSet() {
    return &m_Nodes;
}

int Graph::printGraph(const string& outprefix, Graph* true_graph) {
    stringstream adj_stream;
    for(int i=0; i<m_adjacencyList.size(); i++){
        adj_stream << m_Index2Label[i] << "\t";
        for(int node: *m_adjacencyList[i]){
            adj_stream << m_Index2Label[node] << "\t";
        }
        adj_stream << "\n";
    }

    ofstream fout;
    fout.open(outprefix+"_adjacency_list.txt");
    if(!fout.is_open()){
        cerr<<"Error: Unable to open graph file."<<endl;
        return -99;
    }

    fout << adj_stream.str();
    fout.close();

    stringstream node_stream;
    for(const string& node:*getNodeSet()){
        node_stream << node << "\n";
    }

    fout.open(outprefix+"_node_list.txt");
    if(!fout.is_open()){
        cerr<<"Error: Unable to open graph file."<<endl;
        return -99;
    }
    fout << node_stream.str();
    fout.close();
    
    /*
    stringstream node_stream2;
    int i=1;
    for(const string& node:*getNodeSet()){
        node_stream2 << i << " " << node << "\n";
        i++;
    }
    
     
    fout.open(outprefix+"_node_index.txt");
    if(!fout.is_open()){
        cerr<<"Error: Unable to open graph file."<<endl;
        return -99;
    }
    fout << node_stream2.str();
    fout.close();
    */

    stringstream edge_stream;
    set<pair<string, string>> edges;
    for(int i=0; i<m_adjacencyList.size(); i++){
        for(int node2: *m_adjacencyList[i]){
            pair<string, string> edge= pair<string,string>(m_Index2Label[i], m_Index2Label[node2]);
            if(edges.find(edge)==edges.end()) {
                edges.insert(edge);
                edges.insert(pair<string,string>(m_Index2Label[node2], m_Index2Label[i]));
                if(true_graph->getEdgeDirection(m_Index2Label[i], m_Index2Label[node2])) {
                    edge_stream << m_Index2Label[i] << "\t" << m_Index2Label[node2] << "\n";
                }
                if(true_graph->getEdgeDirection(m_Index2Label[node2], m_Index2Label[i])){
                    edge_stream << m_Index2Label[node2] << "\t" << m_Index2Label[i] << "\n";
                }
            }
        }
    }

    fout.open(outprefix+"_edge_list_names.txt");
    if(!fout.is_open()){
        cerr<<"Error: Unable to open graph file."<<endl;
        return -99;
    }
    fout << edge_stream.str();
    fout.close();

    /*
    stringstream edge_stream2;
    set<pair<int, int>> edges2;
    for(int i=0; i<m_adjacencyList.size(); i++){
        for(int node2: *m_adjacencyList[i]){
            pair<int, int> edge= pair<int, int>(i+1, node2+1);
            if(edges2.find(edge)==edges2.end()) {
                edges2.insert(edge);
                edges2.insert(pair<int, int>(node2+1, i+1));
                if(true_graph->getEdgeDirection(m_Index2Label[i], m_Index2Label[node2])) {
                    edge_stream2 << i+1 << " " << node2+1 << " 1" << "\n";
                }else{
                    edge_stream2 << node2+1 << " " << i+1 << " 1" << "\n";
                }
            }
        }
    }
    
    fout.open(outprefix+"_edge_list.txt");
    if(!fout.is_open()){
        cerr<<"Error: Unable to open graph file."<<endl;
        return -99;
    }
    fout << edge_stream2.str();
    fout.close();
    return 0;
    */
}

bool Graph::getEdgeDirection(string& node1, string& node2) {
    if(m_edgeDirection.find(pair<string, string>(node1, node2))!=m_edgeDirection.end()){
        return true;
    }else{
        return false;
    }
}

set<pair<string, string>> Graph::getEdgeDirections() {
    return m_edgeDirection;
}

int Graph::getNumNodes() {
    return m_numNodes;
}




