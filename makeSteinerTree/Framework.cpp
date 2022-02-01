//
// Created by Spencer Halberg on 9/23/21.
//


#include "Framework.h"
#define INF 99999;

Framework::Framework()
:f_numNodes(0)
{}

Framework::~Framework() {
    set<list<string>*> deleted_ptrs;
    auto iter = f_paths.begin();
    while (iter != f_paths.end()){
        list<string>* list_ptr=iter->second;
        if(deleted_ptrs.find(list_ptr) == deleted_ptrs.end()){
            deleted_ptrs.insert(list_ptr);
            delete list_ptr;
        }
        iter++;
    }
}


int Framework::init(string &g_fname, string &nodeList) {
    list<string> edges;
    set<string> nodes;
    string line;
    ifstream graphFile;
    graphFile.open(g_fname);
    if(!graphFile.is_open()){
        cerr<<"Error: Unable to open graph file."<<endl;
        return -99;
    }
    while(getline(graphFile,line)){
        edges.push_back(line);
        istringstream iss(line);
        while(iss >> line){
            nodes.insert(line);
        }
    }
    Graph g;
    g.init(edges, nodes);
    graphFile.close();
    cout << "Network contains " << g.getNumNodes() << " nodes." <<endl;

    f_numNodes=nodes.size();
    findLargestConnectedComp(g);
    f_numNodes=f_network.getNodeSet()->size();
    cout << "Largest connected component contains " << f_network.getNumNodes() <<endl;
    ifstream nodeFile;
    nodeFile.open(nodeList);
    if(!nodeFile.is_open()){
        cerr << "Error: Unable to open node list file." << endl;
        return -99;
    }

    int num_genes_interest=0;
    int found_genes=0;
    while(getline(nodeFile, line)){
        num_genes_interest++;
        if(f_network.checkNode(line)) {
            found_genes++;
            f_geneList.insert(line);
            f_startNodes.push(line);
        }
    }
    cout << "Searching for " << num_genes_interest <<" nodes in largest connected component. Found " << found_genes << " genes." <<endl;
    return 0;
}

int Framework::construct_tree() {
    set<string> nodes2add=f_geneList;
    for(const string& startNode:f_geneList) {
        find_shortest_paths(startNode, nodes2add);
    }


    string start_node;
    string end_node;
    int dist = INF;
    auto start_iter=nodes2add.begin();
    while(start_iter != nodes2add.end()--) {
        auto target_iter=start_iter;
        target_iter++;
        while (target_iter != nodes2add.end()){
            if (f_distance[pair<string, string>(*start_iter, *target_iter)] < dist) {
                start_node=*start_iter;
                end_node=*target_iter;
                dist = f_distance[pair<string, string>(*start_iter, *target_iter)];
            }
            target_iter++;
        }
        start_iter++;
    }

    add_path(start_node, end_node);
    nodes2add.erase(start_node);
    nodes2add.erase(end_node);


    while (!nodes2add.empty()) {
        dist = INF;
        auto start_iter = nodes2add.begin();
        auto target_iter = f_steinerTree.getNodeSet()->begin();
        while(start_iter != nodes2add.end()) {
            while (target_iter != f_steinerTree.getNodeSet()->end()){
                if(*start_iter != *target_iter) {
                    if (f_distance[pair<string, string>(*start_iter, *target_iter)] < dist) {
                        start_node = *start_iter;
                        end_node = *target_iter;
                        dist = f_distance[pair<string, string>(*start_iter, *target_iter)];
                    }
                }
                target_iter++;
            }
            start_iter++;
        }
        add_path(start_node, end_node);
        nodes2add.erase(start_node);
    }
    return 0;
}

int Framework::add_path(const string &string1, const string &string2){
    list<string>* path=f_paths[pair<string, string> (string1, string2)];
    auto fIter=path->cbegin();
    auto sIter=path->cbegin();
    sIter++;
    string first_node=*fIter;
    string second_node=*sIter;

    if(!f_steinerTree.checkNode(first_node)){
        f_steinerTree.addNode(first_node);
        find_shortest_paths(first_node, f_geneList);
    }

    while(sIter != path->cend()){
        first_node=*fIter;
        second_node=*sIter;
        if(!f_steinerTree.checkNode(second_node)){
            f_steinerTree.addNode(second_node);
            find_shortest_paths(second_node, f_geneList);
        }
        f_steinerTree.addEdge(first_node, second_node);
        fIter++;
        sIter++;
    }
    return 0;
}


int Framework::find_shortest_paths(const string& startNode, set<string> &targets) {
    vector<char> color ( f_numNodes, 'W');
    vector<int> distance(f_numNodes, numeric_limits<int>::max());
    vector<int> predecessor(f_numNodes, numeric_limits<int>::min());
    int start_node_index=f_network.getIndex(startNode);
    queue<int> traverseList;
    int pred_node=0;
    set<string> foundNodes;


    color[start_node_index]='G';
    distance[start_node_index]=0;
    traverseList.push(start_node_index);

    while(!traverseList.empty() && foundNodes != targets){
        pred_node=traverseList.front();
        traverseList.pop();
        for(int neighbor:*f_network.getNeighbors(pred_node)){
            if(color[neighbor]=='W'){
                color[neighbor]='G';
                distance[neighbor]=distance[pred_node]+1;
                predecessor[neighbor]=pred_node;
                traverseList.push(neighbor);
                if(f_geneList.find(f_network.getLabel(neighbor)) != f_geneList.end()){
                    foundNodes.insert(f_network.getLabel(neighbor));
                }
            }
        }
        color[pred_node]='B';
    }

    for(const string& node: foundNodes) {
        if(node != startNode) {
            pair<string, string> key1;
            pair<string, string> key2;
            key1 = make_pair(startNode, node);
            key2 = make_pair(node, startNode);
            if (color[f_network.getIndex(node)] == 'W') {
                int dist=INF
                if(f_distance.find(key1) == f_distance.end()) {
                    f_distance.insert(pair<pair<string, string>, int>(key1, dist));
                }
                if(f_distance.find(key2) == f_distance.end()) {
                    f_distance.insert(pair<pair<string, string>, int>(key2, dist));
                }
            } else {
                if(f_distance.find(key1) == f_distance.end()){
                    f_distance.insert(pair<pair<string,string>, int>(key1, distance[f_network.getIndex(node)]));
                }
                if(f_distance.find(key2) == f_distance.end()) {
                    f_distance.insert(pair<pair<string, string>, int>(key2, distance[f_network.getIndex(node)]));
                }
                list<string> *path;
                path = new list<string>;
                path->push_front(node);
                int curr_node = f_network.getIndex(node);
                while (curr_node != f_network.getIndex(startNode)) {
                    curr_node = predecessor[curr_node];
                    path->push_front(f_network.getLabel(curr_node));
                }
                if(f_paths.find(key1) != f_paths.end() && f_paths.find(key2) != f_paths.end()){
                    delete path;
                }
                if(f_paths.find(key1) == f_paths.end()) {
                    f_paths.insert(pair<pair<string, string>, list<string> *>(key1, path));
                }
                if(f_paths.find(key2) == f_paths.end()) {
                    f_paths.insert(pair<pair<string, string>, list<string> *>(key2, path));
                }
            }
        }
    }


    return 0;
}

int Framework::find_shortest_paths(const string  &startNode, set<string> &targets, set<string> *nodes2add) {
    vector<char> color ( f_numNodes, 'W');
    vector<int> distance(f_numNodes, numeric_limits<int>::max());
    vector<int> predecessor(f_numNodes, numeric_limits<int>::min());
    int start_node_index=f_network.getIndex(startNode);
    queue<int> traverseList;
    int pred_node=0;
    set<string> foundNodes;


    color[start_node_index]='G';
    distance[start_node_index]=0;
    traverseList.push(start_node_index);

    while(!traverseList.empty() && foundNodes != targets){
        pred_node=traverseList.front();
        traverseList.pop();
        for(int neighbor:*f_network.getNeighbors(pred_node)){
            if(color[neighbor]=='W'){
                color[neighbor]='G';
                distance[neighbor]=distance[pred_node]+1;
                predecessor[neighbor]=pred_node;
                traverseList.push(neighbor);
                if(f_geneList.find(f_network.getLabel(neighbor)) != f_geneList.end()){
                    foundNodes.insert(f_network.getLabel(neighbor));
                }
            }
        }
        color[pred_node]='B';
    }
    if (foundNodes.empty()){
        nodes2add->erase(startNode);
        return 0;
    }

    for(const string& node: targets) {
        if(node != startNode) {
            pair<string, string> key1;
            pair<string, string> key2;
            key1 = make_pair(startNode, node);
            key2 = make_pair(node, startNode);
            if (color[f_network.getIndex(node)] == 'W') {
                int dist = INF
                if(f_distance.find(key1) == f_distance.end()) {
                    f_distance.insert(pair<pair<string, string>, int>(key1, dist));
                }
                if(f_distance.find(key2) == f_distance.end()) {
                    f_distance.insert(pair<pair<string, string>, int>(key2, dist));
                }
            } else {
                if(f_distance.find(key1) == f_distance.end()){
                    f_distance.insert(pair<pair<string,string>, int>(key1, distance[f_network.getIndex(node)]));
                }
                if(f_distance.find(key2) == f_distance.end()) {
                    f_distance.insert(pair<pair<string, string>, int>(key2, distance[f_network.getIndex(node)]));
                }
                list<string> *path;
                path = new list<string>;
                path->push_front(node);
                int curr_node = f_network.getIndex(node);
                while (curr_node != f_network.getIndex(startNode)) {
                    curr_node = predecessor[curr_node];
                    path->push_front(f_network.getLabel(curr_node));
                }
                if(f_paths.find(key1) != f_paths.end() && f_paths.find(key2) != f_paths.end()){
                    delete path;
                }
                if(f_paths.find(key1) == f_paths.end()) {
                    f_paths.insert(pair<pair<string, string>, list<string> *>(key1, path));
                }
                if(f_paths.find(key2) == f_paths.end()) {
                    f_paths.insert(pair<pair<string, string>, list<string> *>(key2, path));
                }
            }
        }
    }
    return 0;
}

int Framework::printSteinerTree(string& outfile, Graph* trueGraph) {
    f_steinerTree.printGraph(outfile, trueGraph);
    return 0;
}

int Framework::printOrginalGraph(string& outfile, Graph* trueGraph){
    f_network.printGraph(outfile, trueGraph);
    return 0;
}

Graph* Framework::getGraph() {
    return &f_network;
}

int Framework::findLargestConnectedComp(Graph &g) {
    bool found_max_comp=false;
    Graph connect_comp;
    set<string>::iterator start_node=g.getNodeSet()->begin();
    while(!found_max_comp) {
        vector<char> color(f_numNodes, 'W');
        vector<int> distance(f_numNodes, numeric_limits<int>::max());
        vector<int> predecessor(f_numNodes, numeric_limits<int>::min());
        int start_node_index = g.getIndex(*start_node);
        queue<int> traverseList;
        int pred_node = 0;
        set<string> foundNodes;


        color[start_node_index] = 'G';
        distance[start_node_index] = 0;
        traverseList.push(start_node_index);
        int found_genes = 0;
        while (!traverseList.empty()) {
            pred_node = traverseList.front();
            traverseList.pop();
            for (int neighbor:*g.getNeighbors(pred_node)) {
                if (color[neighbor] == 'W') {
                    color[neighbor] = 'G';
                    distance[neighbor] = distance[pred_node] + 1;
                    predecessor[neighbor] = pred_node;
                    traverseList.push(neighbor);
                    found_genes++;
                }
            }
            color[pred_node] = 'B';
        }
        if(found_genes>=g.getNodeSet()->size()/2){
            found_max_comp = true;
            string gene_name;
            string neighbor;
            for(int i=0; i<color.size(); i++){
                if(color[i] != 'W'){
                    gene_name=g.getLabel(i);
                    f_network.addNode(gene_name);
                }
            }
        }else{
            start_node++;
        }
    }
    f_network.addEdges(g.getEdgeDirections());

    return 0;
}


int main (int argc, char* argv[])
{
    string g_fname=argv[1];
    string nodeList=argv[2];
    string outFile=argv[3];
    string outFile2 = "aspergillus_network_I02";

    Framework fw;
    fw.init(g_fname, nodeList);
    cout << "Constructing Steiner Tree" << endl;
    fw.construct_tree();

    cout << "Printing Steiner Tree" << endl;
    fw.printSteinerTree(outFile, fw.getGraph());
    //fw.printOrginalGraph(outFile2, fw.getGraph());
    return 0;
}