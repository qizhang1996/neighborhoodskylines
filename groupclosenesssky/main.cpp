#include <networkit/graph/Graph.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <iostream>
#include<fstream>
#include "SkylineGraph.hpp"
#include "imprGroupCloseness.hpp"

#define mp(a,b) make_pair(a,b)
using namespace std;

vector<NetworKit::node> skyvec;
int k;
bool usesky;
string filename;

void OutputRes(){
    int cnt=0;
    printf("-------------results-------------\n");
    for(int i=0; i<skyvec.size(); i++){
        if(skyvec[i]==1){
            cnt++;
            cout<<i<<"\n";
        }
    }
    cout<<"skyres.size()="<<cnt<<endl;
}

void readskyvec(const char* skyfile, NetworKit::node N){
    ifstream ifs;
    ifs.open(skyfile, ifstream::in);
    if(!ifs.is_open()){
        cout<<"error opening skyline result file\n";
        return;
    }
    NetworKit::node a;
    skyvec.resize(N, 0);
    while(ifs>>a){
        skyvec[a]=1;
    }
    //OutputRes();
}

void runalg(const NetworKit::Graph &g){
    INFO("-----------------------------");
    INFO("file = ", filename, ", k = ", k, ", usesky = ", usesky);
    imprGroupCloseness gc(g, k);
    gc.run(usesky, skyvec);
    gc.outputTime();
    auto apx = gc.groupMaxCloseness();
    sort(apx.begin(), apx.end());
    for (NetworKit::count i = 0; i < k; ++i) {
        INFO(apx[i]);
    }
    INFO("-----------------------------");
}

int main(int argc, char** argv){

	//read graph
    ifstream ifs;
    ifs.open(argv[1], ifstream::in);
    if(!ifs.is_open()){
        cout<<"error opening file\n";
    }
    NetworKit::node a, b;
    vector<pair<NetworKit::node, NetworKit::node>>belongs;
    NetworKit::node N=0;
    while(ifs>>a>>b){
        belongs.push_back(mp(a,b));
        if(a>N) N=a;
        if(b>N) N=b;
    }
    N++;
    NetworKit::Graph g(N, false, false);
    for(auto pr:belongs){
        g.addEdge(pr.first, pr.second);
    }
    INFO("N=", N);
    filename=argv[1];
    readskyvec(argv[2], N);
    k=atoi(argv[3]);
    usesky=atoi(argv[4]);
    cout<<"k="<<k<<endl;
    cout<<"usesky="<<usesky<<endl;
    omp_set_num_threads(1);
    runalg(g);
 	return 0;

}