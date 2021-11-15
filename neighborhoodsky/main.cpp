#include <networkit/graph/Graph.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <iostream>
#include <cstring>
#include<fstream>
#include "SkylineGraph.hpp"
//#include "imprGroupCloseness.hpp"

#define mp(a,b) make_pair(a,b)
using namespace std;

vector<NetworKit::node> skyvec;
vector<NetworKit::node> skyvec2;
vector<NetworKit::node> skyvec3;
int k;
bool usesky;
int thread;
int step;
string filename;

void OutputRes(const char* str){
    ofstream ofs(str);
    int cnt=0;
    for(int i=0; i<skyvec.size(); i++){
        if(skyvec[i]==1){
            ofs<<i<<endl;
            cnt++;
        }
    }
    ofs.close();
    cout<<"skyvec.size()="<<cnt<<endl;
}

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


void runalg(const NetworKit::Graph &g){
    // INFO("-----------------------------");
    // INFO("file = ", filename, ", k = ", k, ", usesky = ", usesky);
    // imprGroupCloseness gc(g, k);
    // gc.run(usesky, skyvec, step);
    // gc.outputTime();
    // auto apx = gc.groupMaxCloseness();
    // sort(apx.begin(), apx.end());
    // for (NetworKit::count i = 0; i < k; ++i) {
    //     INFO(apx[i]);
    // }
    // INFO("-----------------------------");
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
    ifs.close();
    N++;
    NetworKit::Graph g(N, false, false);
    for(auto pr:belongs){
        g.addEdge(pr.first, pr.second);
    }
    INFO("N=", N);
    filename=argv[1];
    int lf=filename.find_last_of('/');
    int rg=filename.find_last_of('_');
    string outf= filename.substr(lf+1, rg-2);
    outf.append("_");
    outf.append(argv[2]);
    outf.append(".txt");
    omp_set_num_threads(1);
    SkylineGraph sg(g);
    sg.ReadGraphFile(argv[1]);
    Aux::Timer t;
    t.start();
    if(strcmp(argv[2],"basesky")==0)
      sg.ImprNoEdgeSky(g, skyvec);
    if(strcmp(argv[2],"filterrefinesky")==0)
      sg.PrunedCnbrNoEdgeSky4(skyvec);
    t.stop();
    INFO("skyline runtime=", t.elapsedMilliseconds(), "ms");
    OutputRes(outf.data());
 	return 0;

}