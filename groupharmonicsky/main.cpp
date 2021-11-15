#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <iostream>
#include<fstream>
#include <omp.h>

#include "SkylineGraph.hpp"
#include "imprGroupHarmonic.hpp"

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
        cout<<"error opening file\n";
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
    imprGroupHarmonic ghc(g, k);
    ghc.run(usesky, skyvec);
    auto group = ghc.groupMaxHarmonicCloseness();
    sort(group.begin(), group.end());
    for (NetworKit::count i = 0; i < k; ++i) {
        INFO(group[i]);
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

    const auto computeOpt = [&](const NetworKit::Graph &G, NetworKit::count k) -> double {
        std::vector<bool> inGroup(G.upperNodeIdBound());
        std::fill(inGroup.begin(), inGroup.begin() + k, true);
        double opt = -std::numeric_limits<double>::max();
        std::vector<NetworKit::node> group;
        group.reserve(k);
        do {
            group.clear();
            G.forNodes([&](NetworKit::node u) {
                if (inGroup[u])
                    group.push_back(u);
            });
            opt = std::max(opt, imprGroupHarmonic::scoreOfGroup(G, group.begin(), group.end()));
        } while (std::prev_permutation(inGroup.begin(), inGroup.end()));
        return opt;
    };
    runalg(g);
 	return 0;

}