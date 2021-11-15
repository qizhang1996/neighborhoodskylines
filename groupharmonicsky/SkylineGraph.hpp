#ifndef SKYLINEGRAPH_HPP
#define SKYLINEGRAPH_HPP

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <iostream>

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <map>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

using namespace std;

class SkylineGraph{
protected:
    NetworKit::Graph G;
    vector<NetworKit::node> order;
public:
    vector<NetworKit::node> rtnsky;

    vector<NetworKit::node> B;
    vector<NetworKit::node> T;
    vector<NetworKit::node> tmpnbr_u;
    vector<NetworKit::node> tmpnbr_v;

public:
    SkylineGraph(const NetworKit::Graph &G);
    void Skyline_GraphOrder(const NetworKit::Graph &G, vector<NetworKit::node> &skyres);

    // void OutputRes(){
    //     int cnt=0;
    //     printf("-------------results-------------\n");
    //     for(int i=0; i<skyres.size(); i++){
    //         if(skyres[i]==1){
    //             cnt++;
    //             cout<<i<<"\n";
    //         }
    //     }
    //     cout<<"skyres.size()="<<cnt<<endl;
    // }

    // void OutputRes(vector<int> R){
    //     printf("-------------results-------------\n");
    //     for(int i=0; i<R.size(); i++){
    //         if(R[i]==1)
    //             cout<<i<<"\n";
    //     }
    // }

    // void OutputRes(const char* str){
    //     ofstream ofs(str);
    //     int cnt=0;
    //     for(int i=0; i<skyres.size(); i++){
    //         if(skyres[i]==1){
    //             ofs<<i<<endl;
    //             cnt++;
    //         }
    //     }
    //     cout<<"skyres.size()="<<cnt<<endl;
    //     ofs.close();
    // }
};

SkylineGraph::SkylineGraph(const NetworKit::Graph &G): G(G){};

void SkylineGraph::Skyline_GraphOrder(const NetworKit::Graph &G, vector<NetworKit::node> &skyres){

    //omp_lock_t lock;
    //omp_init_lock(&lock);
    //omp_set_num_threads(1);
    //INFO("thread num = ", omp_get_num_threads());
    NetworKit::count n=G.upperNodeIdBound();
    skyres.resize(n);
    order.resize(n);
    B.resize(n, 0);
    T.resize(n, 0);
    //#pragma omp parallel for
    for(NetworKit::node i=0; i<n; i++){
        order[i]=i;
        skyres[i]=0;
        //INFO("thread num = ", omp_get_num_threads());
    } 
    //cout<<"init ok"<<endl;

    //#pragma omp parallel for
    for(NetworKit::node u=0; u<n; u++){
        if(u!=order[u]){
        //     //cout<<u<<"can skip"<<endl;
            continue;
        }
        for(NetworKit::node i=0; i<n; i++){
            B[i]=0;
            T[i]=0;
        //INFO("thread num = ", omp_get_num_threads());
        } 
        // B.clear();
        // T.clear();
        tmpnbr_u.clear();
        //bool flag=false;
        G.forNeighborsOf(u, [&](NetworKit::node v) {
            tmpnbr_u.push_back(v);
        });
        for(NetworKit::node v=0; v<tmpnbr_u.size(); v++){
            tmpnbr_v.clear();
            NetworKit::node vid=tmpnbr_u[v];
            tmpnbr_v.push_back(vid);
            G.forNeighborsOf(vid, [&](NetworKit::node w) {
                tmpnbr_v.push_back(w);
            });
            //sort(tmpnbr_v.begin(), tmpnbr_v.end());
            for(NetworKit::node w=0; w<tmpnbr_v.size(); w++){
                NetworKit::node wid=tmpnbr_v[w];
                if(wid==u){
                    continue;
                }
                if(B[wid]==0){
                    B[wid]=1;
                    T[wid]=0;
                }
                T[wid]++;
                if(T[wid]==G.degree(u)){
                    if(G.degree(wid)==G.degree(u)){
                        int min=u<wid?u:wid;
                        int max=u<wid?wid:u;
                        //omp_set_lock(&lock);
                        if(order[max]==max){
                            order[max]=min;
                        }
                        //omp_unset_lock(&lock);
                        //cout<<u<<" is equal to "<<wid<<endl;
                    }else{
                        if(order[u]==u)
                        {
                            //cout<<u<<" is dominated by "<<wid<<endl;
                            order[u]=wid;
                        }
                        //flag=true;
                        //break;
                    }   
                }
            }
        }
        if(u==order[u]){
            skyres[u]=1;
            //cout<<u<<" is a skyline"<<endl;
        }
        // else{
        //  cout<<u<<" is dominated by "<<order[u]<<endl;
        // }
    }

    // for(NetworKit::node i=0; i<n; i++){
    //     //cout<<i<<", cnt = "<< cnt[i] << ", fa = "<< fa[i] <<", order = " << order[i] <<endl;
    //     if(cnt[i]!=1 && fa[i]==i){
    //         //cout<<i<<" is skyline"<<endl;
    //         skyres[i]=1;
    //     }
    // }
    //OutputRes();
    return;
}

#endif