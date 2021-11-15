#ifndef IMPRGROUPCLOSENESS_HPP
#define IMPRGROUPCLOSENESS_HPP

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/centrality/TopCloseness.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <sys/time.h>

#include "SkylineGraph.hpp"

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <map>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include <numeric>
#include <sstream>
#include <atomic>
#include <memory>


using namespace std;

struct timeval tstart;
struct timeval tend;

class imprGroupCloseness{
protected:
    NetworKit::Graph G;
    NetworKit::count k = 1;
    vector<NetworKit::count> D;
    NetworKit::count iters;
    NetworKit::count maxD;
    vector<NetworKit::count> d;
    vector<NetworKit::count> d1;
    vector<NetworKit::node> S;
    NetworKit::count H = 0;

public:
    imprGroupCloseness(const NetworKit::Graph &G, NetworKit::count k, NetworKit::count H);
    void run(bool usesky, vector<NetworKit::node>& skyvec);
    vector<NetworKit::node> groupMaxCloseness();
    double computeFarness(vector<NetworKit::node> S, NetworKit::count H = numeric_limits<NetworKit::count>::max());
    double scoreOfGroup(const vector<NetworKit::node> &group) const;
    NetworKit::edgeweight computeImprovement(NetworKit::node u, NetworKit::count n, NetworKit::Graph &G, NetworKit::count h);
    vector<NetworKit::count> newDistances(NetworKit::node u, NetworKit::count n, NetworKit::Graph &G, NetworKit::count h);
    static bool pairCompare(const pair<NetworKit::node, NetworKit::count> &firstElem,
    const pair<NetworKit::node, NetworKit::count> &secondElem);
    void checkGroup(const vector<NetworKit::node> &group) const;

    void outputTime()
    {
        printf("exec_time=%lf(ms)\n", timearray);
    }

public:
    double timearray;
};

imprGroupCloseness::imprGroupCloseness(const NetworKit::Graph &G, NetworKit::count k = 1, NetworKit::count H = 0)
    : G(G), k(k), H(H) {}

void imprGroupCloseness::run(bool usesky, vector<NetworKit::node> &skyvec) {

    //INFO("thread num = ", omp_get_num_threads());
    //omp_set_num_threads(1);

    timearray=0.0;
    gettimeofday(&tstart,NULL);

    NetworKit::count n = G.upperNodeIdBound();
    NetworKit::node top = 0;
    iters = 0;
    vector<bool> visited(n, false);
    vector<NetworKit::node> pred(n);
    vector<NetworKit::count> distances(n);
    D.clear();
    D.resize(n, 0);
    vector<pair<NetworKit::node, NetworKit::count>> degPerNode(n);

    // compute degrees per node
    G.parallelForNodes([&](NetworKit::node v) {
        D[v] = G.degree(v);
        degPerNode[v] = make_pair(v, D[v]);
    });
    omp_lock_t lock;
    omp_init_lock(&lock);
    // sort by degree (in descending order) and retrieve max and argmax
    std::sort(degPerNode.begin(), degPerNode.end(), pairCompare);
    NetworKit::node nodeMaxDeg = degPerNode[0].first;

    if (H == 0) {
        NetworKit::TopCloseness topcc(G, 1, true, false);
        topcc.run();
        top = topcc.topkNodesList()[0];
    } else {
        top = nodeMaxDeg;
    }

    // first, we store the distances between each node and the top node
    d.clear();
    d.resize(n);
    NetworKit::BFS bfs(G, top);
    bfs.run();

    G.parallelForNodes([&](NetworKit::node v) { d[v] = bfs.distance(v); });

    // get max distance
    maxD = 0;
    NetworKit::count sumD = 0;
    // TODO: actually, we could have more generic parallel reduction iterators in
    // the Graph class
    G.forNodes([&](NetworKit::node v) {
        if (d[v] > maxD) {
            maxD = d[v];
        }
        sumD += d[v];
    });
    // INFO("maxD = ", maxD);
    // INFO("top-1 = ", top);
    // INFO("sumD = ", sumD);
    // init S
    S.clear();
    S.resize(k, 0);
    S[0] = top;

    vector<int64_t> prevBound(n, 0);
    d1.clear();
    NetworKit::count currentImpr = sumD + 1; // TODO change
    NetworKit::count maxNode = 0;

    vector<NetworKit::count> S2(n, sumD);
    S2[top] = 0;
    vector<int64_t> prios(n);


    // vector<NetworKit::node> skyvec;

    // if(usesky){
    //     SkylineGraph sg(G);
    //     // Aux::Timer ts;
    //     // ts.start();
    //     skyvec=sg.Skyline_GraphOrder(G);
    //     // ts.stop();
    // }
    //INFO("runtime=", ts.elapsedMilliseconds(), "ms");

    // vector<NetworKit::node> rtnsky;
    // for(NetworKit::index i=0; i<n; i++){
    //     if(skyvec[i]==1){
    //         //cout<<i<<" is skyline"<<endl;
    //         rtnsky.push_back(i);
    //     }
    // }
    //cout<<"sky size = "<<rtnsky.size()<<endl;
    NetworKit::count actnum=0;
    //ts.start();
    // loop to find k group members
    for (NetworKit::index i = 1; i < k; i++) {
        //INFO("k = ", i);
        G.parallelForNodes([&](NetworKit::node v) { 
            prios[v] = -prevBound[v]; 
            //INFO("prios ", v, " = ", -prevBound[v]);
        });
        // Aux::BucketPQ Q(prios, currentImpr + 1);
        Aux::BucketPQ Q(n, -currentImpr - 1, 0);
        G.forNodes([&](NetworKit::node v) { Q.insert(prios[v], v); });
        currentImpr = 0;
        maxNode = 0;
        d1.resize(G.upperNodeIdBound());

        std::atomic<bool> toInterrupt{false};
#pragma omp parallel // Shared variables:
        // cc: synchronized write, read leads to a positive race condition;
        // Q: fully synchronized;
        {
            //INFO("thread num = ", omp_get_num_threads());
            while (!toInterrupt.load(std::memory_order_relaxed)) {
                omp_set_lock(&lock);
                if (Q.empty()) {
                    omp_unset_lock(&lock);
                    toInterrupt.store(true, std::memory_order_relaxed);
                    break;
                }

                auto topPair = Q.extractMin();
                NetworKit::node v = topPair.second;
                omp_unset_lock(&lock);
                // INFO("Extracted node ", v, " with prio ", prevBound[v]);
                if (i > 1 && prevBound[v] <= static_cast<int64_t>(currentImpr)) {
                    //INFO("Interrupting! currentImpr = ", currentImpr,
                        //", previous bound = ", prevBound[v]);
                    toInterrupt.store(true, std::memory_order_relaxed);
                    break;
                }
                // if(skyvec[v]==0){
                //     //INFO("Interrupting! ", v,
                //          //" is not a skyline vertex with skyres = ", skyvec[v]);
                //     prevBound[v] = 0;
                //     continue;
                //     // toInterrupt.store(true, std::memory_order_relaxed);
                //     // break;
                // }
                bool vvalid=false;
                if(usesky && skyvec[v]>0){
                    vvalid=true;
                }else if(!usesky){
                    vvalid=true;
                }
                if(vvalid && D[v] > 1 && !(d[v] == 1 && D[v] == 2) && d[v] > 0 &&
                    (i == 1 || prevBound[v] > static_cast<int64_t>(currentImpr))) {
                    
                    NetworKit::count imp = computeImprovement(v, n, G, H);
                    omp_set_lock(&lock);
                    actnum++;
                    if (imp > currentImpr) {
                        currentImpr = imp;
                        maxNode = v;
                        //INFO("New currentImpr = ", imp);
                        //INFO("New maxNode = ", maxNode);
                    }
                    omp_unset_lock(&lock);
                    //INFO("New bound for v = ", imp);
                    prevBound[v] = imp; // TODO use later
                }else {
                    prevBound[v] = 0;
                }
            }
        }
        S[i] = maxNode;
        // INFO("choose ", maxNode, " as ", i, "-th element");

        d1 = newDistances(S[i], n, G, 0);
        G.parallelForNodes([&](NetworKit::node v) { d[v] = d1[v]; });
    }

    gettimeofday(&tend,NULL);
    double exec_time_stop = ((tend.tv_sec-tstart.tv_sec)*1000000+(tend.tv_usec-tstart.tv_usec))/1000;
    timearray=exec_time_stop;

    // ts.stop();
    // INFO("runtime=", ts.elapsedMilliseconds(), "ms");
    //INFO("actnum=", actnum);
}


vector<NetworKit::node> imprGroupCloseness::groupMaxCloseness(){
    return S;
}


double imprGroupCloseness::computeFarness(vector<NetworKit::node> S, NetworKit::count H) {
    // we run a BFS from S up to distance H (if H > 0) and sum the distances
    double farness = 0;
    NetworKit::count k = S.size();
    vector<double> d1(G.upperNodeIdBound(), 0);
    vector<bool> visited(G.upperNodeIdBound(), false);
    queue<NetworKit::node> Q;

    for (NetworKit::node i = 0; i < k; i++) {
        Q.push(S[i]);
        visited[S[i]] = true;
    }

    while (!Q.empty()) {
        NetworKit::node u = Q.front();
        Q.pop();
        if (H > 0 && d1[u] > H) {
            break;
        }
        farness += d1[u];
        G.forNeighborsOf(u, [&](NetworKit::node w) {
            if (!visited[w]) {
                visited[w] = true;
                d1[w] = d1[u] + 1;
                Q.push(w);
            }
        });
    }
    return farness;
}


double imprGroupCloseness::scoreOfGroup(const vector<NetworKit::node> &group) const {
    vector<bool> explored(G.upperNodeIdBound(), false);
    vector<NetworKit::count> distance(G.upperNodeIdBound(), 0);

    for (NetworKit::count i = 0; i < group.size(); ++i) {
        explored[group[i]] = true;
    }

    vector<NetworKit::node> queue;
    auto exploreNode = [&](NetworKit::node w, NetworKit::count d) {
        explored[w] = true;
        queue.push_back(w);
        distance[w] = d;
    };

    NetworKit::count d = 1;
    for (auto u : group) {
        G.forNeighborsOf(u, [&](NetworKit::node v) {
            if (!explored[v]) {
                exploreNode(v, d);
            }
        });
    }

    while (queue.size() > 0) {
        ++d;
        NetworKit::node u = queue.front();
        queue.erase(queue.begin());
        G.forNeighborsOf(u, [&](NetworKit::node v) {
            if (!explored[v]) {
                exploreNode(v, d);
            }
        });
    }

    double dSum = std::accumulate(distance.begin(), distance.end(), 0.0);
    return dSum == 0
               ? 0.
               : ((double)G.upperNodeIdBound() - (double)group.size()) / dSum;
}


NetworKit::edgeweight imprGroupCloseness::computeImprovement(NetworKit::node u, NetworKit::count n, NetworKit::Graph &G, NetworKit::count h) {
    // computes the marginal gain due to adding u to S
    vector<NetworKit::count> d1(n);
    G.forNodes([&](NetworKit::node v) { d1[v] = d[v]; });

    d1[u] = 0;
    NetworKit::count improvement = d[u]; // old distance of u
    queue<NetworKit::node> Q;
    Q.push(u);

    NetworKit::index level = 0;
    while (!Q.empty() && (h == 0 || level <= h)) {
        NetworKit::node v = Q.front();
        Q.pop();
        level = d1[v];
        G.forNeighborsOf(v, [&](NetworKit::node w) {
            if (d1[w] > d1[v] + 1) {
                d1[w] = d1[v] + 1;
                improvement += d[w] - d1[w];
                Q.push(w);
            }
        });
    }
    return improvement;
}

vector<NetworKit::count> imprGroupCloseness::newDistances(NetworKit::node u, NetworKit::count n, NetworKit::Graph &G, NetworKit::count h) {
    vector<NetworKit::count> d1(n);
    G.forNodes([&](NetworKit::node v) { d1[v] = d[v]; });

    d1[u] = 0;
    NetworKit::count improvement = d[u]; // old distance of u
    queue<NetworKit::node> Q;
    Q.push(u);

    NetworKit::index level = 0;
    while (!Q.empty() && (h == 0 || level <= h)) {
        NetworKit::node v = Q.front();
        Q.pop();
        level = d1[v];
        G.forNeighborsOf(v, [&](NetworKit::node w) {
            if (d1[w] > d1[v] + 1) {
                d1[w] = d1[v] + 1;
                improvement += d[w] - d1[w];
                Q.push(w);
            }
        });
    }
    return d1;
}

bool imprGroupCloseness::pairCompare(const pair<NetworKit::node, NetworKit::count> &firstElem,
                 const pair<NetworKit::node, NetworKit::count> &secondElem) {
    return firstElem.second > secondElem.second;
}



void imprGroupCloseness::checkGroup(const vector<NetworKit::node> &group) const {
    const NetworKit::count z = G.upperNodeIdBound();
    vector<bool> check(z, false);
#pragma omp parallel for
    for (NetworKit::omp_index i = 0; i < static_cast<NetworKit::omp_index>(group.size()); ++i) {
        NetworKit::node u = group[i];
        if (u >= z) {
            stringstream ss;
            ss << "Error: node " << u << " is not in the graph.\n";
            throw std::runtime_error(ss.str());
        }
        if (check[u]) {
            stringstream ss;
            ss << "Error: the group contains duplicates of node " << u << ".\n";
            throw std::runtime_error(ss.str());
        }
        check[u] = true;
    }
}

#endif