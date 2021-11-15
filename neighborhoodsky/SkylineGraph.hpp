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


#include <bitset> //输出二进制的头文件

#define mp(a,b) make_pair(a,b)

using namespace std;

int BK;
int BD;

vector<NetworKit::node> Deg;
vector<NetworKit::node> Dlt;

struct edge{
    NetworKit::node s;
    NetworKit::node t;
};

int cmpdeg(NetworKit::node u, NetworKit::node v)
{
	if(Deg[u]==Deg[v]){
		return u<v;
	}
	return Deg[u]>Deg[v];
}

struct bfnode {
    int* bf_arr;
#if BK > 8
    unsigned int h_n;
#else
    unsigned char h_n;
#endif
    int* bf_2hop;
    bool flag=false;
};



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

    vector<vector<NetworKit::node> > Edge_map;
    vector<vector<NetworKit::node> > Edge_map_directed;
    NetworKit::count GN;
    NetworKit::count GE;
    vector<edge> DirEdge;

    vector<NetworKit::node> vsup;
    map<pair<NetworKit::node, NetworKit::node>, int> emap;

    vector<vector<NetworKit::node> > Edge_flag;

    //bloom filter
    vector<bfnode> bfnodes;
    
    NetworKit::count max_deg;

    vector<vector<NetworKit::node> > Attr_map;


public:
    SkylineGraph(const NetworKit::Graph &G);

    void ImprNoEdgeSky(const NetworKit::Graph &G, vector<NetworKit::node> &skyres);

    void WithEdgeSky(vector<NetworKit::node> &skyres);
    void PrunedCnbrNoEdgeSky4(vector<NetworKit::node> &skyres);

    void ReadGraphFile(const char*);//file name ,mode

    int BinarySearch(NetworKit::node e, NetworKit::node s);
    int BinarySearchLoc(NetworKit::node t, NetworKit::node s);

	int h();
  
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


void SkylineGraph::ReadGraphFile(const char* filename){
    ifstream ifs;
    ifs.open(filename, ifstream::in);
    if(!ifs.is_open()){
        cout<<"error opening file\n";
    }
    NetworKit::node a, b;
    vector<pair<NetworKit::node, NetworKit::node>>belongs;
    GN=0;
    GE=0;
    while(ifs>>a>>b){
        belongs.push_back(mp(a,b));
        GE++;
        if(a>GN) GN=a;
        if(b>GN) GN=b;
    }
    ifs.close();

    GN++;
    GE++;
    Edge_map.resize(GN);
    for(auto pr:belongs){
        Edge_map[pr.first].push_back(pr.second);
        Edge_map[pr.second].push_back(pr.first);
    }

    Attr_map.resize(GN);
    for(int i=0; i<GN; i++){
        Attr_map[i]=Edge_map[i];
        Attr_map[i].push_back(i);
        sort(Attr_map[i].begin(), Attr_map[i].end());
        //cout<<"deg "<<i<<" = "<<Deg[i]<<endl;
    }

    Edge_flag.resize(GN);
    for(int i=0; i<GN; i++){
        Edge_flag[i].resize(Edge_map[i].size(), 1);
        //cout<<"deg "<<i<<" = "<<Deg[i]<<endl;
    }
    Deg.resize(GN, 0);
    max_deg=0;
    for(int i=0; i<GN; i++){
        Deg[i]=Edge_map[i].size();
        if(max_deg<Deg[i]){
            max_deg=Deg[i];
        }
        //cout<<"deg "<<i<<" = "<<Deg[i]<<endl;
    }
    INFO("N=", GN);
    INFO("max_deg=", max_deg);
}

int SkylineGraph::BinarySearch(NetworKit::node u1, NetworKit::node u2){
    
    NetworKit::node s, t;
    if(Deg[u1]<=Deg[u2]){
        s=u1;
        t=u2;
    }else{
        s=u2;
        t=u1;
    }

    int st=0;
    int ed=Edge_map[s].size()-1;
    int mid;
    while(st<=ed){
        mid=st+(ed-st)/2;
        if(Edge_map[s][mid]<t)
            st=mid+1;
        else if(Edge_map[s][mid]>t)
            ed=mid-1;
        else
            return mid;
    }
    return -1;
}

int SkylineGraph::BinarySearchLoc(NetworKit::node t, NetworKit::node s){

    int st=0;
    int ed=Edge_map[s].size()-1;
    int mid;
    while(st<=ed){
        mid=st+(ed-st)/2;
        if(Edge_map[s][mid]<t)
            st=mid+1;
        else if(Edge_map[s][mid]>t)
            ed=mid-1;
        else
            return mid;
    }
    return -1;
}

void SkylineGraph::ImprNoEdgeSky(const NetworKit::Graph &G, vector<NetworKit::node> &skyres){

    NetworKit::count n=G.upperNodeIdBound();
    skyres.resize(n);
    order.resize(n);
    T.resize(n, 0);

    for(NetworKit::node i=0; i<n; i++){
        order[i]=i;
        skyres[i]=0;
    } 
    //cout<<"init ok"<<endl;

    for(NetworKit::node u=0; u<n; u++){
        if(u!=order[u]){
        //     //cout<<u<<"can skip"<<endl;
            continue;
        }
        for(NetworKit::node i=0; i<n; i++){
            T[i]=0;
        }
        bool flag=true;

        for(NetworKit::node vid=0; vid<Edge_map[u].size(); vid++){
			NetworKit::node v=Edge_map[u][vid];
			T[v]++;
            if(T[v]==G.degree(u)){
                if(G.degree(v)==G.degree(u)){
                    // int min=u<v?u:v;
                    // int max=u<v?v:u;
                    // if(order[max]==max){
                    //     order[max]=min;
                    // 	cout<<u<<" is equal to "<<v<<endl;
                    // 	cout<<"max = "<<max<<", min = "<<min<<endl;
                    // }
                    if(order[v]==v)
                    {
                        //cout<<v<<" is equal by "<<u<<endl;
                        order[v]=u;
                    }
                }else{
                    if(order[u]==u)
                    {
                        //cout<<u<<" is dominated by "<<v<<endl;
                        order[u]=v;
                        break;
                    }
                }
                // break;  
            }
			for(NetworKit::node wid=0; wid<Edge_map[v].size(); wid++){
				NetworKit::node w=Edge_map[v][wid];
				if(w!=u)
                {
                    T[w]++;
                    if(T[w]==G.degree(u)){
                        if(G.degree(w)==G.degree(u)){
                       //      int min=u<w?u:w;
                       //      int max=u<w?w:u;
                       //      if(order[max]==max){
                       //          order[max]=min;
		                    	// cout<<u<<" is equal to "<<w<<endl;
		                    	// cout<<"max = "<<max<<", min = "<<min<<endl;
                       //      }
                            if(order[w]==w)
		                    {
		                        //cout<<w<<" is equal by "<<u<<endl;
		                        order[w]=u; 

		                    }

                            //cout<<u<<" is equal to "<<wid<<endl;
                        }else{
                            if(order[u]==u)
                            {
                                //cout<<u<<" is dominated by "<<w<<endl;
                                order[u]=w;
                                flag=false;
                        		break; 
                               
                            }
                        }
                          
                    }
                }
			}
			if(!flag){
				break;
			}

		}
        // if(u==order[u]){
        //     skyres[u]=1;
        //     //cout<<u<<" is a skyline"<<endl;
        // }
        // else{
        //  cout<<u<<" is dominated by "<<order[u]<<endl;
        // }
    }
    for(NetworKit::node u=0; u<GN; u++){
    	if(u==order[u]){
            skyres[u]=1;
            // cout<<u<<" is a skyline"<<endl;
        }else{
        	skyres[u]=0;
            // cout<<u<<" is dominated by "<<order[u]<<endl;
        }
    }
    return;
}

void SkylineGraph::WithEdgeSky(vector<NetworKit::node> &skyres){
    skyres.resize(GN);
    order.resize(GN);
    int cnt=0;
    for(NetworKit::node i=0; i<GN; i++){
        order[i]=i;
        skyres[i]=0;
    } 
    for(NetworKit::node u=0; u<GN; u++){
        if(u!=order[u]){
            continue;
        }
        for(NetworKit::node vid=0; vid<Edge_map[u].size(); vid++){
			NetworKit::node v=Edge_map[u][vid];
			if(Deg[v]>=Deg[u] && order[v]==v){
				bool flag=true;
				for(NetworKit::node wid=0; wid<Edge_map[u].size(); wid++){
					NetworKit::node w=Edge_map[u][wid];
					if(w!=v && BinarySearch(v, w)==-1){
                        //cout<<"false"<<endl;
                        flag=false;
                        break;
                    }
				}
				if(flag){
                    if(Deg[v]==Deg[u]){
                        int min=u<v?u:v;
                        int max=u<v?v:u;
                        if(order[max]==max){
                            order[max]=min;
                        }
                        //cout<<u<<" is equal to "<<wid<<endl;
                    }else{
                        if(order[u]==u)
                        {
                            //cout<<u<<" is dominated by "<<wid<<endl;
                            order[u]=v;
                        }
                    }   
                }
			}
		}
        if(u==order[u]){
            skyres[u]=1;
            cnt++;
            // cout<<u<<" is a skyline"<<endl;
        }
        // else{
        // 	skyres[u]=0;
        //  	// cout<<u<<" is dominated by "<<order[u]<<endl;
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
    //cout<<"cnt1="<<cnt<<endl;
    // cnt=0;
    // for(NetworKit::node u=0; u<GN; u++){
    // 	if(u==order[u]){
    //         // cnt++;
    //         cout<<u<<" is a skyline"<<endl;
    //     }
    //     else{
    //         cout<<u<<" is dominated by "<< order[u] <<endl;
    //     }
    // }
    // cout<<"cnt2="<<cnt<<endl;

    return;

}

void SkylineGraph::PrunedCnbrNoEdgeSky4(vector<NetworKit::node> &skyres){


	if(max_deg%32==0){
        BK=max_deg/32;
    }else{
        BK=max_deg/32+1;
    }

    BD=320*BK;

    bfnodes.resize(GN);

    //1-hop bloom filter
    for(NetworKit::node u=0; u<GN; u++){
        if (bfnodes[u].flag == false) {
            bfnodes[u].h_n = h() % (BK * 32);
            bfnodes[u].flag = true;
        }
        bfnodes[u].bf_arr=(int *)malloc(sizeof(int)*BK);
        //cout<<"hash "<<u<<" = "<< hex << bfnodes[u].h_n<<endl;
        //cout<<bitset<sizeof(bfnodes[u].h_n)*8>(bfnodes[u].h_n)<<endl;
        for (int i = 0; i < BK; i++) {
            bfnodes[u].bf_arr[i] = 0;
        }
        for(NetworKit::node vid=0; vid<Edge_map[u].size(); vid++){
            
            NetworKit::node v=Edge_map[u][vid];
            if (bfnodes[v].flag == false) {
                bfnodes[v].h_n = h() % (BK * 32);
                bfnodes[v].flag = true;
            }
            ///cout<<"nbr"<<v<<" = "<< hex << bfnodes[v].h_n<<endl;
            int hu = bfnodes[v].h_n;
            bfnodes[u].bf_arr[(hu >> 5) % BK] |= 1 << (hu & 31);
            //cout<<"set "<<u<<" byteloc: "<< (hu >> 5) % BK << ", bitloc: " << (hu & 31)<<endl;
        }
    }

    // for(NetworKit::node u=0; u<GN; u++){
    //     cout<<"hash "<<u<<" = "<< hex << bfnodes[u].h_n<<endl;
    //     cout<<bitset<sizeof(bfnodes[u].h_n)*8>(bfnodes[u].h_n)<<endl;
    //  //    for (int i = 0; i < BK; i++) {
    //  //         cout << bfnodes[u].bf_arr[i]<<endl;
    //  //         cout<<bitset<sizeof(bfnodes[u].bf_arr[i])*8>(bfnodes[u].bf_arr[i])<<endl;
    //     // }
    // }

    skyres.resize(GN);  
	WithEdgeSky(skyres);

	for(NetworKit::node u=0; u<GN; u++){
		//cout<<u<<"--"<<order[u]<<"--"<<skyres[u]<<endl;
        if(skyres[u]!=1){
            //cout<<u<<"can skip"<<endl;
            continue;
        }
        bool flag=true;
        for(NetworKit::node vid=0; vid<Edge_map[u].size() && flag; vid++){
			NetworKit::node v=Edge_map[u][vid];
			for(NetworKit::node wid=0; wid<Edge_map[v].size() && flag; wid++){
				NetworKit::node w=Edge_map[v][wid];
				if(w!=u && Deg[w]>=Deg[u] && skyres[w]==1){
					bool cmpflag=true;
		            for(NetworKit::node cid=0; cid<Edge_map[u].size(); cid++){
		                NetworKit::node c=Edge_map[u][cid];
		                if((bfnodes[w].bf_arr[(bfnodes[c].h_n >> 5) % BK] & (1 << (bfnodes[c].h_n & 31))) == 0){
		                    cmpflag=false;
		                    break;
		                }else{
		                	if(c!=v && BinarySearch(w, c)==-1){
	                        	//cout<<"false"<<endl;
		                        flag=false;
		                        break;
		                    }
	                    }
		            }
		            if(!cmpflag || !flag){
		            	flag=true;
		                continue;
		            }
					if(flag){
	                    if(Deg[w]==Deg[u]){
	                        int min=u<w?u:w;
	                        int max=u<w?w:u;
	                        if(skyres[max]==1){
	                            skyres[max]=0;
	                        }
                            flag=true;
	                        //cout<<u<<" is equal to "<<wid<<endl;
	                    }else{
	                        if(skyres[u]==1)
	                        {
	                            //cout<<u<<" is dominated by "<<wid<<endl;
	                            skyres[u]=0;
	                        }
                            flag=false;
	                    } 
	                }
				}
			}
		}
    }
    // for(NetworKit::node u=0; u<GN; u++){
    // 	if(u==order[u]){
    //         skyres[u]=1;
    //         //cout<<u<<" is a skyline"<<endl;
    //     }else{
    //     	skyres[u]=0;
    //     }
    // }
    return;
}



int SkylineGraph::h() {
  static int c = 0, r = rand();
  if (c >= (int) GN / BD) {
    c = 0;
    r = rand();
  }
  c++;
  return r;
}

#endif