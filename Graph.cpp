#include "Utility.h"
#include "Graph.h"

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = m = 0;

}
Graph::Graph(){
	n=m=0;
}

Graph::~Graph() {
}

void Graph::add_edge(ui v,ui u){
	adj[v].emplace_back(u);
	adj[u].emplace_back(v);
}
void Graph::remove_edge(ui v,ui u){
	adj[v].remove(u);
	adj[u].remove(v);
}
void Graph::add_soft_edge(ui v,ui u,int weight){
	generalized_adj[v].emplace_back(make_pair(u,weight));
	generalized_adj[u].emplace_back(make_pair(v,weight));
	dynamic_weight[v]+=max(0,-weight);
	dynamic_weight[u]+=max(0,-weight);
	// all_soft_edges.emplace_back(make_pair(u,v));
}
void Graph::remove_soft_edge(ui v,ui u,int weight){
	generalized_adj[v].remove(make_pair(u,weight));
	generalized_adj[u].remove(make_pair(v,weight));
	dynamic_weight[v]-=max(0,-weight);
	dynamic_weight[u]-=max(0,-weight);
}
void Graph::read_graph_GIS() {
	FILE *f = open_file(dir.c_str(), "rb");
    // Skip lines starting with 'c'
    char line[1024];
    char first_char;
    do {
        if (fscanf(f, "%c", &first_char) != 1) break;
        
        if (first_char == 'c') {
            char* tmp=fgets(line, sizeof(line), f);
        } else {
            fseek(f, -1, SEEK_CUR);
            break;
        }
    } while (true);
    char word1[100], word2[100];
    int b=fscanf(f, "%s %s %d %d %d", word1, word2, &n, &m, &generalized_m);
    printf("vertices: %d, Hard edge: %d, Soft edge: %d\n", n, m, generalized_m);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("WA input graph\n");
#endif


	adj.resize(n);
	generalized_adj.resize(n);
	hard_degreeInC.resize(n,0);
	soft_degreeInC.resize(n,0);
	weight_node.resize(n,0);
	dynamic_weight.resize(n,0);
	Cand_reduce.reserve(n);
    char prefix[10];
    int u, v;
	int num_ep=0;
	int tmp;
	int min_hard=INT32_MAX;
	int min_soft=INT32_MAX;
	int max_hard=0;
	int max_soft=0;
    for (int i = 0; i < m; ++i) {
        tmp=fscanf(f, "%s %d %d", prefix, &u, &v);
		u--;
		v--;
		hard_degreeInC[u]++;
		hard_degreeInC[v]++;
		add_edge(u,v);
    }
	int w;
    for (int i = 0; i < n; ++i) {
        tmp=fscanf(f, "%s %d %d", prefix, &u, &w);
			u--;
            weight_node[u]=w;
			dynamic_weight[u]=w;
    }
    for (int i = 0; i < generalized_m; ++i) {
        tmp=fscanf(f, "%s %d %d %d", prefix, &u, &v,&w);
		u--;
		v--;
		soft_degreeInC[u]++;
		soft_degreeInC[v]++;
		add_soft_edge(u,v,w);
    }
	fclose(f);
	soft_neighbors.init(n);
	hard_neighbors.init(n);
	hard_degree_ones.init(n);
	hard_degree_twos.init(n);
	total_degree_ones.init(n);
	total_degree_twos.init(n);
	temp_delete_hard_ones.init(n);
	temp_delete_hard_twos.init(n);
	temp_delete_total_ones.init(n);
	temp_delete_total_twos.init(n);
	S.init(n);
}
void Graph::soft_degree_decreaseOne(char* is, ui v, int &res){
	if(is[v]==0)  return;
	soft_degreeInC[v]--;
	if(is[v]==2 || is[v]==3) return;
	if(soft_degreeInC[v]==2 && hard_degreeInC[v]==0){
		temp_add_total_twos.emplace_back(v);
	}else if(soft_degreeInC[v]==1 && hard_degreeInC[v]==1){
		temp_add_total_twos.emplace_back(v);
	}else if(soft_degreeInC[v]==1 && hard_degreeInC[v]==0){
		temp_delete_total_twos.add(v);
		temp_add_total_ones.emplace_back(v);
	}else if(soft_degreeInC[v]==0 && hard_degreeInC[v]==2){
		temp_add_total_twos.emplace_back(v);
	}else if(soft_degreeInC[v]==0 && hard_degreeInC[v]==1){
		temp_delete_total_twos.add(v);
		temp_add_total_ones.emplace_back(v);
	}else if(soft_degreeInC[v]==0 && hard_degreeInC[v]==0){
			d0.emplace_back(v);
	}
}
void Graph::add_0d(char* is,int& res){
	for(ui u:d0){
		if(is[u]!=1) continue;
		if(weight_node[u]>0)
		add_vertex(u,res,is);
		else
		delete_vertex(u,res,is);
	}
	d0.clear();
}
void Graph::soft_degree_increaseOne(char* is, ui v, int &res){
	if(is[v]==0) return;
	soft_degreeInC[v]++;
	if(is[v]==2 ||is[v]==3 ) return;
	if(soft_degreeInC[v]==2 && hard_degreeInC[v]==0){
		temp_delete_total_ones.add(v);
		temp_add_total_twos.emplace_back(v);
	}else if(soft_degreeInC[v]==1 && hard_degreeInC[v]==1){
		temp_delete_total_ones.add(v);
		temp_add_total_twos.emplace_back(v);
	}else if(soft_degreeInC[v]==2 && hard_degreeInC[v]==1){
		temp_delete_total_twos.add(v);
	}else if(soft_degreeInC[v]==1 && hard_degreeInC[v]==2){
		temp_delete_total_twos.add(v);
	}else if(soft_degreeInC[v]==3 && hard_degreeInC[v]==0){
		temp_delete_total_twos.add(v);
	}
}
void Graph::hard_degree_increaseOne(char* is, ui v, int &res){
	if(is[v]==0)  return;
	hard_degreeInC[v]++;
		if(is[v]==2 || is[v]==3) return;
	if(hard_degreeInC[v]==2){
		temp_delete_hard_ones.add(v);
		temp_add_hard_twos.emplace_back(v);
		if(soft_degreeInC[v]==0){
		temp_delete_total_ones.add(v);
		temp_add_total_twos.emplace_back(v);
		}else if(soft_degreeInC[v]==1){
			temp_delete_total_twos.add(v);
		}
	}else if(hard_degreeInC[v]==1){
		temp_add_hard_ones.emplace_back(v);
		if(soft_degreeInC[v]==1){
			temp_delete_total_ones.add(v);
			temp_add_total_twos.emplace_back(v);
		}else if(soft_degreeInC[v]==2){
			temp_delete_total_twos.add(v);
		}
	}else if(hard_degreeInC[v]==3){
		temp_delete_hard_twos.add(v);
		if(soft_degreeInC[v]==0){
			temp_delete_total_twos.add(v);
		}
	}
}
void Graph::hard_degree_decreaseOne(char* is, ui v, int &res){
	if(is[v]==0)  return;
	hard_degreeInC[v]--;
		if(is[v]==2 || is[v]==3) return;
	if(hard_degreeInC[v]==2){
		temp_add_hard_twos.emplace_back(v);
		if(soft_degreeInC[v]==0){
		temp_add_total_twos.emplace_back(v);	
		}
	}else if(hard_degreeInC[v]==1){
		temp_delete_hard_twos.add(v);
		temp_add_hard_ones.emplace_back(v);
		if(soft_degreeInC[v]==1){
			temp_add_total_twos.emplace_back(v);
		}else if(soft_degreeInC[v]==0){
			temp_delete_total_twos.add(v);
			temp_add_total_ones.emplace_back(v);
		}
	}else if(hard_degreeInC[v]==0){
		temp_delete_hard_ones.add(v);
		if(soft_degreeInC[v]==2){
			temp_add_total_twos.emplace_back(v);
		}else if(soft_degreeInC[v]==1){
			temp_delete_total_twos.add(v);
			temp_add_total_ones.emplace_back(v);
		}else if(soft_degreeInC[v]==0){
			d0.emplace_back(v);
		}
	}

}
void Graph::add_vertex4fold(ui u,int& res,char* is){
		is[u]=2;
		vector<pair<ui,ui>> n2r;
		for(auto& v:adj[u]){
			if(is[v]==1){
				is[v]=0;
				vector<pair<ui,int>> need2remove;
				for(auto& w:generalized_adj[v]){
					if(is[w.first]!=0){
						need2remove.emplace_back(w);
					}
				}
				for(auto& w:need2remove){
					soft_degree_decreaseOne(is,w.first,res);
					remove_soft_edge(w.first,v,w.second);
				}
				for(ui w:adj[v]){
					if(is[w]!=0){
						n2r.emplace_back(make_pair(w,v));
						hard_degree_decreaseOne(is,v,res);
					}
				}
			}
		}
		for(auto& w:n2r){
			remove_edge(w.first,w.second);
			hard_degree_decreaseOne(is,w.first,res);
		}
}
int Graph::local_search(int old_res,char* is,vector<ui>& k2e,vector<ui>& e2k, double& best_time,int& best_size){
	int run_time_seconds = 30;
	int best_res=old_res;
	char* old_is=new char[n];
	int max_depth=0;
	int res=old_res;
	int flag=0;
	Graph* new_kernel02=copy_kernel();
	for(int i=0;i<S.vnum;i++){
		new_kernel02->S.add(S.vlist[i]);
	}
	for(auto& myturple:exclusive_pairs){
		new_kernel02->exclusive_pairs.emplace_back(myturple);
	}
	auto start_time = std::chrono::high_resolution_clock::now();
	int num=0;
	while(true){
		num++;
		for(int i=0;i<n;i++){
			old_is[i]=is[i];
		}
		new_kernel02->S.clear();
		for(int i=0;i<S.vnum;i++){
			new_kernel02->S.add(S.vlist[i]);
		}
		res=old_res;
		vector<ui> c2k,k2c;
		Graph* search_kernel=copy_Cand(is,res,c2k,k2c);
		if(c2k.size()!=k2c.size() && search_kernel->n <100000){
			vector<ui> rc2k,k2rc;
			Graph* Cand_Copy=copy_Cand(is,res,rc2k,k2rc);
			Cand_Copy->Random_Peeling(old_is,res);
			for(int i=0;i<Cand_Copy->S.vnum;i++){
				search_kernel->S.add(Cand_Copy->S.vlist[i]);
			}
		}
		ui depth;
		depth=(int)n*10;
		bool need2reduction=false;
		need2reduction=true;
		int true_depth=search_kernel->origin_search(res,k2c,depth,start_time,run_time_seconds,best_time,best_res,need2reduction);
		max_depth=max(true_depth,max_depth);
		for(int i=0;i<search_kernel->S.vnum;i++){
			ui u=search_kernel->S.vlist[i];
			new_kernel02->S.add(c2k[u]);
		}
		if(!new_kernel02->exclusive_pairs.empty())
		new_kernel02->unfold();
		if(res>best_res){
			best_res=res;
			best_size=S.vnum;
		}
		auto current_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_time = current_time - start_time;
		if (elapsed_time.count() >= run_time_seconds) {
			cout<<"cost time: "<<elapsed_time.count()<<endl;
			break; 
		}
	}
	cout<<"res:"<<best_res<<" search_time: "<<best_time<<" size: "<<best_size<<endl;
	cout<<"search times: "<<num<<endl;
	return best_res;
}
void Graph::add_vertex(ui u,int& res,char* is){
				if(is[u]!=1) return; 
				res+=weight_node[u];
				is[u]=2;
				S.add(u);
				vector<pair<ui,int>> removed_edges;
				for(auto& v:generalized_adj[u]){
					if(is[v.first]==1){
						weight_node[v.first]-=v.second;
						dynamic_weight[v.first]-=v.second;
						removed_edges.emplace_back(v);
						soft_degree_decreaseOne(is,u,res);
					}
				}
				for(auto& v:removed_edges){
					soft_degree_decreaseOne(is,v.first,res);
					remove_soft_edge(u,v.first,v.second);
				}
				vector<pair<ui,ui>> n2r;
				for(auto& v:adj[u]){
					if(is[v]==1){
						is[v]=0;
						vector<pair<ui,int>> need2remove;
						for(auto& w:generalized_adj[v]){
							//
							if(is[w.first]==1){
								need2remove.emplace_back(w);
							}
						}
						for(auto& w:need2remove){
							bool f=false;
							int old_dy=dynamic_weight[w.first];
							int old_wei=weight_node[w.first];
							if(dynamic_weight[w.first]>=weight_node[w.first]) f=true;
							soft_degree_decreaseOne(is,w.first,res);
							remove_soft_edge(w.first,v,w.second);
						}
						for(ui w:adj[v]){
							if(is[w]==1){
								n2r.emplace_back(make_pair(w,v));
								hard_degree_decreaseOne(is,v,res);
							}
						}
					}
				}
				for(auto& w:n2r){
					remove_edge(w.first,w.second);
					hard_degree_decreaseOne(is,w.first,res);
				}
}
int Graph::update_penalty(ui u,ui v, int add_num,char* is,int& res){
	int puv=0;
	bool has_edge=false;
	if(add_num==0) return 0;
	for(auto& w:generalized_adj[u]){
		if(w.first==v){
			puv=w.second;
			has_edge=true;
			break;
		}
	}
	if(has_edge){
		remove_soft_edge(u,v,puv);
		if(puv+add_num==0){
			soft_degree_decreaseOne(is,u,res);
			soft_degree_decreaseOne(is,v,res);
			return puv;
		}
		add_soft_edge(u,v,puv+add_num);
	}
	else{
		puv=add_num;
		add_soft_edge(u,v,add_num);
		soft_degree_increaseOne(is,u,res);
		soft_degree_increaseOne(is,v,res);
	}
	return puv;
}
void Graph::delete4fold(ui u,int& res,char* is){
		if(is[u]!=1) return;
		is[u]=3;
		for(auto& v:generalized_adj[u]){
			if(is[v.first]!=0){
				soft_degree_decreaseOne(is,v.first,res);
				dynamic_weight[u]-=max(0,-v.second);
				dynamic_weight[v.first]-=max(0,-v.second);
			}
		}
		for(ui v:adj[u]){
			if(is[v]!=0){
				hard_degree_decreaseOne(is,v,res);
			}
		}
}
void Graph::delete_vertex(ui u,int& res,char* is){
		if(is[u]!=1) return;
		is[u]=0;
		vector<pair<ui,int>> need2remove;
		for(auto& v:generalized_adj[u]){
			if(is[v.first]==1){
				need2remove.emplace_back(v);
			}
		}
		for(auto& v:need2remove){
			soft_degree_decreaseOne(is,v.first,res);
			remove_soft_edge(v.first,u,v.second);
		}
		need2remove.clear();
		for(ui v:adj[u]){
			if(is[v]==1){
				need2remove.emplace_back(make_pair(v,0));
			}
		}
		for(auto& v:need2remove){
			remove_edge(v.first,u);
			hard_degree_decreaseOne(is,v.first,res);
		}
}
void Graph::soft2hard(char* is,int& res){
	int num=0;
	vector<pair<ui,int>> need2trans;
	for(int i=0;i<n;i++){
		if(is[i]!=1) continue;
		need2trans.clear();
		bool need2delete=false;
		for(auto& v:generalized_adj[i]){
			if(is[v.first]!=1 || v.first<=i || v.second<min(dynamic_weight[i],dynamic_weight[v.first])) continue;
			num++;
			if(is[v.first]==2) {
				need2delete=true;
				break;
				}
			need2trans.emplace_back(v);
		}
		if(need2delete) {
			delete_vertex(i,res,is);
			continue;
			}
		for(auto& v:need2trans){
			if(is[i]==2){
				delete_vertex(v.first,res,is);
				continue;
			}
			add_edge(i,v.first);
			hard_degree_increaseOne(is,v.first,res);
			hard_degree_increaseOne(is,i,res);
			soft_degree_decreaseOne(is,v.first,res);
			soft_degree_decreaseOne(is,i,res);
			remove_soft_edge(v.first,i,v.second);
		}
	}
}
int Graph::update_edges(char* is, int &res,int& valid_n){
	ui ids_n = 0;
	soft2hard(is,res);
	ui new_m = 0;
	ui new_generalized_m = 0;
	ui k=0;
	valid_n=0;
	for(ui i=0;i<n;i++){
		if(is[i]==0) continue;
		if(is[i]==1) valid_n++;
		ids_n++;
	}
	return ids_n;
}
int Graph::update_edges(char* is, int &res){
	ui ids_n = 0;
	soft2hard(is,res);
	ui new_m = 0;
	ui new_generalized_m = 0;
	ui k=0;
	for(ui i=0;i<n;i++){
		if(is[i]==0) continue;
		ids_n++;
	}
	return ids_n;
}
void Graph::update_sets(char* is){
	for(ui u:temp_add_total_ones){
		if(is[u]==1 && hard_degreeInC[u]+soft_degreeInC[u]==1 && !total_degree_ones.contains(u))
		total_degree_ones.add(u);
	}
	for(ui u:temp_add_total_twos){
		if(is[u]==1 && hard_degreeInC[u]+soft_degreeInC[u]==2 && !total_degree_twos.contains(u))
		total_degree_twos.add(u);
	}
	for(ui u:temp_add_hard_ones){
		if(is[u]==1 && hard_degreeInC[u]==1 && !hard_degree_ones.contains(u))
		hard_degree_ones.add(u);
	}
	for(ui u:temp_add_total_twos){
		if(is[u]==1 && hard_degreeInC[u]==2 & !hard_degree_twos.contains(u))
		hard_degree_twos.add(u);
	}
	temp_add_hard_ones.clear();
	temp_add_hard_twos.clear();
	temp_add_total_ones.clear();
	temp_add_total_twos.clear();

	temp_delete_total_ones.clear();
	temp_delete_total_twos.clear();
	temp_delete_hard_ones.clear();
	temp_delete_hard_twos.clear();
	for(int i=total_degree_ones.vnum-1;i>=0;i--){
		ui u=total_degree_ones.vlist[i];
		if(is[u]!=1 || hard_degreeInC[u]+soft_degreeInC[u]!=1)
		total_degree_ones.remove(u);
	}
	for(int i=total_degree_twos.vnum-1;i>=0;i--){
		ui u=total_degree_twos.vlist[i];
		if(is[u]!=1 || hard_degreeInC[u]+soft_degreeInC[u]!=2)
		total_degree_twos.remove(u);
	}
	for(int i=hard_degree_ones.vnum-1;i>=0;i--){
		ui u=hard_degree_ones.vlist[i];
		if(is[u]!=1 ||hard_degreeInC[u]!=1)
		hard_degree_ones.remove(u);
	}
	for(int i=hard_degree_twos.vnum-1;i>=0;i--){
		ui u=hard_degree_twos.vlist[i];
		if(is[u]!=1 ||hard_degreeInC[u]!=2)
		hard_degree_twos.remove(u);
	}
}
void Graph::update_sets(){
	for(ui u:temp_add_total_ones){
		if(hard_degreeInC[u]+soft_degreeInC[u]==1 && !total_degree_ones.contains(u))
		total_degree_ones.add(u);
	}
	for(ui u:temp_add_total_twos){
		if(hard_degreeInC[u]+soft_degreeInC[u]==2 && !total_degree_twos.contains(u))
		total_degree_twos.add(u);
	}
	for(ui u:temp_add_hard_ones){
		if(hard_degreeInC[u]==1 && !hard_degree_ones.contains(u))
		hard_degree_ones.add(u);
	}
	for(ui u:temp_add_total_twos){
		if(hard_degreeInC[u]==2 & !hard_degree_twos.contains(u))
		hard_degree_twos.add(u);
	}
	temp_add_hard_ones.clear();
	temp_add_hard_twos.clear();
	temp_add_total_ones.clear();
	temp_add_total_twos.clear();
	for(int i=0;i<temp_delete_total_ones.vnum;i++){
		ui u=temp_delete_total_ones.vlist[i];
		if(hard_degreeInC[u]+soft_degreeInC[u]!=1 && total_degree_ones.contains(u))
		total_degree_ones.remove(u);
	}
	for(int i=0;i<temp_delete_total_twos.vnum;i++){
		ui u=temp_delete_total_twos.vlist[i];
		if(hard_degreeInC[u]+soft_degreeInC[u]!=2 && total_degree_twos.contains(u))
		total_degree_twos.remove(u);
	}
	for(int i=0;i<temp_delete_hard_ones.vnum;i++){
		ui u=temp_delete_hard_ones.vlist[i];
		if(hard_degreeInC[u]!=1 && hard_degree_ones.contains(u))
		hard_degree_ones.remove(u);
	}
	for(int i=0;i<temp_delete_hard_twos.vnum;i++){
		ui u=temp_delete_hard_twos.vlist[i];
		if(hard_degreeInC[u]!=2 && hard_degree_twos.contains(u))
		hard_degree_twos.remove(u);
	}
	temp_delete_total_ones.clear();
	temp_delete_total_twos.clear();
	temp_delete_hard_ones.clear();
	temp_delete_hard_twos.clear();
}
void Graph::reduction_neighborhood(char* is,int& res,bool init){
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui i:Cand_reduce){
			if(is[i]!=1) continue;
			int sum=0;
			for(auto& v :adj[i]){
				if(is[v]==1){
					sum+=max(0,dynamic_weight[v]);
					if(sum>weight_node[i]) break ;
				}
			}
			if(sum>weight_node[i]) continue;
			for(auto& v :generalized_adj[i]){
				if(is[v.first]==1){
					sum+=max(0,dynamic_weight[v.first]);
					if(sum>weight_node[i]) break ;
				}
			}
			if(sum>weight_node[i]) continue;
			if(sum<=weight_node[i]){
				add_vertex(i,res,is);
				for(auto& w:generalized_adj[i]){
					if(is[w.first]==1)
					temp_cand.emplace_back(w.first);
				}
			}
	}
	Cand_reduce.clear();
	Cand_reduce=temp_cand;
}
void Graph::has_hard(ui v,ui x,char* is,int puv,int &res_w){
	bool res=false;
	for(ui w:adj[v]){
		if(w==x){
			res=true;
			break;
		}
	}
	bool has=false;
	if(!res){
		int pe=0;
		for(auto& w:generalized_adj[v]){
			if(w.first==x){
				pe=w.second;
				has=true;
				break;
			}
		}
		if(!has){
			soft_degree_increaseOne(is,x,res_w);
			soft_degree_increaseOne(is,v,res_w);
			add_soft_edge(x,v,-puv);
		}else{
			remove_soft_edge(x,v,pe);
			if(pe-puv==0){
			soft_degree_decreaseOne(is,x,res_w);
			soft_degree_decreaseOne(is,v,res_w);
				return;
			}
			add_soft_edge(x,v,pe-puv);
		}
	}
}
bool Graph::has_hard_edge(char* is,ui u,ui v){
	for(ui w:adj[u]){
		if(is[w]!=1) continue;
		if(w==v){
			return true;
		}
	}
	return false;
}
void Graph::has_hard(ui v,ui x,ui y,char* is,int puv,int& res_w){
	bool resx=false;
	bool resy=false;
	for(ui w:adj[v]){
		if(is[w]!=1) continue;
		if(w==x){
			resx=true;
		}else if(w==y){
			resy=true;
		}
		if(resx && resy) break;
	}
	bool hasx=false;
	bool hasy=false;
	if(resx &&resy) return;
	if(!resx){
		int pex=0;
		for(auto& w:generalized_adj[v]){
			if(w.first==x){
				pex=w.second;
				hasx=true;
				break;
			}
		}
		if(!hasx){
			soft_degree_increaseOne(is,x,res_w);
			soft_degree_increaseOne(is,v,res_w);
			add_soft_edge(x,v,-puv);
		}else{
			remove_soft_edge(x,v,pex);
			if(pex==puv){
			soft_degree_decreaseOne(is,x,res_w);
			soft_degree_decreaseOne(is,v,res_w);
				return ;
			}
			add_soft_edge(x,v,pex-puv);
		}
	}
	if(!resy){
		int pey=0;
		for(auto& w:generalized_adj[v]){
			if(w.first==y){
				pey=w.second;
				hasy=true;
				break;
			}
		}
		if(!hasy){
			soft_degree_increaseOne(is,y,res_w);
			soft_degree_increaseOne(is,v,res_w);
			add_soft_edge(y,v,-puv);
		}else{
			remove_soft_edge(y,v,pey);
			if(pey==puv){
			soft_degree_decreaseOne(is,y,res_w);
			soft_degree_decreaseOne(is,v,res_w);
				return;
			}
			add_soft_edge(y,v,pey-puv);
		}
	}
}
void Graph::reduction_hard_degree_one(char* is,int& res){
	for(int i=hard_degree_ones.vnum-1;i>=0;i--){
		ui u=hard_degree_ones.vlist[i];
		if(is[u]!=1 || hard_degreeInC[u]!=1) continue;
		if(weight_node[u]<=0) continue;
		int sum=0;
		for(auto& v:generalized_adj[u]){
			if(is[v.first]!=1) continue;
			sum+=max(0,v.second);
			if(sum>=weight_node[u]){
				break;
			}
		}
		if(sum>=weight_node[u]) continue;
		ui x=n+1;
		for(ui v:adj[u]){
			if(is[v]!=1) continue;
			x=v;
			break;
		}
		res+=weight_node[u];
		exclusive_pairs.emplace_back(Myturple(u,x,0));
		weight_node[x]-=weight_node[u];
		dynamic_weight[x]-=weight_node[u];
		hard_neighbors.clear();
		for(ui v:adj[x]){
			if(is[v]!=1) continue;
			hard_neighbors.add(v);
		}
		for(auto& v:generalized_adj[u]){
			if(is[v.first]!=1) continue;
			weight_node[v.first]-=v.second;
			dynamic_weight[v.first]-=v.second;
			if(!hard_neighbors.contains(v.first)){
				update_penalty(v.first,x,-v.second,is,res);
			}
		}
		delete4fold(u,res,is);
	}
	hard_degree_ones.clear();
	update_sets(is);
	 
}
void Graph::reduction_hard_degree_two(char* is,int& res){
	for(int i=hard_degree_twos.vnum-1;i>=0;i--){
		ui u=hard_degree_twos.vlist[i];
		if(is[u]!=1 || hard_degreeInC[u]!=2) continue;
		if(weight_node[u]<=0) continue;
		int sum=0;
		for(auto& v:generalized_adj[u]){
			if(is[v.first]!=1) continue;
			sum+=max(0,v.second);
			if(sum>weight_node[u]) break;
		}
		if(sum>weight_node[u]) continue;
		ui x=n+1,y=n+1;
		for(ui v:adj[u]){
			if(is[v]!=1) continue;
			if(x==n+1){
				x=v;
			}else{
				y=v;
				break;
			}
		}
			bool has_triangle=has_hard_edge(is,x,y);
			if(dynamic_weight[x]<dynamic_weight[y]){
				swap(x,y);
			}
			if(has_triangle){
				if(weight_node[u]>=dynamic_weight[x]){
					add_vertex(u,res,is);
					delete_vertex(x,res,is);
					delete_vertex(y,res,is);
				}
				else if(weight_node[u]>=dynamic_weight[y]){
					delete_vertex(y,res,is);
					delete4fold(u,res,is);
					res+=weight_node[u];
					exclusive_pairs.emplace_back(Myturple(u,x,0));
					weight_node[x]-=weight_node[u];
					dynamic_weight[x]-=weight_node[u];
					for(auto& v:generalized_adj[u]){
						if(is[v.first]!=1) continue;
						weight_node[v.first]-=v.second;
						dynamic_weight[v.first]-=v.second;
						has_hard(v.first,x,is,v.second,res) ;
					}
				}
				//case3:w(u)<w(y)
				else if(weight_node[u]<dynamic_weight[y]){
					delete4fold(u,res,is);
					res+=weight_node[u];
					exclusive_pairs.emplace_back(u,x,y,1);
					weight_node[x]-=weight_node[u];
					weight_node[y]-=weight_node[u];
					dynamic_weight[x]-=weight_node[u];
					dynamic_weight[y]-=weight_node[u];
					for(auto& v:generalized_adj[u]){
						if(is[v.first]!=1) continue;
						weight_node[v.first]-=v.second;
						dynamic_weight[v.first]-=v.second;
						has_hard(v.first,x,y,is,v.second,res);
					}
				}
			}
	}
	hard_degree_twos.clear();
	update_sets(is);
}
void Graph::update_swap_target(RandList& swap_set,vector<ui>& swap_target){
	for(int i=0;i<swap_set.vnum;i++){
		ui v=swap_set.vlist[i];
		if(swap_target[v]==-1 || !S.contains(swap_target[v])){
			for(ui u:adj[v]){
				if(S.contains(u)){
					swap_target[v]=u;
					break;
				}
			}
		}
	}
}
int Graph::select_drop(vector<int>& utility,int& gain,int res){
	int min_utility=INT32_MAX;
	int min_u=-1;
	vector<ui> min_list;
		if(false){
			int r=S.vlist[get_rand_n(S.vnum)];
			gain=-utility[r];
			return r;
		}else{
			for(int i=0;i<S.vnum;i++){
				ui u=S.vlist[i];
				int current_utility=utility[u];
				if(current_utility<min_utility){
					min_list.clear();
					min_list.emplace_back(u);
					min_utility=current_utility;
				}else if(current_utility==min_utility){
					min_list.emplace_back(u);
				}
			}
			if(min_list.empty()) return -1;
			int r=get_rand_n(min_list.size());
			gain=-min_utility;
			double stand=(double)min_utility/((double)res/S.vnum);
			if(get_rand_n(10)<10*stand){
				r=S.vlist[get_rand_n(S.vnum)];
				gain=-utility[r];
				return r;
			}
			return min_list[r];
		}
}
int Graph::select_add(RandList& add_set,vector<ui>& is_tabu,vector<int>& utility,int& res,int best_res,int& gain){
		vector<ui> in_tabu;
		vector<ui> not_tabu;
		int max_in_tabu=INT32_MIN;
		int max_not_tabu=INT32_MIN;
		for(int i=0;i<add_set.vnum;i++){
			ui u=add_set.vlist[i];
			if(is_tabu[u]){
				if(utility[u]>max_in_tabu){
					in_tabu.clear();
					max_in_tabu=utility[u];
					in_tabu.emplace_back(u);
				}else if(utility[u]==max_in_tabu){
					in_tabu.emplace_back(u);
				}
			}else{
				if(utility[u]>max_not_tabu){
					not_tabu.clear();
					max_not_tabu=utility[u];
					not_tabu.emplace_back(u);
				}else if(utility[u]==max_not_tabu){
					not_tabu.emplace_back(u);
				}
			}
		}
		if(!in_tabu.empty() && max_in_tabu>max_not_tabu && ((max_in_tabu + res) > best_res) && (max_not_tabu + res)<=best_res){
			gain=max_in_tabu;
			return in_tabu[get_rand_n(in_tabu.size())];
		}else if(!not_tabu.empty()){
			gain=max_not_tabu;
			return not_tabu[get_rand_n(not_tabu.size())];
		}else{
			return -1;
		}
}
int Graph::select_swap(RandList& swap_set,vector<ui>& swap_target,vector<ui>& is_tabu,vector<int>& utility,int& res,int best_res,int& gain){
		vector<ui> in_tabu;
		vector<ui> not_tabu;
		int max_in_tabu=INT32_MIN;
		int max_not_tabu=INT32_MIN;
		update_swap_target(swap_set,swap_target);
		for(int i=0;i<swap_set.vnum;i++){
			ui v=swap_set.vlist[i];
			ui u=swap_target[v];
			int current_utility=get_swap_utility(u,v,utility);
			if(is_tabu[v]){
				if(current_utility>max_in_tabu){
					in_tabu.clear();
					max_in_tabu=current_utility;
					in_tabu.emplace_back(v);
				}else if(current_utility==max_in_tabu){
					in_tabu.emplace_back(v);
				}
			}else{
				if(current_utility>max_not_tabu){
					not_tabu.clear();
					max_not_tabu=current_utility;
					not_tabu.emplace_back(v);
				}else if(current_utility==max_not_tabu){
					not_tabu.emplace_back(v);
				}
			}
		}
		if(!in_tabu.empty() && max_in_tabu>max_not_tabu && max_in_tabu + res > best_res && max_not_tabu + res <= best_res){
			gain=max_in_tabu;
			int r=get_rand_n(in_tabu.size());
			return in_tabu[r];
		}else if(!not_tabu.empty()){
			gain=max_not_tabu;
			int r=get_rand_n(not_tabu.size());
			return not_tabu[r];
		}else{
			return -1;
		}
}
int Graph::origin_search(int& res,vector<ui>& k2c,int depth,std::chrono::high_resolution_clock::time_point start_time,int run_time,double & best_time,int& global_best,bool need2reduction){
	RandList add_set;
	priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>> tabu;
	add_set.init(n);
	RandList swap_set;
	swap_set.init(n);
	int current;
	vector<ui> swap_target(n,n+1);
	vector<int> utility_add;
	vector<ui> is_tabu;
	vector<ui> hard_degree_in_S;
	hard_degree_in_S.resize(n,0);
	is_tabu.resize(n,0);
	utility_add.resize(n);
	char* is=new char[n];
	for(int i=0;i<n;i++) is[i]=1;
	for(int i=0;i<n;i++){
		if(!S.contains(i)){
			add_set.add(i);
			is[i]=1;
		}
		utility_add[i]=get_add_utility(i);
	}
	for(int i=0;i<S.vnum;i++){
		ui u=S.vlist[i];
		is[u]=2;
		for(ui v:adj[u]){
			hard_degree_in_S[v]++;
			if(hard_degree_in_S[v]==1){
				add_set.remove(v);
				swap_set.add(v);
				swap_target[v]=u;
			}else if(hard_degree_in_S[v]==2){
				swap_set.remove(v);
				swap_target[v]=-1;
			}
		}
	}
	int best_res=res;
	current=0;
	int l=0;
	RandList best_S;
	best_S.init(n);
	best_S=S;
	while(add_set.vnum){
		ui u=add_set.vlist[get_rand_n(add_set.vnum)];
		int gain=0;
		Add_Set(u,res,current,tabu,is_tabu,add_set,swap_set,swap_target,utility_add,hard_degree_in_S,is,best_res);
		l++;
		current++;
	}
	bool can_drop=false;
	bool first_swap=true;
	int not_increase=0;
	if(res>best_res){
		best_res=res;
		best_S=S;
		if(res>global_best){
			auto current_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_time = current_time - start_time;
			best_time=elapsed_time.count();
		}
	}
	int not_work=0;
	while(l<depth || !add_set.empty()){
		auto current_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_time = current_time - start_time;
		if (elapsed_time.count() >= run_time) {
			break; 
		}
	current++;
	if(!tabu.empty())
	while(!tabu.empty() && tabu.top().first<current){
		ui min=tabu.top().second;
		is_tabu[min]=0;
		tabu.pop();
	}
	if(l>=depth){
		int old_current=current;
		current+=10*(S.vnum+add_set.vnum);
		int num=0;
		while(add_set.vnum){
			num++;
			int gain=0;
			int u=select_add(add_set,is_tabu,utility_add,res,best_res,gain);
			if(gain<=0 || u==-1) {
				current=old_current;
				break;
				}
			Add_Set(u,res,current,tabu,is_tabu,add_set,swap_set,swap_target,utility_add,hard_degree_in_S,is,best_res);
			if(res>best_res){
				best_res=res;
				best_S=S;
		if(res>global_best){
			auto current_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_time = current_time - start_time;
			best_time=elapsed_time.count();
		}
			}
		}
		current=old_current;
		break;
	}
	if(not_increase>S.vnum/2){
		int delete_num=0;
			delete_num=S.vnum/5;
			delete_num=max(delete_num,1);
		for(int i=0;i<delete_num;i++){
			int gain;
			int u=select_drop(utility_add,gain,res);
			Drop_Set(u,res,current,tabu,is_tabu,add_set,swap_set,swap_target,utility_add,hard_degree_in_S,can_drop,is);
		}
		not_increase=0;
		not_work++;
		continue;
	}
	int add_pos=-1,swap_pos=-1,drop_pos=-1;
	int gain_add=INT32_MIN,gain_swap=INT32_MIN,gain_drop=INT32_MIN;
		l++;
		not_increase++;
		add_pos=select_add(add_set,is_tabu,utility_add,res,best_res,gain_add);
		swap_pos=select_swap(swap_set,swap_target,is_tabu,utility_add,res,best_res,gain_swap);
			if(add_pos!=-1){
				if(swap_pos!=-1){
					if(gain_add>gain_swap){
						Add_Set(add_pos,res,current,tabu,is_tabu,add_set,swap_set,swap_target,utility_add,hard_degree_in_S,is,best_res);
					}else{
						Swap_Set(swap_pos,res,current,tabu,is_tabu,swap_set,swap_target,add_set,utility_add,hard_degree_in_S,is,best_res);
					}	
				}else{
					Add_Set(add_pos,res,current,tabu,is_tabu,add_set,swap_set,swap_target,utility_add,hard_degree_in_S,is,best_res);
				}
			}else{
				if(swap_pos!=-1){
					Swap_Set(swap_pos,res,current,tabu,is_tabu,swap_set,swap_target,add_set,utility_add,hard_degree_in_S,is,best_res);
				}else{
					int delete_num=S.vnum/5;
					delete_num=max(delete_num,1);
					for(int i=0;i<delete_num;i++){
						int gain;
						int u=select_drop(utility_add,gain,res);
						Drop_Set(u,res,current,tabu,is_tabu,add_set,swap_set,swap_target,utility_add,hard_degree_in_S,can_drop,is);
					}
					not_increase=0;
					not_work++;
					continue;
					;
					}
			}
			if(res>best_res){
		if(res>global_best){
			auto current_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_time = current_time - start_time;
			best_time=elapsed_time.count();
		}
				best_res=res;
				best_S=S;
				not_increase=0;
				not_work=0;
			}
		int old_num;
			if(res>best_res){
		if(res>global_best){
			auto current_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_time = current_time - start_time;
			best_time=elapsed_time.count();
		}
				best_res=res;
				best_S=S;
				not_increase=0;
				not_work=0;
			}
	}
	S=best_S;
	res=best_res;
	return current;
}
void Graph::unfold(){
		RandList fold_set;
		fold_set.init(n);
		for(auto& myturple:exclusive_pairs){
			switch (myturple.state)
			{
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
				fold_set.add(myturple.u);
				break;
			default:
			break;
			}
		}
		int Sturn=0;
		vector<int> processed(n,0);
		bool flag=true;
		//
		int num=0;
		while(flag){
			num++;
			flag=false;
			int case0=0,case1=0,case2=0,case3=0,case4=0,case5=0,case6=0,case7=0;
			for(int i=exclusive_pairs.size()-1;i>=0;i--){
				if(processed[i]==1) continue; 
				auto& myturple=exclusive_pairs[i];
				switch (myturple.state)
				{
				case 0:
					if(fold_set.contains(myturple.x)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 &&!S.contains(myturple.x)){
						S.add(myturple.u);
					}else if(myturple.u==n+1 && !S.contains(myturple.x) && myturple.x!=n+1){
						S.add(myturple.x);
					}
					case0++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 1:
					if(fold_set.contains(myturple.x) || fold_set.contains(myturple.y)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && (!S.contains(myturple.x) && !S.contains(myturple.y))){
						S.add(myturple.u);
					}else if(myturple.u==n+1 && !S.contains(myturple.x) && !S.contains(myturple.y) && myturple.x==n+1 && myturple.y!=n+1){
						S.add(myturple.y);
					}else if(myturple.u==n+1 && !S.contains(myturple.x) && !S.contains(myturple.y) && myturple.x!=n+1 && myturple.y==n+1){
						S.add(myturple.x);
					}else if(myturple.u==n+1 && !S.contains(myturple.x) && !S.contains(myturple.y) && myturple.x!=n+1 && myturple.y!=n+1){
						S.add(myturple.x);
					}
					case1++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 2:
					if(fold_set.contains(myturple.x) || fold_set.contains(myturple.y)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && !(S.contains(myturple.x) && S.contains(myturple.y))){
						S.add(myturple.u);
					}else if(myturple.u==n+1 && myturple.x!=n+1 && myturple.y!=n+1){
						S.add(myturple.x);
						S.add(myturple.y);
					}
					case2++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 3:
					if(fold_set.contains(myturple.x) || fold_set.contains(myturple.y)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && (!S.contains(myturple.x) || S.contains(myturple.y))){
						S.add(myturple.u);
					}else if(myturple.u==n+1 && myturple.y!=n+1){
						S.add(myturple.y);
					}
					case3++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 4:
					if(fold_set.contains(myturple.x)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && S.contains(myturple.x)){
						S.add(myturple.u);
					}
					case4++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 5:
					if(fold_set.contains(myturple.x) || fold_set.contains(myturple.y)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && (!S.contains(myturple.x) && S.contains(myturple.y))){
						S.add(myturple.u);
					}else if(myturple.u==n+1 && !S.contains(myturple.y) && myturple.x!=n+1){
						S.add(myturple.x);
					}
					case5++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 6:
					if(fold_set.contains(myturple.x) || fold_set.contains(myturple.y)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && (S.contains(myturple.x) || S.contains(myturple.y))){
						S.add(myturple.u);
					}
					case6++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				case 7:
					if(fold_set.contains(myturple.x) || fold_set.contains(myturple.y)){
						flag=true;
						break;
					}
					if(myturple.u!=n+1 && (S.contains(myturple.x) && S.contains(myturple.y))){
						S.add(myturple.u);
					}
					case7++;
					processed[i]=1;
					fold_set.remove(myturple.u);
					break;
				default: cout<<"state:"<<myturple.state<<endl;
				}
			}
		}
}
void Graph::twin(char* is,int &res,bool init){
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	if(Cand_reduce.size()==1) {
		Cand_reduce.clear();
		return;
	}
	vector<ui> add_list;
	temp_cand.clear();
	for(int u=0;u<Cand_reduce.size()-1;u++){
		ui i=Cand_reduce[u];
		if(is[i]!=1) continue;
		hard_neighbors.clear();
		for(ui v:adj[i]){
			if(is[v]!=1) continue;
			hard_neighbors.add(v);
		}
		int sum=0;
		vector<int> neighbor;
		neighbor.resize(n);
		for(auto& v:generalized_adj[i]){
			if(is[v.first]!=1) continue;
			sum+=max(0,v.second);
			neighbor[v.first]=v.second;
			if(sum>weight_node[i]){
				break;
			}
		}
		if(sum>weight_node[i]) continue;
		if(hard_degreeInC[i]==0) {
			add_list.emplace_back(i);
			continue;
			}
		for(int s=u+1;s<Cand_reduce.size();s++){
			ui j=Cand_reduce[s];
			if(is[j]!=1 || hard_neighbors.contains(j) || j<=i || hard_degreeInC[j]==0 || hard_degreeInC[j]!=hard_degreeInC[i]) continue;
			if(hard_degreeInC[i]==hard_degreeInC[j]){
				bool is_equal=true;
				for(ui v:adj[j]){
					if(is[v]!=1) continue;
					if(!hard_neighbors.contains(v)){
						is_equal=false;
						break;
					}
				}
				if(!is_equal){
					break;
				}
				sum=0;
				int puv=0;
				for(auto& v:generalized_adj[j]){
					if(is[v.first]!=1) continue;
					sum+=max(0,v.second);
					if(v.first==i) {puv=v.second;}
					if(sum>weight_node[j]){
						break;
					}
				}	
				if(sum>weight_node[j]) break;
				for(auto& v:generalized_adj[j]){
					if(is[v.first]!=1 || v.first==i) continue;
					if(neighbor[v.first]!=0){
						remove_soft_edge(v.first,i,neighbor[v.first]);
						add_soft_edge(v.first,i,v.second+neighbor[v.first]);
					}else{
						add_soft_edge(v.first,i,v.second);
						soft_degree_increaseOne(is,v.first,res);
						soft_degree_increaseOne(is,i,res);
					}
				}
				weight_node[i]+=weight_node[j]-puv;
				dynamic_weight[i]+=weight_node[j]-puv;
				exclusive_pairs.emplace_back(Myturple(j,i,4));
				delete4fold(j,res,is);
				for(ui w:adj[j]){
					if(is[w]==1)
					temp_cand.emplace_back(w);
				}
				for(auto& w:generalized_adj[j]){
					if(is[w.first]==1 && w.second<0){
						temp_cand.emplace_back(w.first);
					}
				}
			}
		}
	}
	for(ui u:add_list){
		add_vertex(u,res,is);
	}
	Cand_reduce.clear();
	if(!temp_cand.empty())
	Cand_reduce=temp_cand;
}

void Graph::basic_single_edge(char* is,int& res,bool init){
	vector<ui> need2remove;
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n-1;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui u:Cand_reduce){
		if(is[u]!=1 || S.contains(u)) continue;
		if(weight_node[u]<0) continue;
		need2remove.clear();
		need2remove.reserve(hard_degreeInC[u]);
		int old_weight_hard=0;
		hard_neighbors.clear();
		for(ui v:adj[u]){
			if(is[v]!=1) continue;
			old_weight_hard+=max(0,dynamic_weight[v]);
			hard_neighbors.add(v);
		}
		int old_weight_soft=0;
		int old_weight_edge=0;
		soft_neighbors.clear();
		vector<int> pux;
		pux.resize(n,0);
		for(auto& w:generalized_adj[u]){
			if(is[w.first]!=1) continue;
			soft_neighbors.add(w.first);
			pux[w.first]=max(0,w.second);
			old_weight_soft+=max(0,dynamic_weight[w.first]);
			old_weight_edge+=max(0,w.second);
		}
		for(ui v:adj[u]){
			if(is[v]!=1 || v<=u) continue;
			if(weight_node[u]<dynamic_weight[v]) continue;
			int weight_soft=old_weight_soft;
			int weight_edge=old_weight_edge;
			int weight_hard=old_weight_hard;
			for(ui w:adj[v]){
				if(is[w]!=1) continue;
				if(soft_neighbors.contains(w)){
					weight_soft-=max(0,dynamic_weight[w]);
					weight_edge-=max(0,pux[w]);
				}
				if(hard_neighbors.contains(w)){
					weight_hard-=max(0,dynamic_weight[w]);
				}
			}
			int min_=min(weight_soft,weight_edge);
			if(weight_node[u]>=weight_hard+min_){
				need2remove.emplace_back(v);
				is[v]=0;
			}
		}
		for(ui v:need2remove){
			is[v]=1;
			delete_vertex(v,res,is);
			for(ui w:adj[v]){
				if(is[w]==1)
				temp_cand.emplace_back(w);
			}
			for(auto& w:generalized_adj[v]){
				if(is[w.first]==1)
				temp_cand.emplace_back(w.first);
			}
		}
	}
	Cand_reduce.clear();
	Cand_reduce=temp_cand;
}
void Graph::extend_single_edge(char* is,int& res,bool init){
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui u:Cand_reduce){
		if(is[u]!=1) continue;
		if(weight_node[u]<0 || S.contains(u)) continue;
		vector<ui> need2remove;
		need2remove.reserve(hard_degreeInC[u]);
		hard_neighbors.clear();
		int sum=0;
		for(ui v:adj[u]){
			if(is[v]!=1) continue;
			sum+=max(0,dynamic_weight[v]);
			hard_neighbors.add(v);
		}
		for(auto& v:generalized_adj[u]){
			if(is[v.first]!=1) continue;
			sum+=max(0,dynamic_weight[v.first]);
		}
		for(ui v:adj[u]){
			if(weight_node[v]<0 || S.contains(v)) continue;
			if(is[v]!=1 || v<u) continue;
			if(weight_node[u]>=sum-weight_node[v]){
				for(ui w:adj[v]){
					if(is[w]!=1) continue;
					if(hard_neighbors.contains(w)){
						need2remove.emplace_back(w);
						is[w]=0;
					}
				}
			}
		}
		for(ui v:need2remove){
			is[v]=1;
			delete_vertex(v,res,is);
			for(ui w:adj[v]){
				if(is[w]==1)
				temp_cand.emplace_back(w);
			}
			for(auto& w:generalized_adj[v]){
				if(is[w.first]==1)
				temp_cand.emplace_back(w.first);
			}
		}
	}
	Cand_reduce.clear();
	Cand_reduce=temp_cand; 
}
void Graph::edge_domination(char* is, int &res,bool init){
	vector<ui> need2remove;
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui u:Cand_reduce){
		if(is[u]!=1) continue;
		hard_neighbors.clear();
		for(ui v:adj[u]){
			if(is[v]!=1) continue;
			hard_neighbors.add(v);
		}
		unordered_map<ui,int> soft;
		for(auto& v:generalized_adj[u]){
			if(is[v.first]!=1) continue;
			soft.insert(v);
		}
		bool u_d_v=true;
		for(ui v:adj[u]){
			if(is[v]!=1 || v<u || hard_degreeInC[v]>hard_degreeInC[u] ||weight_node[v]<weight_node[u]) continue;
			for(ui w:adj[v]){
				if(is[w]!=1) continue;
				if(!hard_neighbors.contains(w)){
					u_d_v=false;
					break;
				}
			}
			if(!u_d_v) continue;
			bool sat=true;
			for(auto& w:generalized_adj[v]){
				if(is[w.first]!=1) continue;
				if(soft.find(w.first)!=soft.end()){
					if(w.second>soft[w.first]){
						sat=false;
						break;
					}
				}else if(soft.find(w.first)==soft.end()){
					if(w.second>0){
						sat=false;
						break;
					}
				}
			}
			if(sat){
				delete_vertex(u,res,is);
				for(ui w:adj[u]){
					if(is[w]==1)
					temp_cand.emplace_back(w);
				}
				for(auto& w:generalized_adj[u]){
					if(is[w.first]==1)
					temp_cand.emplace_back(w.first);
				}
				break;
			}
		}
	}
	 Cand_reduce.clear();
	 Cand_reduce=temp_cand;

}
void Graph::reduction_degree_one(char* is,int& res){
	for(int i=total_degree_ones.vnum-1;i>=0;i--){
		ui u=total_degree_ones.vlist[i];
		if(is[u]!=1 || hard_degreeInC[u]+soft_degreeInC[u]!=1) continue;
			pair<ui,int> v;
			v.first=n+1;
			for(auto& w:generalized_adj[u]){
				if(is[w.first]==1 && v.first==n+1) {
					v=w;
					break;
					}
			}
			if(v.first==n+1)
			for(auto& w:adj[u]){
				if(is[w]==1 && v.first==n+1) {
					v.first=w;
					v.second=INT32_MAX;
					break;
					}
			}
			if(weight_node[u]>=v.second){
				if(weight_node[u]>=0){
					add_vertex(u,res,is);
				}else{
					delete4fold(u,res,is);
					weight_node[v.first]=weight_node[v.first]-v.second+weight_node[u];
					dynamic_weight[v.first]=dynamic_weight[v.first]-v.second+weight_node[u];
					exclusive_pairs.emplace_back(Myturple(u,v.first,4));
				}
			}else if(weight_node[u]>=dynamic_weight[v.first]){
				if(weight_node[u]>=0){
					add_vertex(u,res,is);
					delete_vertex(v.first,res,is);
				}else{
					delete_vertex(u,res,is);
				}
			}else if(weight_node[u]<dynamic_weight[v.first]){
				if(weight_node[u]>=0){
					delete4fold(u,res,is);
					res+=weight_node[u];
					exclusive_pairs.emplace_back(Myturple(u,v.first,0));
					weight_node[v.first]-=weight_node[u];
					dynamic_weight[v.first]-=weight_node[u];
				}else{
					delete_vertex(u,res,is);
				}
			}
		// }
	}
	total_degree_ones.clear();
	update_sets(is);
} 
void Graph::reduction_triangle(char* is,int& res){
	for(int i=total_degree_twos.vnum-1;i>=0;i--){
		ui u=total_degree_twos.vlist[i];
		if(is[u]!=1 || hard_degreeInC[u]+soft_degreeInC[u]!=2) continue;
			pair<ui,int> x,y;
			x.first=n+1;
			y.first=n+1;
			for(auto& v:generalized_adj[u]){
				if(is[v.first]!=1) continue;
				if(x.first==n+1){
					x=v;
				}else{
					y=v;
				}
			}
			for(ui v:adj[u]){
				if(is[v]!=1) continue;
				if(x.first==n+1){
					x.first=v;
					x.second=INT32_MAX/2;
				}else{
					y.first=v;
					y.second=INT32_MAX/2;
				}
			}
			bool has_triangle=has_hard_edge(is,x.first,y.first);
			if(x.second<y.second){
				swap(x,y);
			}
			if(has_triangle) {
				if(weight_node[u]>=x.second){
					if(weight_node[u]>=0){
						add_vertex(u,res,is);
					}else{
						delete4fold(u,res,is);
						weight_node[x.first]+=weight_node[u]-x.second;
						dynamic_weight[x.first]+=weight_node[u]-x.second;
						weight_node[y.first]+=weight_node[u]-y.second;
						dynamic_weight[y.first]+=weight_node[u]-y.second;
						exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,6));
					}
				}
				else if(weight_node[u]>=y.second && weight_node[u]<x.second){
					if(weight_node[u]>=0){
						if(weight_node[u]>=dynamic_weight[x.first]){
							add_vertex(u,res,is);
							delete_vertex(x.first,res,is);
						}else{
							delete4fold(u,res,is);
							exclusive_pairs.emplace_back(Myturple(u,x.first,0));
							res+=weight_node[u];
							weight_node[x.first]-=weight_node[u];
							dynamic_weight[x.first]-=weight_node[u];
							weight_node[y.first]-=y.second;
							dynamic_weight[y.first]-=y.second;
						}
					}else{
						delete4fold(u,res,is);
						weight_node[y.first]+=(weight_node[u]-y.second);
						dynamic_weight[y.first]+=(weight_node[u]-y.second);
						exclusive_pairs.emplace_back(Myturple(u,y.first,4));
					}
				}
				else if(weight_node[u]<x.second && weight_node[u]<y.second){
					if(dynamic_weight[x.first]<dynamic_weight[y.first]){
						swap(x,y);
					}
					if(weight_node[u]>=0){
						if(weight_node[u]>=dynamic_weight[x.first]){
							add_vertex(u,res,is);
							delete_vertex(x.first,res,is);
							delete_vertex(y.first,res,is);
						}else if(weight_node[u]>=dynamic_weight[y.first]){
							delete_vertex(y.first,res,is);
							delete4fold(u,res,is);
							res+=weight_node[u];
							exclusive_pairs.emplace_back(Myturple(u,x.first,0));
							weight_node[x.first]-=weight_node[u];
							dynamic_weight[x.first]-=weight_node[u];
						}else{
							delete4fold(u,res,is);
							res+=weight_node[u];
							exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,1));
							weight_node[x.first]-=weight_node[u];
							weight_node[y.first]-=weight_node[u];
							dynamic_weight[x.first]-=weight_node[u];
							dynamic_weight[y.first]-=weight_node[u];
						}
					}else{
						delete_vertex(u,res,is);
					}
				}	
			}else{
				if(weight_node[u]<y.second){
					if(weight_node[u]>=0){
							res+=weight_node[u];
							exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,1));
							weight_node[x.first]-=weight_node[u];
							weight_node[y.first]-=weight_node[u];
							dynamic_weight[x.first]-=weight_node[u];
							dynamic_weight[y.first]-=weight_node[u];
							update_penalty(x.first,y.first,-weight_node[u],is,res);
							delete4fold(u,res,is);
						}else{ 
							if(weight_node[u]<y.second+x.second){
								delete_vertex(u,res,is);
							}else{						
								exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,7));
								update_penalty(x.first,y.first,x.second+y.second-weight_node[u],is,res);
								delete4fold(u,res,is);
							}
						
						}
				}else if(weight_node[u]<x.second){
					if(weight_node[u]>=0){
						if(weight_node[u]>=x.second+y.second){
							res+=weight_node[u];
							exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,3));
							weight_node[y.first]-=y.second;
							dynamic_weight[y.first]-=y.second;
							weight_node[x.first]-=weight_node[u];
							dynamic_weight[x.first]-=weight_node[u];
							update_penalty(x.first,y.first,x.second-weight_node[u],is,res);
							delete4fold(u,res,is);
						}else{
							res+=weight_node[u];
							exclusive_pairs.emplace_back(Myturple(u,x.first,0));
							weight_node[y.first]-=y.second;
							dynamic_weight[y.first]-=y.second;
							weight_node[x.first]-=weight_node[u];
							dynamic_weight[x.first]-=weight_node[u];
							update_penalty(x.first,y.first,-y.second,is,res);
							delete4fold(u,res,is);
						}
					}else{
						if(weight_node[u]>=x.second+y.second){
							exclusive_pairs.emplace_back(Myturple(u,y.first,4));
							weight_node[y.first]+=(weight_node[u] - y.second);
							dynamic_weight[y.first]+=(weight_node[u] - y.second);
							update_penalty(x.first,y.first,x.second,is,res);
							delete4fold(u,res,is);
						}else{
							exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,5));
							weight_node[y.first]+=(weight_node[u] - y.second);
							dynamic_weight[y.first]+=(weight_node[u] - y.second);
							update_penalty(x.first,y.first,weight_node[u]-y.second,is,res);
							delete4fold(u,res,is);
						}
					}
			}else{
				if(weight_node[u]>=0){
						if(weight_node[u]>=x.second+y.second){
							add_vertex(u,res,is);
						}else{
							res+=weight_node[u];
							exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,2));
							weight_node[y.first]-=y.second;
							dynamic_weight[y.first]-=y.second;
							weight_node[x.first]-=x.second;
							dynamic_weight[x.first]-=x.second;
							update_penalty(x.first,y.first,weight_node[u]-x.second-y.second,is,res);
							delete4fold(u,res,is);
						}
				}else{
					exclusive_pairs.emplace_back(Myturple(u,x.first,y.first,6));
					weight_node[y.first]+=(weight_node[u]-y.second);
					dynamic_weight[y.first]+=(weight_node[u]-y.second);
					weight_node[x.first]+=(weight_node[u]-x.second);
					dynamic_weight[x.first]+=(weight_node[u]-x.second);
					if(weight_node[u]!=0)
					update_penalty(x.first,y.first,weight_node[u],is,res);
					delete4fold(u,res,is);
				}
			}
			}
			
	}
	total_degree_twos.clear();
	update_sets(is);
	 
}
void Graph::isolated(char* is,int& res,bool init){
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui i:Cand_reduce){
		if(is[i]!=1) continue;
		if(weight_node[i]<=0) continue;
		bool is_clique=true;
		int sum=0;
		if(soft_degreeInC[i]>0)
		for(auto& v:generalized_adj[i]){
			if(is[v.first]!=1) continue;
			sum+=max(0,v.second);
			if(sum>weight_node[i]){
				is_clique=false;
				break;
			}
		}
		if(!is_clique) continue;
		for(ui v:adj[i]){
			if(is[v]!=1) continue;
			if(weight_node[i]<sum+dynamic_weight[v]){
				is_clique=false;
				break;
			}
		}
		if(!is_clique) continue;
		for(ui v:adj[i]){
			if(is[v]!=1) continue;
			if(hard_degreeInC[v]<hard_degreeInC[i]){
				is_clique=false;
				break;
			}
			hard_neighbors.clear();
			for(ui w:adj[v]){
				if(is[w]!=1) continue;
				hard_neighbors.add(w);
			}
			for(ui w:adj[i]){
				if(is[w]!=1 || w==v) continue;
				if(!hard_neighbors.contains(w)){
					is_clique=false;
					break;
				}
			}
			if(!is_clique) break;
		}
		if(is_clique){
			add_vertex(i,res,is);
			for(auto& w:generalized_adj[i]){
				if(is[w.first]==1)
				temp_cand.emplace_back(w.first);
			}
		}
	}
	Cand_reduce.clear();
	Cand_reduce=temp_cand;
}
Graph* Graph::copy_Cand(char* is,int& res,vector<ui>& c2k,vector<ui>& k2c){
	Graph* kernel=new Graph();
	kernel->dir=dir;
	int new_n=0;
	kernel->m=0;
	kernel->generalized_m=0;
	k2c.resize(n,n+1);
	for(int i=0;i<n;i++){
		if(is[i]!=1) continue;
		c2k.emplace_back(i);
		k2c[i]=new_n++;
	}
	kernel->n=new_n;
	kernel->hard_degreeInC.resize(new_n,0);
	kernel->soft_degreeInC.resize(new_n,0);
	kernel->weight_node.resize(new_n,0);
	kernel->dynamic_weight.resize(new_n,0);
	kernel->adj.resize(new_n);
	kernel->generalized_adj.resize(new_n);
	for(int i=0;i<new_n;i++){
		for(ui v:adj[c2k[i]]){
			if((is[v]!=1)|| k2c[v]<=i) continue;
			kernel->add_edge(k2c[v],i);
			kernel->m++;
			kernel->hard_degreeInC[i]++;
			if(is[c2k[i]]==1)
			kernel->hard_degreeInC[k2c[v]]++;
		}
		for(auto& v:generalized_adj[c2k[i]]){
			if((is[v.first]!=1)||k2c[v.first]<=i) continue;
			kernel->add_soft_edge(k2c[v.first],i,v.second);
			kernel->generalized_m++;
			kernel->soft_degreeInC[i]++;
			if(is[c2k[i]]==1)
			kernel->soft_degreeInC[k2c[v.first]]++;
		}
		kernel->weight_node[i]=weight_node[c2k[i]];
		kernel->dynamic_weight[i]=dynamic_weight[c2k[i]];
	}
	kernel->soft_neighbors.init(new_n);
	kernel->hard_neighbors.init(new_n);
	kernel->hard_degree_ones.init(new_n);
	kernel->hard_degree_twos.init(new_n);
	kernel->total_degree_ones.init(new_n);
	kernel->total_degree_twos.init(new_n);
	kernel->temp_delete_hard_ones.init(new_n);
	kernel->temp_delete_hard_twos.init(new_n);
	kernel->temp_delete_total_ones.init(new_n);
	kernel->temp_delete_total_twos.init(new_n);
	kernel->S.init(new_n);
	return kernel;
}
Graph* Graph::generate_kernel(char* is,int& res,vector<ui>& k2e,vector<ui>& e2k){
	Graph* kernel=new Graph();
	kernel->dir=dir;
	int new_n=0;
	for(int i=0;i<n;i++){
		if(is[i]==0) continue;
		k2e.emplace_back(i);
		e2k[i]=new_n++;
	}
	kernel->n=new_n;
	kernel->hard_degreeInC.resize(new_n,0);
	kernel->soft_degreeInC.resize(new_n,0);
	kernel->weight_node.resize(new_n,0);
	kernel->dynamic_weight.resize(new_n,0);
	kernel->adj.resize(new_n);
	kernel->generalized_adj.resize(new_n);
	kernel->exclusive_pairs.reserve(exclusive_pairs.size()+1);
	for(auto& myturple:exclusive_pairs){
		if(myturple.state==-1) continue;
		if(myturple.u!=-1) myturple.u=e2k[myturple.u];
		if(myturple.x!=-1) myturple.x=e2k[myturple.x];
		if(myturple.y!=-1) myturple.y=e2k[myturple.y];
		if(myturple.fold_vertex!=-1) myturple.fold_vertex=e2k[myturple.fold_vertex];
		kernel->exclusive_pairs.emplace_back(myturple);
	}
	for(int i=0;i<new_n;i++){
		for(ui v:adj[k2e[i]]){
			if((is[v]==0)|| e2k[v]<=i) continue;
			kernel->add_edge(e2k[v],i);
			if(is[v]==1)
			kernel->hard_degreeInC[i]++;
			if(is[k2e[i]]==1)
			kernel->hard_degreeInC[e2k[v]]++;
		}
		for(auto& v:generalized_adj[k2e[i]]){
			if((is[v.first]==0)||e2k[v.first]<=i) continue;
			kernel->add_soft_edge(e2k[v.first],i,v.second);
			if(is[v.first]==1)
			kernel->soft_degreeInC[i]++;
			if(is[k2e[i]]==1)
			kernel->soft_degreeInC[e2k[v.first]]++;
		}
		kernel->weight_node[i]=weight_node[k2e[i]];
		kernel->dynamic_weight[i]=dynamic_weight[k2e[i]];
	}
	kernel->soft_neighbors.init(new_n);
	kernel->hard_neighbors.init(new_n);
	kernel->hard_degree_ones.init(new_n);
	kernel->hard_degree_twos.init(new_n);
	kernel->total_degree_ones.init(new_n);
	kernel->total_degree_twos.init(new_n);
	kernel->temp_delete_hard_ones.init(new_n);
	kernel->temp_delete_hard_twos.init(new_n);
	kernel->temp_delete_total_ones.init(new_n);
	kernel->temp_delete_total_twos.init(new_n);
	kernel->S.init(new_n);
	return kernel;
}
void Graph::reduction_neighborhood4search(int& res,RandList& add_set,RandList& swap_set,vector<ui>& swap_target,vector<int>& utility,vector<ui>& hard_degree_in_S,char* is){
		RandList remove_from_add;
		remove_from_add.init(n);
		for(int j=0;j<add_set.vnum;++j){
			ui i=add_set.vlist[j];
			if(weight_node[i]<0 || remove_from_add.contains(i)) continue;
			int sum=0;
			for(auto& v :adj[i]){
					sum+=max(0,dynamic_weight[v]);
					if(sum>weight_node[i]) break ;
			}
			if(sum>weight_node[i]) continue;
			for(auto& v :generalized_adj[i]){
					sum+=max(0,dynamic_weight[v.first]);
					if(sum>weight_node[i]) break ;
			}
			if(sum>weight_node[i]) continue;
			if(sum<=weight_node[i]){
				int max_vertex=i;
				S.add(max_vertex);
				remove_from_add.add(max_vertex);
				for(ui v:adj[max_vertex]){
					hard_degree_in_S[v]++;
					if(hard_degree_in_S[v]==1){
						remove_from_add.add(v);
						swap_set.add(v);
						swap_target[v]=max_vertex;
					}else if(hard_degree_in_S[v]==2){
						swap_set.remove(v);
						swap_target[v]=-1;
					}
				}
				res+=weight_node[max_vertex];
				for(auto& v:generalized_adj[max_vertex]){
					utility[v.first]-=v.second;
					if(S.contains(v.first)){
						res-=v.second;
					}
				}
			}
	}
	for(int i=0;i<remove_from_add.vnum;i++){
		add_set.remove(remove_from_add.vlist[i]);
	}
}
void Graph::reduction_penalty4search(int& res,RandList& add_set,RandList& swap_set,vector<ui>& swap_target,vector<int>& utility,vector<ui>& hard_degree_in_S,char* is){
		RandList remove_from_add;
		remove_from_add.init(n);
		for(int j=0;j<add_set.vnum;++j){
			ui i=add_set.vlist[j];
			if(weight_node[i]<0 || remove_from_add.contains(i)) continue;
			int sum=0;
				for(ui v:adj[i]){
					sum+=max(0,dynamic_weight[v]);
					if(sum>weight_node[i]) {
						break;
					}
				}
			if(sum>weight_node[i]) continue;
				for(auto& v:generalized_adj[i]){
					sum+=max(0,v.second);
					if(sum>weight_node[i]) {
						break;
					}
				}
			if(sum>weight_node[i]) continue;
			if(sum<=weight_node[i]){
				int max_vertex=i;
				S.add(max_vertex);
				remove_from_add.add(max_vertex);
				for(ui v:adj[max_vertex]){
					hard_degree_in_S[v]++;
					if(hard_degree_in_S[v]==1){
						remove_from_add.add(v);
						swap_set.add(v);
						swap_target[v]=max_vertex;
					}else if(hard_degree_in_S[v]==2){
						swap_set.remove(v);
						swap_target[v]=-1;
					}
				}
				res+=weight_node[max_vertex];
				for(auto& v:generalized_adj[max_vertex]){
					utility[v.first]-=v.second;
					if(S.contains(v.first)){
						res-=v.second;
					}
				}
			}
	}
	for(int i=0;i<remove_from_add.vnum;i++){
		add_set.remove(remove_from_add.vlist[i]);
	}
}
void Graph::reduction_penalty(char* is,int& res,bool init){
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui i:Cand_reduce){
			if(is[i]!=1) continue;
			int sum=0;
				for(ui v:adj[i]){
					if(is[v]!=1) continue;
					sum+=max(0,dynamic_weight[v]);
					if(sum>weight_node[i]) {
						break;
					}
				}
			if(sum>weight_node[i]) continue;
				for(auto& v:generalized_adj[i]){
					if(is[v.first]!=1) continue;
					sum+=max(0,v.second);
					if(sum>weight_node[i]) {
						break;
					}
				}
			if(sum>weight_node[i]) continue;
			if(sum<=weight_node[i]){
				add_vertex(i,res,is);
				for(auto& w:generalized_adj[i]){
					if(is[w.first]==1)
					temp_cand.emplace_back(w.first);
				}
			}
	}
	Cand_reduce.clear();
	Cand_reduce=temp_cand;
}

int Graph::get_rand_n(int valid_n){
	return rand()%valid_n;
}
void Graph::write_solution(string file_name){
    std::ofstream outFile(file_name);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
	for(int i=0;i<S.vnum;i++){
		ui u=S.vlist[i];
		outFile << u+1 << std::endl;
	}
    outFile.close();
	cout<<"writing successfully"<<endl;
}
void Graph::write_txt(string file_name){
    std::ofstream outFile(file_name);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
	m=0;
	generalized_m=0;
	for(int i=0;i<n;i++){
		for(ui v:adj[i]){
			if(v>i) m++;
		}
		for(auto& v:generalized_adj[i]){
			if(v.first>i) generalized_m++;
		}
	}
 	outFile << "p " << "edge " << n << " " <<m<<" "<<generalized_m<< std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (const auto& value : adj[i]) {
			if(value>i)
            outFile << "e " << i+1 << " " << value+1 << std::endl;
        }
    }
    for (size_t i = 0; i < n; ++i) {
            outFile << "n " << i+1 << " " << weight_node[i] << std::endl;
    }
    for (size_t i = 0; i < n; ++i) {
        for (const auto& value : generalized_adj[i]) {
			if(value.first>i)
            outFile << "not_e " << i+1 << " " << value.first+1 << " " << value.second << std::endl;
        }
    }
    outFile.close();
	cout<<"writing successfully"<<endl;
}
void Graph::exert_reductions(char* is,int& res,int& valid_n){
	int last_turn=valid_n;
	int ids_n=valid_n;
	do{
		reduction_neighborhood(is,res,true);
		add_0d(is,res);
		reduction_penalty(is,res,true);
		add_0d(is,res);
		Negative_Profit(is,res,true);
		add_0d(is,res);
		isolated(is,res,true);
		add_0d(is,res);
		extend_single_edge(is,res,true);
		add_0d(is,res);
		last_turn=valid_n;
		ids_n=update_edges(is,res,valid_n);
		add_0d(is,res);
	}while(last_turn-valid_n>n/1000);
}
void Graph::Random_Peeling(char* is,int&res){
	exclusive_pairs.clear();
	int num=0;
	int valid_n=0;
	for(int i=0;i<n;i++) {
		if(is[i]==1){
		valid_n++;
		}
	}
	vector<ui> valid2n;
	valid2n.reserve(valid_n);
	while(valid_n/10>1){
		num++;
		valid2n.clear();
		for(int i=0;i<n;i++) {
			if(is[i]==1){
			valid2n.emplace_back(i);
			}
		}
		set<ui> random_pos;
		while(random_pos.size()<valid_n/10){
			random_pos.insert(get_rand_n(valid_n));
		}
		for(ui u:random_pos){
			delete_vertex(valid2n[u],res,is);
			valid_n--;
		}
		update_sets(is);
		exert_reductions(is,res,valid_n);
		if(valid_n/10<=1) break;
	}
}
int Graph::get_add_utility(ui u){
	int weight_Npu=0;
	int pe_u=0;
	for(auto& v:generalized_adj[u]){
		if(S.contains(v.first)){
			pe_u+=v.second;
		}
	}
	return weight_node[u]-pe_u;
}
int Graph::get_swap_utility(ui u,ui v,vector<int>& utility){
	//utility(u,v)=Bv-Bu
	return utility[v]-utility[u];
}

int Graph::Add_Set(int pos,int& res,int& current, priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>>& tabu,vector<ui>& is_tabu,
RandList& add_set,RandList& swap_set,vector<ui>& swap_target,vector<int>& utility,vector<ui>& hard_degree_in_S,char* is,int best_res){
	vector<ui> in_tabu;
	vector<ui> not_tabu;
	int max_in_tabu=INT32_MIN;
	int max_not_tabu=INT32_MIN;
	int max_vertex=-1;
	if(pos==-1){
		for(int i=0;i<add_set.vnum;i++){
			ui u=add_set.vlist[i];
			if(is_tabu[u]){
				if(utility[u]>max_in_tabu){
					in_tabu.clear();
					max_in_tabu=utility[u];
					in_tabu.emplace_back(u);
				}else if(utility[u]==max_in_tabu){
					in_tabu.emplace_back(u);
				}
			}else{
				if(utility[u]>max_not_tabu){
					not_tabu.clear();
					max_not_tabu=utility[u];
					not_tabu.emplace_back(u);
				}else if(utility[u]==max_not_tabu){
					not_tabu.emplace_back(u);
				}
			}
		}
		if(!in_tabu.empty() && max_in_tabu>max_not_tabu && max_in_tabu + res > best_res){
			max_vertex=in_tabu[get_rand_n(in_tabu.size())];
		}else if(!not_tabu.empty()){
			max_vertex=not_tabu[get_rand_n(not_tabu.size())];
		}else{
			return -1;
		}
	}else{
		max_vertex=pos;
	}
	S.add(max_vertex);
	add_set.remove(max_vertex);
	//更新Np_S
	for(ui v:adj[max_vertex]){
		hard_degree_in_S[v]++;
		if(hard_degree_in_S[v]==1){
			add_set.remove(v);
			swap_set.add(v);
			swap_target[v]=max_vertex;
		}else if(hard_degree_in_S[v]==2){
			swap_set.remove(v);
			swap_target[v]=-1;
		}
	}
	res+=weight_node[max_vertex];
	for(auto& v:generalized_adj[max_vertex]){
		utility[v.first]-=v.second;
		if(S.contains(v.first)){
			res-=v.second;
		}
		dynamic_weight[v.first]-=max(0,-v.second);
	}
	return 0;
}
void Graph::remove_from_S(int& res,ui u){
	S.remove(u);
	res-=weight_node[u];
}
void Graph::remove_from_S(int& res,ui u,vector<int>& utility){
	S.remove(u);
	res-=weight_node[u];
	for(auto& v:generalized_adj[u]){
		utility[v.first]+=v.second;
		if(S.contains(v.first)){
			res+=v.second;
		}
		dynamic_weight[v.first]-=max(0,-v.second);
	}
}
void Graph::Drop_Set(int pos,int& res,int current,priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>>& tabu,vector<ui>& is_tabu,RandList& add_set,RandList& swap_set,vector<ui>& swap_target,
 vector<int>& utility,vector<ui>& hard_degree_in_S, bool& can_drop,char* is){
	if(S.vnum==0) return;
	can_drop=true;
	int min_utility=INT32_MAX;
	int min_u=-1;
	vector<ui> min_list;
	if(pos==-1){
		for(int i=0;i<S.vnum;i++){
			ui u=S.vlist[i];
			int current_utility=utility[u];
			if(current_utility<min_utility){
				min_list.clear();
				min_list.emplace_back(u);
				min_utility=current_utility;
			}else if(current_utility==min_utility){
				min_list.emplace_back(u);
			}
		}
		if(false){
			min_u=S.vlist[get_rand_n(S.vnum)];
		}else{
			if(!min_list.empty())
			min_u=min_list[get_rand_n(min_list.size())];
		}
	}else{
		min_u=pos;
	}
	remove_from_S(res,min_u,utility);
	add_set.add(min_u);
	for(ui v:adj[min_u]){
		hard_degree_in_S[v]--;
		if(hard_degree_in_S[v]==0){
			swap_set.remove(v);
			swap_target[v]=-1;
			add_set.add(v);
		}else if(hard_degree_in_S[v]==1){
			swap_set.add(v);
		}
	}
	tabu.push(make_pair(current+max(1,min((int)S.vnum/10,10)),min_u));
	is_tabu[min_u]=1;
}
void Graph::write_kernel(Graph* kernel,vector<ui>& k2e,vector<ui>& c2k){
	size_t lastSlashPos = kernel->dir.find_last_of("/\\");
	std::string filename = kernel->dir.substr(lastSlashPos + 1);
	    size_t dotPos = filename.rfind(".txt");
    if (dotPos != std::string::npos) {
        filename = filename.substr(0, dotPos);
    }
    std::ofstream outFile("all_kernel/"+filename+"_kernel.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
 	outFile << "p " << "edge " << kernel->n << " " <<kernel->m<<" "<<kernel->generalized_m<< std::endl;
    for (size_t i = 0; i < kernel->n; ++i) {
        for (const auto& value : kernel->adj[i]) {
			if(value>i)
            outFile << "e " << i+1 << " " << value+1 << std::endl;
        }
    }
    for (size_t i = 0; i < kernel->n; ++i) {
            outFile << "n " << i+1 << " " << kernel->weight_node[i] << std::endl;
    }
    for (size_t i = 0; i < kernel->n; ++i) {
        for (const auto& value : kernel->generalized_adj[i]) {
			if(value.first>i)
            outFile << "not_e " << i+1 << " " << value.first+1 << " " << value.second << std::endl;
        }
    }
    outFile.close();
	std::ofstream outFile2("kernel/"+filename+"_mapping.txt");
	for(int i=0;i<kernel->n;i++){
		outFile2<<i<<" "<<c2k[i]<<endl;
	}
	outFile2.close();
	std::ofstream outFile3("kernel/"+filename+"_exc.txt");
	for(auto& myturple: exclusive_pairs){
		outFile3<<myturple.u<<" "<<myturple.x<<" "<<myturple.y<<" "<<myturple.fold_vertex<<" "<<myturple.state<<endl;
	}
	outFile3.close();
	std::ofstream outFile4("kernel/"+filename+"_k2e.txt");
	for(int i=0;i<n;i++){
		outFile4<<i<<" "<<k2e[i]<<endl;
	}
	outFile4.close();
	std::ofstream outFile5("kernel/"+filename+"_oS.txt");
	for(int i=0;i<S.vnum;i++){
		outFile5<<S.vlist[i]<<endl;
	}
	outFile5.close();
	cout<<"writing successfully"<<endl;
}
int Graph::Swap_Set(int pos,int& res,int current,priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>>& tabu,vector<ui>& is_tabu,RandList& swap_set,vector<ui>& swap_target,RandList& add_set,
 vector<int>& utility,vector<ui>& hard_degree_in_S,char* is,int best_res){
	if(swap_set.empty()) return -1;
	vector<ui> in_tabu;
	vector<ui> not_tabu;
	int max_in_tabu=INT32_MIN;
	int max_not_tabu=INT32_MIN;
	int max_u=-1;
	int max_v=-1;
	if(pos==-1){
		for(int i=0;i<swap_set.vnum;i++){
			ui v=swap_set.vlist[i];
			ui u=swap_target[v];
			int current_utility=get_swap_utility(u,v,utility);
			if(is_tabu[v]){
				if(current_utility>max_in_tabu){
					in_tabu.clear();
					max_in_tabu=current_utility;
					in_tabu.emplace_back(v);
				}else if(current_utility==max_in_tabu){
					in_tabu.emplace_back(v);
				}
			}else{
				if(current_utility>max_not_tabu){
					not_tabu.clear();
					max_not_tabu=current_utility;
					not_tabu.emplace_back(v);
				}else if(current_utility==max_not_tabu){
					not_tabu.emplace_back(v);
				}
			}
		}
		if(!in_tabu.empty() && max_in_tabu>max_not_tabu && max_in_tabu + res > best_res){
			int r=get_rand_n(in_tabu.size());
			max_v=in_tabu[r];
			max_u=swap_target[max_v];
		}else if(!not_tabu.empty()){
			int r=get_rand_n(not_tabu.size());
			max_v=not_tabu[r];
			max_u=swap_target[max_v];
		}else{
			return -1;
		}
	}else{
		max_v=pos;
		max_u=swap_target[max_v];
	}
	S.add(max_v);
	swap_set.remove(max_v);
	swap_target[max_v]=-1;
	remove_from_S(res,max_u,utility);
	swap_set.add(max_u);
	swap_target[max_u]=max_v;
	for(ui v:adj[max_v]){
		hard_degree_in_S[v]++;
		if(hard_degree_in_S[v]==1 && !S.contains(v)){
			add_set.remove(v);
			swap_set.add(v);
			swap_target[v]=max_v;
		}else if(hard_degree_in_S[v]==2){
			swap_set.remove(v);
			swap_target[v]=-1;
		}
	}
	for(ui v:adj[max_u]){
		hard_degree_in_S[v]--;
		if(hard_degree_in_S[v]==0 && !S.contains(v)){
			swap_set.remove(v);
			swap_target[v]=-1;
			add_set.add(v);
		}else if(hard_degree_in_S[v]==1){
			swap_set.add(v);
		}
	}
	res+=weight_node[max_v];
	for(auto& v:generalized_adj[max_v]){
		utility[v.first]-=v.second;
		if(S.contains(v.first)){
			res-=v.second;
		}
		dynamic_weight[v.first]-=max(0,-v.second);
	}
	tabu.push(make_pair(current+max(1,min(int(S.vnum)/10,10))+get_rand_n(swap_set.vnum),max_u));
	is_tabu[max_u]=1;
	return 0;
}
void Graph::cal_utility(int& add_pos,int& swap_pos,int& drop_pos, vector<int>& utility,vector<ui>& is_tabu,RandList& add_set, RandList& swap_set,vector<ui>& swap_target,int&res,int best_res){
	int min_utility=INT32_MAX;
	int min_u=-1;
	int max_utility=INT32_MIN;
	int max_u=-1;
	vector<ui> in_tabu;
	vector<ui> not_tabu;
	vector<ui> min_list;
	int max_in_tabu=INT32_MIN;
	int max_not_tabu=INT32_MIN;
	for(int i=0;i<add_set.vnum;i++){
		ui u=add_set.vlist[i];
		if(is_tabu[u]){
			if(utility[u]>max_in_tabu){
				in_tabu.clear();
				max_in_tabu=utility[u];
				in_tabu.emplace_back(u);
			}else if(utility[u]==max_in_tabu){
				in_tabu.emplace_back(u);
			}
		}else{
			if(utility[u]>max_not_tabu){
				not_tabu.clear();
				max_not_tabu=utility[u];
				not_tabu.emplace_back(u);
			}else if(utility[u]==max_not_tabu){
				not_tabu.emplace_back(u);
			}
		}
	}
	if(!in_tabu.empty() && max_in_tabu>max_not_tabu && max_in_tabu + res > best_res){
		max_u=in_tabu[get_rand_n(in_tabu.size())];
		max_utility=max_in_tabu;
	}else if(!not_tabu.empty()){
		max_u=not_tabu[get_rand_n(not_tabu.size())];
		max_utility=max_not_tabu;
	}
	if(add_set.empty() || max_u==-1)
	for(int i=0;i<S.vnum;i++){
		ui u=S.vlist[i];
		if(utility[u]<min_utility){
			min_list.clear();
			min_list.emplace_back(u);
			min_utility=utility[u];
			// min_u=u;
		}else if(utility[u]==min_utility){
			min_list.emplace_back(u);
		}
	}
	if(!min_list.empty())
	min_u=min_list[get_rand_n(min_list.size())];
	max_in_tabu=INT32_MIN;
	max_not_tabu=INT32_MIN;
	in_tabu.clear();
	not_tabu.clear();
	int max_swap=INT32_MIN;
	int max_swap_i=-1;
	for(int i=0;i<swap_set.vnum;i++){
		ui v=swap_set.vlist[i]; 
		ui u=swap_target[v];
		int current_utility=get_swap_utility(u,v,utility);
		if(is_tabu[v]){
			if(current_utility>max_in_tabu){
				in_tabu.clear();
				max_in_tabu=current_utility;
				in_tabu.emplace_back(v);
			}else if(current_utility==max_in_tabu){
				in_tabu.emplace_back(v);
			}
		}else{
			if(current_utility>max_not_tabu){
				not_tabu.clear();
				max_not_tabu=current_utility;
				not_tabu.emplace_back(v);
			}else if(current_utility==max_not_tabu){
				not_tabu.emplace_back(v);
			}
		}
	}

	if(!in_tabu.empty() && max_in_tabu>max_not_tabu && max_in_tabu + res > best_res){
		int r=get_rand_n(in_tabu.size());
		max_swap_i=in_tabu[r];
		max_swap=max_in_tabu;
	}else if(!not_tabu.empty()){
		int r=get_rand_n(not_tabu.size());
		max_swap_i=not_tabu[r];
		max_swap=max_in_tabu;
	}
	if(max_utility>-min_utility){
		if(max_utility>max_swap){
			add_pos=max_u;
		}else{
			swap_pos=max_swap_i;
		}
	}else{
		if(-min_utility>max_swap){
			drop_pos=min_u;
		}else{
			swap_pos=max_swap_i;
		}
	}
}

Graph* Graph::copy_kernel(){
	Graph* kernel=new Graph();
	kernel->dir=dir;
	kernel->n=n;
	kernel->hard_degreeInC.resize(n,0);
	kernel->soft_degreeInC.resize(n,0);
	kernel->weight_node.resize(n,0);
	kernel->dynamic_weight.resize(n,0);
	kernel->adj.resize(n);
	kernel->generalized_adj.resize(n);
	kernel->exclusive_pairs.reserve(exclusive_pairs.size()+1);
	for(int i=0;i<n;i++){
		kernel->hard_degreeInC[i]=hard_degreeInC[i];
		kernel->soft_degreeInC[i]=soft_degreeInC[i];
		for(ui v:adj[i]){
			if(v>i)
			kernel->add_edge(v,i);
		}
		for(auto& v:generalized_adj[i]){
			if(v.first>i)
			kernel->add_soft_edge(v.first,i,v.second);
		}
		kernel->weight_node[i]=weight_node[i];
	}
	kernel->soft_neighbors.init(n);
	kernel->hard_neighbors.init(n);
	kernel->hard_degree_ones.init(n);
	kernel->hard_degree_twos.init(n);
	kernel->total_degree_ones.init(n);
	kernel->total_degree_twos.init(n);
	kernel->temp_delete_hard_ones.init(n);
	kernel->temp_delete_hard_twos.init(n);
	kernel->temp_delete_total_ones.init(n);
	kernel->temp_delete_total_twos.init(n);
	kernel->S.init(n);
	for(int i=0;i<S.vnum;i++){
		kernel->S.add(S.vlist[i]);
	}
	return kernel;
}
Graph* Graph::induced_search_kernel(vector<ui>& sk2k,vector<ui>& k2sk,vector<ui>& k2e,vector<ui>& e2k,RandList& old_S, RandList& Np_S){
	Graph* origin_graph=new Graph();
	origin_graph->dir=dir;
	origin_graph->read_graph_GIS();
	Graph* kernel=new Graph();
	kernel->n=n;
	kernel->hard_degreeInC.resize(n,0);
	kernel->soft_degreeInC.resize(n,0);
	kernel->weight_node.resize(n,0);
	kernel->dynamic_weight.resize(n,0);
	kernel->adj.resize(n);
	kernel->generalized_adj.resize(n);
	for(int nu=0;nu<n;nu++){
		ui ou=k2e[sk2k[nu]];
		for(ui ov:origin_graph->adj[ou]){
			ui nv=k2sk[e2k[ov]];
			if(nv==-1) continue;
			if(nv>nu){
				kernel->add_edge(nu,nv);
			}
		}
		for(auto& ov:origin_graph->generalized_adj[ou]){
			ui nv=k2sk[e2k[ov.first]];
			if(nv==-1) continue;
			if(nv>nu){
				kernel->add_soft_edge(nu,nv,ov.second);
			}
		}
		kernel->weight_node[nu]=origin_graph->weight_node[ou];
	}
	for(int i=0;i<old_S.vnum;i++){
		ui ou=k2e[old_S.vlist[i]];
		for(auto& ov:origin_graph->generalized_adj[ou]){
			ui nv=k2sk[e2k[ov.first]];
			ui kv=e2k[ov.first];
			if(!Np_S.contains(kv)){
				kernel->weight_node[nv]-=ov.second;
			}
		}

	}
	kernel->soft_neighbors.init(n);
	kernel->hard_neighbors.init(n);
	kernel->hard_degree_ones.init(n);
	kernel->hard_degree_twos.init(n);
	kernel->total_degree_ones.init(n);
	kernel->total_degree_twos.init(n);
	kernel->temp_delete_hard_ones.init(n);
	kernel->temp_delete_hard_twos.init(n);
	kernel->temp_delete_total_ones.init(n);
	kernel->temp_delete_total_twos.init(n);
	kernel->S.init(n);
	for(int i=0;i<S.vnum;i++){
		kernel->S.add(S.vlist[i]);
	}
	return kernel;
}
Graph* Graph::induced_graph(vector<ui>& k2e,vector<ui>& e2k){
	Graph* origin_graph=new Graph();
	origin_graph->dir=dir;
	origin_graph->read_graph_GIS();
	Graph* kernel=new Graph();
	kernel->n=n;
	kernel->hard_degreeInC.resize(n,0);
	kernel->soft_degreeInC.resize(n,0);
	kernel->weight_node.resize(n,0);
	kernel->dynamic_weight.resize(n,0);
	kernel->adj.resize(n);
	kernel->generalized_adj.resize(n);
	kernel->exclusive_pairs.reserve(exclusive_pairs.size()+1);
	for(int nu=0;nu<n;nu++){
		ui ou=k2e[nu];
		for(ui ov:origin_graph->adj[ou]){
			ui nv=e2k[ov];
			if(nv>nu){
				kernel->add_edge(nu,nv);
			}
		}
		for(auto& ov:origin_graph->generalized_adj[ou]){
			ui nv=e2k[ov.first];
			if(nv>nu){
				kernel->add_soft_edge(nu,nv,ov.second);
			}
		}
		kernel->weight_node[nu]=origin_graph->weight_node[ou];
	}
	kernel->S.init(n);
	return kernel;
}
void Graph::check_maximal(vector<ui>& e2k){
	hard_neighbors.clear();
	for(int i=0;i<S.vnum;i++){
		ui u=S.vlist[i];
		for(ui v:adj[u]) hard_neighbors.add(v);
	}
	for(int i=0;i<n;i++){
		if(!S.contains(i) && !hard_neighbors.contains(i)){
			int sum=0;
			for(auto& w:generalized_adj[i]){
				if(S.contains(w.first)) sum+=w.second;
			}
		}
	}
}
void Graph::check(int& res_w){
	int res=0;
	bool is_IS=true;
	for(int i=0;i<S.vnum;i++){
		ui u=S.vlist[i];
		res+=weight_node[u];
		for(ui v:adj[u]){
			if(S.contains(v)) {
				is_IS=false;
				cout<<"u:"<<u<<" v:"<<v<<endl;
				}
		}
		if(!is_IS) break;
		for(auto& v:generalized_adj[u]){
			if(S.contains(v.first) && v.first>u) res-=v.second;
		}
	}
	if(!is_IS){
		cout<<"NOT AN IS"<<endl;
	}
}
int Graph::check(){
	int res=0;
	bool is_IS=true;
	for(int i=0;i<S.vnum;i++){
		ui u=S.vlist[i];
		res+=weight_node[u];
		for(ui v:adj[u]){
			if(S.contains(v)) {
				is_IS=false;
				cout<<"u:"<<u<<" v:"<<v<<endl;
				}
		}
		if(!is_IS) break;
		for(auto& v:generalized_adj[u]){
			if(S.contains(v.first) && v.first>u) res-=v.second;
		}
	}
	if(!is_IS){
		cout<<"NOT AN IS"<<endl;
	}
	return res;
}
void Graph::Negative_Profit(char* is,int& res,bool init){
	if(init){
		Cand_reduce.clear();
		for(ui i=0;i<n;i++){
			Cand_reduce.emplace_back(i);
		}
	}
	temp_cand.clear();
	for(ui i:Cand_reduce){
		if(is[i]!=1) continue;
		if(weight_node[i]<0){
			bool all_large=true;
			int sum=0;
			for(auto& v:generalized_adj[i]){
				if(is[v.first]!=1) continue;
				sum+=min(0,v.second);
				if(sum<=weight_node[i]){
						all_large=false;
						break;
				}
			}
			if(all_large){
				delete_vertex(i,res,is);
			for(ui w:adj[i]){
				if(is[w]==1)
				temp_cand.emplace_back(w);
			}
			for(auto& w:generalized_adj[i]){
				if(is[w.first]==1)
				temp_cand.emplace_back(w.first);
			}
			}
		}
	}

	Cand_reduce.clear();
	Cand_reduce=temp_cand;
}
void Graph::RLS(vector<ui>& k2e,vector<ui>& e2k, char* old_is,int& res){
	char *is = new char[n];
	for(ui i = 0;i < n;i ++) {
		if(S.contains(i)) is[i]=2;
		else if(old_is[k2e[i]]==3)
				is[i] = 3;
		else if(old_is[k2e[i]]==1)
				is[i] = 1;
		}
	for(ui i=0;i<n;i++){
		if(is[i]!=1) continue;
		if(hard_degreeInC[i]==2){
			hard_degree_twos.add(i);
			if(soft_degreeInC[i]==0){
				total_degree_twos.add(i);
			}
		}else if(hard_degreeInC[i]==1){
			hard_degree_ones.add(i);
			if(soft_degreeInC[i]==1){
				total_degree_twos.add(i);
			}else if(soft_degreeInC[i]==0){
				total_degree_ones.add(i);
			}
		}else if(hard_degreeInC[i]==0){
			if(soft_degreeInC[i]==2){
				total_degree_twos.add(i);
			}else if(soft_degreeInC[i]==1){
				total_degree_ones.add(i);
			}else if(soft_degreeInC[i]==0){
					add_vertex(i,res,is);
			}
		}
	}
	int best_res=res;
	srand(0);
	int search_times=10;
	double sum=0;
	int best_size=S.vnum;
	double average=0;
	double average_size=0;
	int best_turn=0;
	cout<<"search_times:"<<search_times<<endl;
	for(int i=0;i<search_times;i++){
		cout<<"turn:"<<i+1<<endl;
		double best_time=0;
		best_size=S.vnum;
		int temp = local_search(res,is,k2e,e2k,best_time,best_size);
		if(temp>best_res) {
			best_res=temp;
			best_turn=i+1;
			}
		sum+=best_time;
		average_size+=best_size;
		average+=temp;
	}
	cout<<"average_time: "<<sum/search_times<<endl;
	cout<<"average_res: "<<average/search_times<<endl;
	cout<<"average_size: "<<average_size/search_times<<endl;
	cout<<"best_trun: "<<best_turn<<endl;
	cout<<"best res: "<<best_res<<endl;
	cout<<"size of best res: "<<best_size<<endl;
	res=best_res;
	// S=best_S;
}
void Graph::recover_and_unfold(){
	size_t lastSlashPos = dir.find_last_of("/\\");
	std::string filename = dir.substr(lastSlashPos + 1);
	size_t dotPos = filename.rfind(".txt");
    if (dotPos != std::string::npos) {
        filename = filename.substr(0, dotPos);
    }
	std::ifstream file("kernel/"+filename+"_res.txt");
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return ;
    }
	vector<ui> temp_res;
	ui u;
	while(file>>u){
		temp_res.emplace_back(u);
	}
	std::sort(temp_res.begin(),temp_res.end());
	file.close();
	ui u_c,u_k;
	std::ifstream file2("kernel/"+filename+"_mapping.txt");
	int idx=0;
	while(file2>>u_c>>u_k){
		if(temp_res[idx]==u_c){
			idx++;
			S.add(u_k);
		}
	}
	file2.close(); 
	std::ifstream file3("kernel/"+filename+"_oS.txt");
	while(file3>>u){
		S.add(u);
	}
	file3.close();
	std::ifstream file4("kernel/"+filename+"_exc.txt");
	Myturple* myturple=new Myturple();
	while(file4>>myturple->u>>myturple->x>>myturple->y>>myturple->fold_vertex>>myturple->state){
		exclusive_pairs.emplace_back(myturple);
	}
	file4.close();
	unfold();
	std::ifstream file5("kernel/"+filename+"_k2e.txt");
	ui u_e;
	vector<ui> k2e;
	vector<ui> e2k;
	e2k.resize(n,n+1); 
	while(file5>>u_k>>u_e){
		k2e.emplace_back(u_e);
		e2k[u_e]=u_k;
	}
	file5.close();
	temp_res.clear();
	for(int i=0;i<S.vnum;i++) temp_res.emplace_back(S.vlist[i]);
	S.clear();
	for(ui v:temp_res){
		S.add(k2e[v]);
	}
	int res =check();
 	std::ofstream outFile("kernel_result/"+filename+".txt");
	outFile<<filename<<":"<<res<<endl;
}
void Graph::GIS() {
	char *is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;
	int res = 0;
	for(ui i=0;i<n;i++){
		if(hard_degreeInC[i]==2){
			hard_degree_twos.add(i);
			if(soft_degreeInC[i]==0){
				total_degree_twos.add(i);
			}
		}else if(hard_degreeInC[i]==1){
			hard_degree_ones.add(i);
			if(soft_degreeInC[i]==1){
				total_degree_twos.add(i);
			}else if(soft_degreeInC[i]==0){
				total_degree_ones.add(i);
			}
		}else if(hard_degreeInC[i]==0){
			if(soft_degreeInC[i]==2){
				total_degree_twos.add(i);
			}else if(soft_degreeInC[i]==1){
				total_degree_ones.add(i);
			}else if(soft_degreeInC[i]==0){
				d0.emplace_back(i);
			}
		}
	}
	add_0d(is,res);
	int last_turn=n;
	int ids_n=n;
	bool too_dense=false;
	if((m+generalized_m)/n>100) too_dense=true;
	auto start_time = std::chrono::high_resolution_clock::now();
	do{
		reduction_neighborhood(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			reduction_neighborhood(is,res,false); 
			add_0d(is,res);
		}
		reduction_penalty(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			reduction_penalty(is,res,false); 
			add_0d(is,res);
		}
		isolated(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			isolated(is,res,false); 
			add_0d(is,res);
		}
		Negative_Profit(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			Negative_Profit(is,res,false); 
			add_0d(is,res);
		}
		update_sets(is);
		while(!hard_degree_ones.empty()){
			reduction_hard_degree_one(is,res);
			add_0d(is,res);
		}
		reduction_hard_degree_two(is,res);
		add_0d(is,res);
		while(!total_degree_ones.empty()){
			reduction_degree_one(is,res);
			add_0d(is,res);
		}
		while(!total_degree_twos.empty()){
			reduction_triangle(is,res);
			add_0d(is,res);
		}
		if(!too_dense)
		basic_single_edge(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			basic_single_edge(is,res,false); 
			add_0d(is,res);
		}
		extend_single_edge(is,res,true); 
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			extend_single_edge(is,res,false); 
			add_0d(is,res);
		}
		twin(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			twin(is,res,false); 
			add_0d(is,res);
		}
		edge_domination(is,res,true);
		add_0d(is,res);
		while(!Cand_reduce.empty()){
			edge_domination(is,res,false); 
			add_0d(is,res);
		}
		last_turn=ids_n;
		ids_n=update_edges(is,res);
	}while(last_turn - ids_n>0);
	int temp=0;
	vector<ui> k2e,e2k;
	e2k.resize(n,n+1);
	Graph* kernel=new Graph();
	kernel->dir=dir;
	cout<<kernel->dir<<endl;
	kernel=generate_kernel(is,res,k2e,e2k); 
	int num_=0;
	int num_fold=0;
	// vector<ui> C;
	for(int i=0;i<kernel->n;i++){
		if(is[k2e[i]]==2){
			kernel->S.add(i);
			// cout<<k2e[i]<<" ";
		}else if(is[k2e[i]]==1){
			num_++;
			// C.emplace_back(k2e[i]);
		}else if(is[k2e[i]]==3){
			num_fold++;
		}
	}
	auto current_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = current_time - start_time;
	cout<<"genrating time："<<elapsed_time.count()<<endl;
	cout<<"S:"<<kernel->S.vnum<<endl;
	cout<<"kernel_size:"<<num_<<endl;
	cout<<"fold and exc:"<<num_fold<<endl;
	cout<<"current res after unfold:"<<res<<endl;
	kernel->RLS(k2e,e2k,is,res);
}

