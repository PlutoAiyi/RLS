#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include <list>
#include <chrono>
#include "RandList.hpp"
#include "Myturple.hpp"
#include <fstream>
#include <limits>
#include <random>
using namespace std;

struct Edge {
	int id, duplicate;
	int next;
};

class Graph {
public:
	string dir; //input graph directory
	ui n, m, generalized_m; //number of nodes and edges of the graph
	ui test;
	vector<int> weight_node;
	//dyamic_weight=weight_node+ sum max(0,-p(v,x))
	vector<int> dynamic_weight;
	vector<ui> d0;
	vector<ui> Cand_reduce;
	vector<ui> temp_cand;
	// list<tuple<ui,ui,int>> all_soft_edges;
	// ui *weight_edge;
	// ui *pstart; //offset of neighbors of nodes
	// ui *edges; //adjacent ids of edges
	vector<list<ui>> adj;
	vector<list<pair<ui,int>>> generalized_adj;
	// vector<list<pair<ui,int>>> virtual_adj;
	// vector<list<ui>> weight_edge;
	// vector<ui> hard_degree;
	// vector<ui> soft_degree;
	vector<ui> hard_degreeInC;
	vector<ui> soft_degreeInC;
	RandList hard_degree_ones, hard_degree_twos, total_degree_ones, total_degree_twos;
	RandList soft_neighbors;
	RandList hard_neighbors;
	RandList temp_delete_hard_ones, temp_delete_hard_twos, temp_delete_total_ones, temp_delete_total_twos;
	vector<ui> temp_add_hard_ones, temp_add_hard_twos, temp_add_total_ones, temp_add_total_twos;
	vector<Myturple> exclusive_pairs;
	RandList S;
	RandList temp_NP;
	// ui *is_generalized_edges;
	// ui *generalized_pstart; 
	// ui *generalized_edges; 
// public:
	Graph(const char *_dir) ;
	~Graph() ;
	Graph();
	void read_graph_GIS() ;
	void GIS() ;

private:
	void reduction_penalty(char* is,int& res,bool init);
	void recover_and_unfold();
	void write_kernel(Graph* kernel,vector<ui>& k2e,vector<ui>& c2k);
	void write_txt(string filename);
	void reduction_penalty4search(int& res,RandList& add_set,RandList& swap_set,vector<ui>& swap_target,vector<int>& utility,vector<ui>& hard_degree_in_S,char* is);
	void reduction_neighborhood4search(int& res,RandList& add_set,RandList& swap_set,vector<ui>& swap_target,vector<int>& utility,vector<ui>& hard_degree_in_S,char* is);
	void reduction_neighborhood(char* is,int& res,bool init);
	void reduction_degree_one(char* is,int& res);
	void reduction_hard_degree_one(char* is,int& res);
	void reduction_triangle(char* is,int& res);
	void add_0d(char* is,int&res);
	void check_maximal(vector<ui>& e2k);
	void unfold();
	void detect();
	Graph* copy_Cand(char* is,int& res,vector<ui>& c2k,vector<ui>& k2c);
	void reduction_hard_degree_two(char* is,int& res);
	void twin(char* is, int &res,bool init);
	int select_drop(vector<int>& utility,int& gain,int res);
	int select_add(RandList& add_set,vector<ui>& is_tabu,vector<int>& utility,int& res,int best_res,int& gain);
	int select_swap(RandList& swap_set,vector<ui>& swap_target,vector<ui>& is_tabu,vector<int>& utility,int& res,int best_res,int& gain);
	void cal_utility(int& add_pos,int& swap_pos,int& drop_pos, vector<int>& utility,vector<ui>& is_tabu,RandList& add_set, RandList& swap_set,vector<ui>& swap_target,int&res,int best_res);
	void basic_single_edge(char* is, int &res,bool init);
	void extend_single_edge(char* is, int &res,bool init);
	void edge_domination(char* is, int &res,bool init);
	void isolated(char* is,int& res,bool init);
	void RLS(vector<ui>& k2e,vector<ui>& e2k, char* old_is,int& res);
	int check();
	void check4java();
	void write_solution(string filename);
	void extend_S(char* is,int& res);
	void Negative_Profit(char* is,int& res,bool init);
	void check(int& res);
	bool has_hard_edge(char* is,ui u,ui v);
	int get_rand_n(int valid_n);
	void Random_Peeling(char* is,int& res);
	Graph* search_kernel(vector<ui>& sk2k,vector<ui>& k2sk,char* is);
	void exact_search(int& res,vector<ui>& k2e);
	void branch(int& res,ui u,RandList Np_S,RandList S,vector<ui>& k2e);
	int origin_search(int& res,vector<ui>& k2e,int depth,std::chrono::high_resolution_clock::time_point start_time,int run_time,double & best_time,int& global_best,bool need2reduction);
	int Add_Set(int pos,int& res,int& current,priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>>& tabu,vector<ui>& is_tabu,RandList& add_set,RandList& swap_set,vector<ui>& swap_target,vector<int>& utility,vector<ui>& hard_degree_in_S,char* is,int best_res);
	int get_add_utility(ui u);
	int get_swap_utility(ui u,ui v,vector<int>& utility);
	void update_swap_target(RandList& swap_set,vector<ui>& swap_target);
	void delete4fold(ui u,int& res,char* is);
	void Drop_Set(int pos,int& res,int current,priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>>& tabu,vector<ui>& is_tabu,RandList& add_set,RandList& swap_set,vector<ui>& swap_target, vector<int>& utility, vector<ui>& hard_degree_in_S, bool& can_drop,char* is);
	int Swap_Set(int pos,int& res,int current,priority_queue<pair<ui,ui>,vector<pair<ui,ui>>,greater<pair<ui,ui>>>& tabu,vector<ui>& is_tabu,RandList& swap_set,vector<ui>& swap_target,RandList& add_set, vector<int>& utility, vector<ui>& hard_degree_in_S,char* is,int best_res);
	void remove_from_S(int& res,ui u,vector<int>& utility);
	void remove_from_S(int& res,ui u);
	Graph* generate_kernel(char* is,int& res,vector<ui>& k2e,vector<ui>& e2k);
	int local_search(int res,char* is,vector<ui>& k2e,vector<ui>& e2k,double& best_time,int& best_size);
	void exert_reductions(char* is,int& res,int& valid_n);
	void soft_degree_increaseOne(char* is, ui v, int &res);
	void soft_degree_decreaseOne(char* is, ui v, int &res);
	void hard_degree_decreaseOne(char* is, ui v, int &res);
	void hard_degree_increaseOne(char* is, ui v, int &res);
	void has_hard(ui v,ui x,char* is,int puv,int &res_w);
	void has_hard(ui v,ui x,ui y,char* is,int puv,int& res_w);
	int update_edges(char* is, int &res);
	int update_edges(char* is, int &res,int& valid_n);
	int update_penalty(ui u,ui v, int add_num,char* is,int& res);
	void update_sets();
	void update_sets(char* is);
	Graph* copy_kernel();
	Graph* induced_graph(vector<ui>& k2e,vector<ui>& e2k);
	Graph* induced_search_kernel(vector<ui>& sk2k,vector<ui>& k2sk,vector<ui>& k2e,vector<ui>& e2k,RandList& old_S, RandList& Np_S);
	void add_vertex(ui u,int& res,char* is);
	void add_vertex4fold(ui u,int& res,char* is);
	void delete_vertex(ui u,int& res,char* is);
	void add_edge(ui v,ui u);
	void remove_edge(ui v,ui u);
	void add_soft_edge(ui v,ui u,int weight);
	void remove_soft_edge(ui v,ui u,int weight);
	void soft2hard(char* is,int& res);

};

#endif
