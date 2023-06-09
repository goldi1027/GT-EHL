/*
An Experimental Study on Hub Labeling based Shortest Path Algorithms [Experiments and Analyses]

Authors: Ye Li, Leong Hou U, Man Lung Yiu, Ngai Meng Kou
Contact: yb47438@umac.mo
Affiliation: University of Macau

The MIT License (MIT)

Copyright (c) 2016 University of Macau

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
#pragma once
#ifndef COVER_ORDERING_PATH_H
#define COVER_ORDERING_PATH_H

#include<algorithm>
#include<unordered_set>
#include <cmath>
//#include<unordered_map>
#include<time.h>
#include "graph.h"
#include "graph_search.h"
#include "time_util.h" 
#include "heap.h"
#include "paras.h"
#include "labels.h"

#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define CNT 16
#define cover_value_type long long
#include<unordered_map>
#include<cmath>

using namespace time_util;
  
class COrdering_P {
public:
	vector<NodeID> inv; // Fetch the original vertex id by a given ranking.
	vector<NodeID> rank; // Fetch the ranking of a given vertex id.
	
	void Relabel(Graph& graph) {
		
		for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;
		
		// Array Representation
		vector<EdgeID> new_vertices(numOfVertices + 1);
		vector<NodeID> new_edges;
		new_edges.reserve(graph.edges.size());
		for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
			NodeID originalVertex = inv[ranking];
			for (EdgeID eid = graph.vertices[originalVertex]; eid < graph.vertices[originalVertex + 1]; ++eid)
				new_edges.push_back(rank[graph.edges[eid]]);
			new_vertices[ranking + 1] = new_edges.size();
		}
		graph.vertices.swap(new_vertices);
		graph.edges.swap(new_edges);

		if (DIRECTED_FLAG == true) {
			vector<EdgeID> r_new_vertices(numOfVertices + 1);
			vector<NodeID> r_new_edges;
			r_new_edges.reserve(graph.r_edges.size());
			for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
				NodeID originalVertex = inv[ranking];
				for (EdgeID eid = graph.r_vertices[originalVertex]; eid < graph.r_vertices[originalVertex + 1]; ++eid)
					r_new_edges.push_back(rank[graph.r_edges[eid]]);
				r_new_vertices[ranking + 1] = r_new_edges.size();
			}
			graph.r_vertices.swap(r_new_vertices);
			graph.r_edges.swap(r_new_edges);
		}

	}

	void Relabel(WGraph& wgraph) {

		for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;

		/*vector<vector<NodeID> > new_adj(numOfVertices);
		vector<vector<EdgeWeight> > new_adj_weight(numOfVertices);
		vector<vector<NodeID> > new_r_adj(numOfVertices);
		vector<vector<EdgeWeight> > new_r_adj_weight(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < wgraph.adj[v].size(); ++i) {
				new_adj[rank[v]].push_back(rank[wgraph.adj[v][i]]);
				new_adj_weight[rank[v]].push_back(wgraph.adj_weight[v][i]);
			}
			if (DIRECTED_FLAG == true) {
				for (NodeID i = 0; i < wgraph.r_adj[v].size(); ++i) {
					new_r_adj[rank[v]].push_back(rank[wgraph.r_adj[v][i]]);
					new_r_adj_weight[rank[v]].push_back(wgraph.r_adj_weight[v][i]);
				}
			}
		}
		wgraph.adj.swap(new_adj);
		wgraph.adj_weight.swap(new_adj_weight);
		if (DIRECTED_FLAG == true) {
			wgraph.r_adj.swap(new_r_adj);
			wgraph.r_adj_weight.swap(new_r_adj_weight);
		}*/

		// Array Representation
		vector<EdgeID> new_vertices(numOfVertices + 1);
		vector<NodeEdgeWeightPair> new_edges;
		new_edges.reserve(wgraph.edges.size());
		for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
			NodeID originalVertex = inv[ranking];
			for (EdgeID eid = wgraph.vertices[originalVertex]; eid < wgraph.vertices[originalVertex + 1]; ++eid)
				new_edges.push_back(make_pair(rank[wgraph.edges[eid].first],wgraph.edges[eid].second));
			new_vertices[ranking + 1] = new_edges.size();
		}
		wgraph.vertices.swap(new_vertices);
		wgraph.edges.swap(new_edges);

		if (DIRECTED_FLAG == true) {
			vector<EdgeID> r_new_vertices(numOfVertices + 1);
			vector<NodeEdgeWeightPair> r_new_edges;
			r_new_edges.reserve(wgraph.r_edges.size());
			for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
				NodeID originalVertex = inv[ranking];
				for (EdgeID eid = wgraph.r_vertices[originalVertex]; eid < wgraph.r_vertices[originalVertex + 1]; ++eid)
					r_new_edges.push_back(make_pair(rank[wgraph.r_edges[eid].first],wgraph.r_edges[eid].second));
				r_new_vertices[ranking + 1] = r_new_edges.size();
			}
			wgraph.r_vertices.swap(r_new_vertices);
			wgraph.r_edges.swap(r_new_edges);
		}
	} 

	void ReswapLabel(Graph& graph) {

		vector<vector<NodeID> > new_adj(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < graph.adj[v].size(); ++i) {
				new_adj[v].push_back(inv[graph.adj[rank[v]][i]]);
			}
		}
		graph.adj.swap(new_adj);
	}

	int get_num_vertices(){
	    return numOfVertices;
	}

	double get_avg_degree(WGraph& wgraph, double &avg_degree, double &max_degree){
        vector<NodeID> deg_inv;
        vector<NodeID> deg_rank;
        deg_inv.resize(numOfVertices);
        deg_rank.resize(numOfVertices);
        double num_degree = 0;
        double temp_degree;
        vector<pair<float, NodeID> > deg(numOfVertices);
        for (size_t v = 0; v < numOfVertices; ++v) {
            deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) + float(random()) / RAND_MAX, v);
            num_degree += wgraph.vertices[v + 1] - wgraph.vertices[v];
            temp_degree = wgraph.vertices[v+1] - wgraph.vertices[v];

            if(temp_degree > max_degree){
                max_degree = temp_degree;
            }
        }
        avg_degree = num_degree/numOfVertices;

	}

	void save_rank(const char* order_file) {
		ofstream ofs(order_file);
		for (int i = 0; i < numOfVertices; ++i) {
			ofs << inv[i] << endl;
		}
		ofs.close();
	}

	~COrdering_P() {
		inv.clear();
		rank.clear();
	}

};



class Coverage_Ordering_Path : public COrdering_P {
	typedef	vector<NodeID> tree;
	public:
		PLabel plabels;
		NodeID last_available;
		
		//tree parent_tree: while a vertex v is not in the tree, parent_tree[v] = numOfVertices
	NodeID labeling_source_dij(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID>& visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parents, vector<NodeID>& parents, WGraph& wgraph){
			
		descendants.clear();
		NodeID visited_arcs = 0;
		
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

		//initialise all distance to inf weight
		for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);
		//distance to itself is 0
		distances[source] = 0;
		parent_tree[source] = numOfVertices;
		root_hop[source] = 0;
		parents[source] = source; 
		
		while (!pqueue.empty()) {

			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
			pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parents[v];
			
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;
			for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
				NodeID w = tmp_idx_v.first[i];
				EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
				if (td <= v_d) {
					goto pruned;
				}
			}

			// Traverse
			tmp_idx_v.first.back() = ranking;
			tmp_idx_v.second.back() = v_d;
			tmp_idx_v.first.push_back(numOfVertices);
			tmp_idx_v.second.push_back(INF_WEIGHT);
			
			
			//current source becomes parent of respective nodes
			tmp_idx_parent_v.first.back() = ranking;
			tmp_idx_parent_v.second.back() = parents[v];
			tmp_idx_parent_v.first.push_back(numOfVertices);
			tmp_idx_parent_v.second.push_back(numOfVertices);
			
			descendants.push_back(v);
			visited_arcs++;

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;
				if (!vis[w]) {
					if( distances[w] > w_d ){
						pqueue.update(w, w_d);
						distances[w] = w_d;
						parent_tree[w] = v;
						root_hop[w] = root_hop[v] + 1;
						parents[w] = v;
					}
				}
			}
			pruned: 
				{}
			}

		    //initialise and clear queue for next source node operation
			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				pqueue.clear(vis_v);
				parents[vis_v] = numOfVertices;
			}

			pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[source] = true;
			
			return visited_arcs;
	}
	
	
	NodeID labeling_source_dij_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants, vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID>& visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r,vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parents, vector<NodeID>& parents,vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_parents, vector<NodeID>& r_parents,WGraph& wgraph){
			
		descendants.clear();
		r_descendants.clear();
		NodeID visited_arcs = 0;
		
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

		for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);
		distances[source] = 0;		
		parent_tree[source] = numOfVertices;
		root_hop[source] = 0;
		parents[source] = source; 
		
		while (!pqueue.empty()) {

			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
			pair<vector<NodeID>, vector<NodeID> > &r_tmp_idx_parent_v = r_tmp_idx_parents[v];
			vis[v] = true;
			visited_que.push(v);
 
			if (usd[v]) continue;
			for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
				NodeID w = r_tmp_idx_v.first[i];
				EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
				if (td <= v_d) {
					goto pruned_forward;
				}
			}

			// Traverse
			r_tmp_idx_v.first.back() = ranking;
			r_tmp_idx_v.second.back() = v_d;
			r_tmp_idx_v.first.push_back(numOfVertices);
			r_tmp_idx_v.second.push_back(INF_WEIGHT);
			
			

			r_tmp_idx_parent_v.first.back() = ranking;
			r_tmp_idx_parent_v.second.back() = parents[v];
			r_tmp_idx_parent_v.first.push_back(numOfVertices);
			r_tmp_idx_parent_v.second.push_back(numOfVertices);
			
			descendants.push_back(v);
			visited_arcs++;

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;
				if (!vis[w]) {
					if( distances[w] > w_d ){
						pqueue.update(w, w_d);
						distances[w] = w_d;
						parent_tree[w] = v;
						root_hop[w] = root_hop[v] + 1;
						parents[w] = v;
					}
				}
			}
			pruned_forward: 
				{}
			}

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				pqueue.clear(vis_v); 
				parents[vis_v] = numOfVertices;
			}

			pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			pqueue.update(source, 0);
			distances[source] = 0;		
			r_parents[source] = source;
			
			// reverse search
			while (!pqueue.empty()) {

				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parents[v];
				
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_backward;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);
				
				
				tmp_idx_parent_v.first.back() = ranking;
				tmp_idx_parent_v.second.back() = r_parents[v];
				tmp_idx_parent_v.first.push_back(numOfVertices);
				tmp_idx_parent_v.second.push_back(numOfVertices);
				
				r_descendants.push_back(v);
				visited_arcs++;

				for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
					NodeID w = wgraph.r_edges[eid].first;
					EdgeWeight w_d = wgraph.r_edges[eid].second + v_d;
					if (!vis[w]) {
						if( distances[w] > w_d ){
							pqueue.update(w, w_d);
							distances[w] = w_d;
							r_parent_tree[w] = v;
							r_root_hop[w] = root_hop[v] + 1;
							r_parents[w] = v;
						}
					}
				}
				pruned_backward: 
					{}
			}

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				pqueue.clear(vis_v);
				r_parents[vis_v] = numOfVertices;
			}

			pqueue.clear_n();

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
			
			usd[source] = true;
			
			return visited_arcs;
	}
	
	NodeID labeling_source_bfs(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx,vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parents, vector<NodeID>& parents, Graph& graph) { 
			descendants.clear();
			NodeID visited_arcs = 0;
			
			int max_hop = -1;
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
					
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}
			
			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			parent_tree[source] = numOfVertices;
			root_hop[source] = 0;
			last_hop[source] = 0;
			parents[source] = source; 

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					
					
					NodeID v = que[que_i];
					if (usd[v]) continue; 

					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parents[v];
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}
 
				//	if(source == 28674)
				//		cout << v << " ," << tmp_idx_v.first.size()  << "," << descendants.size() << endl;

					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;

					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					
					
					tmp_idx_parent_v.first.back() = ranking;
					tmp_idx_parent_v.second.back() = parents[v];
					tmp_idx_parent_v.first.push_back(numOfVertices);
					tmp_idx_parent_v.second.push_back(numOfVertices);
					
					descendants.push_back(v);
					visited_arcs++;
					 
				//	if(source == 28674)
				//		cout << "2" << endl;
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
	
						NodeID w = graph.edges[eid];
						if (!vis[w]) {											
							que[que_h++] = w;
							parent_tree[w] = v;
							root_hop[w] = root_hop[v] + 1;
							last_hop[w] = root_hop[w];
							if(max_hop < root_hop[w])
								max_hop = root_hop[w];
							vis[w] = true;
							parents[w] = v;
						}
					}
				//	if(source == 28674)
				//		cout << "3" << endl;
				pruned_forward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			for (size_t i = 0; i < que_h; ++i){
				vis[que[i]] = false;
				parents[que[i]] = numOfVertices;
			}
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			}
			
			usd[source] = true;
			return visited_arcs;
	}
	
	NodeID labeling_source_bfs_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants,vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parents, vector<NodeID>& parents,vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_parents, vector<NodeID>& r_parents, Graph& graph) { 
			descendants.clear();
			r_descendants.clear();
			NodeID visited_arcs = 0;
			
			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			parent_tree[source] = numOfVertices;
			root_hop[source] = 0;
			parents[source] = source; 

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &r_tmp_idx_parent_v = r_tmp_idx_parents[v];
					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}

					// Traverse
					r_tmp_idx_v.first.back() = ranking;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
					
					r_tmp_idx_parent_v.first.back() = ranking;
					r_tmp_idx_parent_v.second.back() = parents[v];
					r_tmp_idx_parent_v.first.push_back(numOfVertices);
					r_tmp_idx_parent_v.second.push_back(numOfVertices);
					
					descendants.push_back(v);
					
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
							parent_tree[w] = v;
							root_hop[w] = root_hop[v] + 1;
							parents[w] = v;
						}
					}
				pruned_forward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			for (size_t i = 0; i < que_h; ++i){
				vis[que[i]] = false;
				parents[que[i]] = numOfVertices;
			}
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;

 
			// Backward search. 
			// Initialize backward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			r_parent_tree[source] = numOfVertices;
			r_root_hop[source] = 0;
			r_parents[source] = source;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parents[v];
					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the backward labels of r and forward labels of v in the backward search from r when reaching v (v->r path).
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
						if (td <= d) {
							goto pruned_backward;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					
					tmp_idx_parent_v.first.back() = ranking;
					tmp_idx_parent_v.second.back() = r_parents[v];
					tmp_idx_parent_v.first.push_back(numOfVertices);
					tmp_idx_parent_v.second.push_back(numOfVertices);
				
					r_descendants.push_back(v);

					// Array Representation
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
						NodeID w = graph.r_edges[eid];

						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
							r_parent_tree[w] = v;
							r_root_hop[w] = r_root_hop[v] + 1;
							r_parents[w] = v;
						}
					}
				pruned_backward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			for (size_t i = 0; i < que_h; ++i){
				vis[que[i]] = false;
				r_parents[que[i]] = numOfVertices;
			}
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
 
			usd[source] = true;
			return visited_arcs;
		}
	
	
		
	void calcover(vector<NodeID>& descendants, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& root_hop, vector<NodeID>& last_alive, vector<NodeID>& depth){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
			coverage[dv]++;
			//acc[dv]++;
			if(parent_tree[dv]!= numOfVertices){
				coverage[parent_tree[dv]] += coverage[dv];	
				//acc[parent_tree[dv]]++;	
				if(depth[dv] < root_hop[dv])
					depth[dv] = root_hop[dv];
				if(depth[dv] > root_hop[parent_tree[dv]]){
					depth[parent_tree[dv]] = depth[dv];
				}
				last_alive[dv] = coverage[dv];
			}
		}
		return;
	}
	
	void calcover(vector<NodeID>& descendants, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& root_hop){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			coverage[dv]++;
			if(parent_tree[dv]!= numOfVertices){
				coverage[parent_tree[dv]] += coverage[dv];	
			}
		}
		return;
	}
	
	
	void clear_tmp(vector<NodeID>& descendants, vector<NodeID>& coverage, tree& parent_tree, vector<NodeID>& root_hop, vector<NodeID>& depth){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
			coverage[dv] = 0;
			if(parent_tree[dv]!= numOfVertices)
				coverage[parent_tree[dv]] = 0;	
			parent_tree[dv] = numOfVertices;
			root_hop[dv] = 0;
			depth[dv] = 0;
		}
		descendants.clear();
		return;
	}
	
	double normalize(NodeID mina, NodeID maxa, NodeID minb, NodeID maxb, NodeID a, NodeID b){
		double a_score = (double) a * (double)mina / (double) maxa;
		double b_score = (double) b * (double)minb / (double) maxb;
		return a_score + b_score;
	}
	
	
		
	
	vector<NodeID> findSigPath(NodeID source, vector<NodeID>& coverage, tree& parent_tree, vector<bool>& usd,vector<NodeID>& upwardsPow, WGraph& wgraph, NodeID& maxdegree){
		vector<NodeID> sigpath;		
		NodeID v = source;
		upwardsPow.clear();
		maxdegree = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			NodeID max_v_degree;
			bool non_leaf_flag = false;
			//iterate over all connected edges to source node, finds the descendant with most visited nodes
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				if (parent_tree[w] == v && usd[w] == false) {
					non_leaf_flag = true;
					if(max_cover < coverage[w]){
						max_cover = coverage[w];
						max_v = w;
						max_v_degree = wgraph.vertices[w+1] - wgraph.vertices[w];
					}
				/*	if(max_cover < abs(coverage[v] - coverage[w]) * (wgraph.vertices[w + 1] - wgraph.vertices[w]) ){
						max_cover = abs(coverage[v] - coverage[w]) * (wgraph.vertices[w + 1] - wgraph.vertices[w]);
						max_v = w;
					}*/
				}					
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			if(maxdegree < max_v_degree)
				maxdegree = max_v_degree;
			sigpath.push_back(v);
		}		
		//cout << "sig path len:" << sigpath.size() << endl;
		return sigpath;
	}
	
	vector<NodeID> findRevSigPath(NodeID source, vector<NodeID>& coverage, vector<NodeID>& r_coverage, tree& parent_tree,tree& r_parent_tree, vector<bool>& usd, WGraph& wgraph, NodeID& maxdegree, bool& reverse_flag){
		vector<NodeID> sigpath;		
		NodeID v = source;
		NodeID max_degree = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree < (wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]))
				max_degree = (wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){
					NodeID w = wgraph.edges[eid].first;
					if (parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						NodeID c_coverage = coverage[w] * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);// + r_coverage[w];//(wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]); //coverage[w] + r_coverage[w]						
						if(max_cover < c_coverage){
							max_cover = c_coverage;
							max_v = w;
						}
					}					
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}
		
		vector<NodeID> rev_sigpath;		 
		v = source;
		NodeID max_degree_rev = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree_rev < (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) * (wgraph.vertices[v + 1] - wgraph.vertices[v]))
				max_degree_rev = (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) * (wgraph.vertices[v + 1] - wgraph.vertices[v]);
			for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){
					NodeID w = wgraph.r_edges[eid].first;
					if (r_parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						NodeID c_coverage = r_coverage[w]* (wgraph.vertices[v + 1] - wgraph.vertices[v]);// + coverage[w];// (wgraph.vertices[v + 1] - wgraph.vertices[v]); //r_coverage[w] + coverage[w];
						if(max_cover < c_coverage){
							max_cover =  c_coverage;
							max_v = w;
						}
					}					
			} 
			if(non_leaf_flag == false) break;
			v = max_v;
			rev_sigpath.push_back(v);
		}
		
		maxdegree = max_degree;
		reverse_flag = false;
		
		if(max_degree < max_degree_rev){
			reverse(rev_sigpath.begin(), rev_sigpath.end());
			sigpath = rev_sigpath;
			maxdegree = max_degree_rev;
			reverse_flag = true;
		} /*
		if(coverage[0] < r_coverage[0]){
			reverse_flag = true;
		}*/
		
		return sigpath;
	}
	
	vector<NodeID> findSigPath(NodeID source, vector<NodeID>& coverage, tree& parent_tree, vector<bool>& usd, vector<NodeID>& upwardsPow, Graph& graph){
		vector<NodeID> sigpath;		
		upwardsPow.clear();
		NodeID v = source;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;			
			bool non_leaf_flag = false;
			NodeID childrenNeighbors = 0;
			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
					NodeID w = graph.edges[eid];
					if (parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						childrenNeighbors++;
						if(max_cover < coverage[w]){
							max_cover = coverage[w];
							max_v = w;
						}/*
						if(max_cover < graph.vertices[w + 1] - graph.vertices[w]){
							max_cover = graph.vertices[w + 1] - graph.vertices[w];
							max_v = w;
						}*/
					}					
			}
			/*if( v != source){
				if(childrenNeighbors == 0)
					upwardsPow.push_back(0);
				else 
					upwardsPow.push_back( (coverage[v] / childrenNeighbors ) * ( graph.vertices[v + 1 ] - graph.vertices[v] - childrenNeighbors)) ;
			}*/
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}		
		//cout << "sig path len:" << sigpath.size() << endl;
		return sigpath;
	}
	
	vector<NodeID> findRevSigPath(NodeID source, vector<NodeID>& coverage, vector<NodeID>& r_coverage, tree& parent_tree,tree& r_parent_tree, vector<bool>& usd, Graph& graph){
		vector<NodeID> sigpath;		
		NodeID v = source;
		NodeID max_degree = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree < (graph.vertices[v + 1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]))
				max_degree = (graph.vertices[v + 1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]);
			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
					NodeID w = graph.edges[eid];
					if (parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						if(max_cover < coverage[w] + r_coverage[w]){
							max_cover = coverage[w] + r_coverage[w];
							max_v = w;
						}
					}					
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}
		
		vector<NodeID> rev_sigpath;		
		v = source;
		NodeID max_degree_rev = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree_rev < (graph.r_vertices[v + 1] - graph.r_vertices[v]) * (graph.vertices[v + 1] - graph.vertices[v]))
				max_degree_rev = (graph.r_vertices[v + 1] - graph.r_vertices[v]) * (graph.vertices[v + 1] - graph.vertices[v]);
			for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid){
					NodeID w = graph.r_edges[eid];
					if (r_parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						if(max_cover < r_coverage[w] + coverage[w]){
							max_cover =  r_coverage[w] + coverage[w];
							max_v = w;
						}
					}					
			} 
			if(non_leaf_flag == false) break;
			v = max_v;
			rev_sigpath.push_back(v);
		}
		/*
		cout << "sig path len:" << sigpath.size() << endl;
		NodeID maxd = -1;
		NodeID maxcd = -1;
		NodeID maxd_choice = -1;
		NodeID maxcd_choice = -1;
		if(sigpath.size()!=0){
			for(NodeID i = 0; i < sigpath.size(); ++i){
			//	cout << i << ":" << (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) << " ";
				if(maxd<(graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]])){
					maxd_choice = i;
					maxd = (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]);
				}
			}
			cout << endl;				
			for(NodeID i = 0; i < sigpath.size() - 1; ++i){
				//cout << i << ":" << (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]]) << " ";
				if(maxcd < (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]])){
					maxcd = (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]]);
					maxcd_choice = i;
				}
			}
		//	cout << endl;
			//cout << maxd_choice << " vs. " << maxcd_choice << endl;

			//for(NodeID i = 0; i < sigpath.size() - 1; ++i)
			//	cout << i << ":" << (graph.vertices[sigpath[i+1]+1] - graph.vertices[sigpath[i+1]]) *(graph.r_vertices[sigpath[i+1]+1] - graph.r_vertices[sigpath[i+1]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]]) << " ";
			//cout << endl;
		}
		
		cout << "rev sig path len:" << rev_sigpath.size() << endl;
		maxd_choice = -1;
		maxcd_choice = -1;
		NodeID rev_maxd = -1;
		NodeID rev_maxcd = -1;
		if(rev_sigpath.size()!=0){			
			for(NodeID i = 0; i < rev_sigpath.size(); ++i){
				//cout << i << ":" << (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) << " ";
				if(rev_maxd <  (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]])){
					rev_maxd =  (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]);
					maxd_choice = i;
				}
			}
			cout << endl;			
			//for(NodeID i = 0; i < rev_sigpath.size() - 1; ++i)
			//	cout << i << ":" << abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]) << " ";
			//cout << endl;		
			for(NodeID i = 0; i < rev_sigpath.size() - 1; ++i){
				//cout << i << ":" << (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) * abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]) << " ";
				if(rev_maxcd <(graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) * abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]) ){
					rev_maxcd = (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) * abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]);
					maxcd_choice = i;
				}
			}
		//	cout << endl;
		//	cout << maxd_choice << " vs. " << maxcd_choice << endl;
		}
		*/
		
		
		if(max_degree < max_degree_rev){
			reverse(rev_sigpath.begin(), rev_sigpath.end());
			sigpath = rev_sigpath;
		} 
		
		return sigpath;
	}
	
	
	NodeID calupwardFanout(vector<NodeID>& descendants, vector<NodeID>& upwardcoverage, vector<NodeID>& coverage, tree& parent_tree, vector<NodeID>& root_hop, Graph& graph){
		for(NodeID i = 0; i < descendants.size(); ++i){
			NodeID v = descendants[i];
			if(root_hop[v] == 0) continue;
			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
				NodeID w = graph.edges[eid];
				if (root_hop[v] == 1 && root_hop[w] != 0 ) {
					upwardcoverage[v]++;
				}

				if(root_hop[v] > 1){
					if(root_hop[w] <= root_hop[v])
						upwardcoverage[v]++;
				}
				
				/*if (root_hop[v] > root_hop[w] && root_hop[w] != 0)
					upwardcoverage[w] += upwardcoverage[v];	
				*/
			}				
		}
	}
	
	void clear_upwardtmp(vector<NodeID>& descendants, vector<NodeID>& upwardcoverage){
		for(NodeID i = 0; i < descendants.size(); ++i){
			upwardcoverage[descendants[i]] = 0;
		}
	}
	
	NodeID calupwardTree(NodeID v, tree& parent_tree, vector<NodeID>& root_hop, Graph& graph){
		NodeID calup = 0;
		for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
			NodeID w = graph.edges[eid];
			if (parent_tree[w] != v) {
				calup += root_hop[w] + 1;
			}
		}
		return calup;
	}
	
	void candidate_gen(vector<NodeID>& candidate,vector<NodeID>& inv, vector<bool>& usd){
		NodeID candidate_size = 300;
		candidate.clear();
		for(NodeID i = last_available; i < numOfVertices; ++i){
			NodeID v = inv[i];
			// Candidate
			if(usd[v] == false)
				candidate.push_back(v);
			if(candidate.size() > candidate_size - 1)
				break;
		}
		
		while(usd[inv[last_available]] == true)
			last_available++;
		
		return;
	}
	
	NodeID topcan(vector<NodeID>& inv, vector<bool>& usd, NodeID& last_available){
		while(true){
			NodeID v = inv[last_available];
			if(usd[v] == false)
				return v;
			else
				last_available++;
		}
	}
	
	void walk_stats(Graph& graph, vector<NodeID> border){
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
				
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
				deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) + float(random()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> last_alive(numOfVertices, 0);		
		vector<NodeID> last_hop(numOfVertices, 0);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		
		vector<NodeID> second_parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> second_root_hop(numOfVertices, 0);	
		vector<NodeID> second_coverage(numOfVertices, 0);
		vector<NodeID> second_depth(numOfVertices, 0);
		vector<NodeID> second_last_alive(numOfVertices, 0);		
		vector<NodeID> second_last_hop(numOfVertices, 0);
		vector<NodeID> second_descendants;
		second_descendants.reserve(numOfVertices);
				
		vector<NodeID> max_hop_set;
		
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
				
		vector<NodeID> acc_count;
		acc_count.resize(numOfVertices);
		
		vector<NodeID> que(numOfVertices);
		vector<NodeID> que2(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> usd(numOfVertices, false);
// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
		NodeID chosen = deg_inv[0];
		
		NodeID taken = numOfVertices; // /10;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<pair<vector<NodeID>, vector<NodeID> > >
			tmp_idx_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices)));
				
				
		vector<NodeID> parents(numOfVertices, numOfVertices);
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		for(NodeID i = 0; i < taken; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, tmp_idx_parents, parents, graph);
			currentsum += actual_cover;
			/*
			if(last_iter_hop > 0){
				actual_stats[last_iter_hop] += actual_cover;
				over_stats[last_iter_hop] += last_iter_cover;
			}
			cout << last_iter_cover << " vs. " << actual_cover << " on " << last_iter_hop << " hops" << endl;
			*/
			if(i == numOfVertices - 1 ) break;  
			
			/*for(NodeID j = descendants.size() - 1; j >= 0; --j){
				if(c_max_hop == root_hop[descendants[j]])
					if(usd[descendants[j]] == false)
						max_hop_set.push_back(descendants[j]);
			}*/
			
		//	NodeID second_source = max_hop_set[rand()%max_hop_set.size()];
			
		//	build_tmp_tree(second_source, second_parent_tree, second_coverage, second_descendants, second_root_hop, second_last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			calcover(descendants, parent_tree, coverage, root_hop, last_alive, depth);
			
			
			
		//	calcover_secondtree(second_descendants, second_parent_tree, parent_tree, second_coverage, second_root_hop, second_last_alive, second_depth);
			 
		//	NodeID topcan = candidate_gen(candidate, deg_inv, usd); 
			
	/*		vector<NodeID> level_howmany(50, 0);
			vector<NodeID> level_counts(50, 0);
		
			for(NodeID j = descendants.size() - 1; j >= 0; --j){
				NodeID tv = descendants[j];
				level_howmany[root_hop[tv]]++;
			}
		*/	
			int iterround = 30;
			//for(int j = 0; j < candidate.size(); ++j){
			//	NodeID tv = candidate[j];
			/*
			calupwardFanout(descendants, upwardcoverage, coverage, parent_tree,root_hop, graph);
			
			for(int j = descendants.size() - 1; j >= 0; --j){
				NodeID tv = descendants[j];
				if(root_hop[tv] == 0 ) continue;
				if(level_counts[root_hop[tv]] >= iterround) continue;
				int rollnum = rand()%level_howmany[root_hop[tv]];
				if(rollnum < iterround){
					//NodeID walk = bfs_walk(tv,  parent_tree, coverage,  root_hop, last_hop, que,que2,vis,  dst_r, usd,  graph);
					//NodeID walk = bfs_one_hop_walk(tv,  parent_tree, coverage,  root_hop, last_hop, que,que2,vis,  dst_r, usd,  graph);
					//NodeID walk = calupwardTree(tv, parent_tree, root_hop, graph);	
					NodeID walk = upwardcoverage[tv];
					level_counts[root_hop[tv]]++; 
					cout <<  tv << " " << root_hop[tv] << " " << coverage[tv] << " " << walk << " " << graph.vertices[tv + 1] - graph.vertices[tv] << " ";
				}
			}  
			*/
			vector<NodeID> upwardsPow;
			vector<NodeID> sigpath = findSigPath(chosen, coverage, parent_tree, usd, upwardsPow, graph);
		
		
		if(sigpath.size()==0){
			chosen = topcan(deg_inv, usd, last_available);
		}
			
		NodeID maxDegree = -1;
		NodeID chosenhop = -1;
	
		for(NodeID j = 0; j < sigpath.size(); ++j){
			NodeID curDegree = graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]];
			if(curDegree > maxDegree){
				maxDegree = curDegree;
				chosen = sigpath[j];
				chosenhop = j;
			}
			//cout << bfs_walk(sigpath[j], parent_tree, coverage,  root_hop, last_hop, que,que2,vis, dst_r, usd, graph) + coverage[sigpath[j]] << "," << curDegree << "," << coverage[sigpath[j]] << " " ;
		} 
		//cout << endl;
		
		/*
		if(sigpath.size()>2){
			NodeID mintrend = INF_WEIGHT;
			NodeID maxtrend = -1;
			NodeID mindegree = INF_WEIGHT;
			NodeID maxdegree = -1;
			NodeID chosentrend = 0;
			//cout << "Coverage Trends:";
			for(NodeID j = 0; j < sigpath.size() - 1; ++j){
				NodeID curtrend = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				NodeID curdegree = (graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]]);
				if(maxtrend <  curtrend){
					 maxtrend = curtrend;
					//chosentrend = j;
					//chosen = sigpath[j];
					//chosenhop = j;
				}
				if(mintrend > curtrend){
					mintrend = curtrend;
				}
				if(maxdegree < curdegree) maxdegree = curdegree;
				if(mindegree > curdegree) mindegree = curdegree;
			} 
			
			double max_score = -1;
			//cout << "Coverage Trends:";
			for(NodeID j = 0; j < sigpath.size() - 1; ++j){
				NodeID curtrend = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				NodeID curdegree = (graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]]);
				double curscore =  normalize(mindegree, maxdegree, mintrend, maxtrend, curdegree, curtrend);
				if(max_score < curscore){
					curscore = max_score;
					chosen = sigpath[j];
					chosenhop = j;
				}
			} 
		}*/
		
	/*	if(sigpath.size()<=1){
			chosen = topcan(deg_inv, usd);
		}
					
	
		if(sigpath.size()>1){
			NodeID maxtrend = -1;
			NodeID chosentrend = 0;
			cout << "Coverage Trends:";
			for(NodeID j = 1; j < sigpath.size() - 1; ++j){
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				//cout << " " << curDegree;
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				cout << " " <<  j << ":" << abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				if(maxtrend <  abs(coverage[sigpath[j+1]] - coverage[sigpath[j]])){
					 maxtrend = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
					chosentrend = j;
				}
			}
		}
		cout << endl;
		*/
	/*	for(NodeID j = 0; j < numOfVertices; ++j){
			if(usd[deg_inv[j]] == false ){
				if(root_hop[deg_inv[j]] == chosenhop){
					chosen = deg_inv[j];
				}
			}
		}*/
		
	//		cout << "iteration:" << i << " - " << currentsum << " - " << usd[chosen] << "-" << chosenhop << " - " << sigpath.size() << endl;
			 
	/*	if(usd[chosen] == true){
			cout << inv[i] << ", shit," << chosen << endl;
			for(NodeID j = 1; j < sigpath.size(); ++j){
				cout << " " << sigpath[j];
			}
			cout << endl;
		} 
		*/
			
	/*		
			for(NodeID j = 1; j < hop_counts.size(); ++j){
				if(hop_counts[j] != 0 )
					if(j < sigpath.size())
						cout <<  sigpath[j] << " " << j << " " << 0 << " " << hop_walks[j] / hop_counts[j] << " " << graph.vertices[sigpath[j] + 1] - graph.vertices[sigpath[j]] << " ";
					else
						cout <<  j << " " << j << " " << 0 << " " << hop_walks[j] / hop_counts[j] << " " << j << " ";

			}
			
			cout << endl; 
			*/
			//clear_upwardtmp(descendants, upwardcoverage);
			
			/*cout << "candidate" << endl;
			for(int j = 0; j < candidate.size(); ++j) 
				cout << candidate[j] << " ";
			cout << endl;*/
			
			 
			//chosen = simplepick(graph, usd, i, candidate,  coverage, root_hop);
		//	cout << "5:" << chosen << "," << graph.vertices[chosen+1] - graph.vertices[chosen] << "," << graph.vertices[candidate[0]+1] - graph.vertices[candidate[0]] << endl;
		//	cout << i << "," << last_available << "," << candidate.size() << endl;
			//chosen = pick(graph, usd, i, candidate, coverage, acc_count, parent_tree, root_hop, last_hop, last_alive, depth);
	//		chosen = pick(graph, usd, i, candidate,  coverage, second_coverage);
			
		/*	
			if(usd[border[i+1]] == false)
				if(acc_count[chosen] - 1 > 0 && acc_count[border[i+1]] - 1 > 0 )
				cout << root_hop[chosen] << "," << coverage[chosen] << "," << (acc[chosen] - coverage[chosen])/(acc_count[chosen] - 1) << "\t" << root_hop[border[i+1]] << "," << coverage[border[i+1]] << "," << (acc[border[i+1]] - coverage[border[i+1]])/(acc_count[border[i+1]]-1) << endl;
			
			chosen = border[i+1];
			*/
			
		//	last_iter_cover = coverage[chosen];
		//	last_iter_hop = root_hop[chosen];
			//if( i == taken )
			//	cout << i << "," << root_hop[chosen] << "," << last_hop[chosen] << endl;
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		//	clear_tmp(second_descendants, second_coverage, second_parent_tree, second_root_hop, second_depth);
		}
		
	/*	for(unordered_map<int, double>::iterator it = over_stats.begin(); it != over_stats.end(); ++it){
			cout << (*it).first << ":" << (double)((*it).second / actual_stats[(*it).first]) << endl;
		}
		*/
		vector<NodeID> remaining;
		remaining.reserve(numOfVertices - taken);
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			if(usd[deg_inv[i]] == false)
				remaining.push_back(deg_inv[i]);
		}
		
		for(NodeID i = taken; i < numOfVertices; ++i){
			chosen = remaining[i-taken];
		if(usd[chosen] == true) cout << "shit" << endl;
			inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, tmp_idx_parents, parents, graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		
		plabels.index_.resize(numOfVertices);
		
		double amount = 0; 
		for (NodeID v = 0; v < numOfVertices; ++v)  {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			plabels.index_[v].spt_v.resize(k);
			plabels.index_[v].spt_d.resize(k);
			plabels.index_[v].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_p[i] = tmp_idx_parents[v].second[i];
			
			
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
 
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
			
			tmp_idx_parents[v].first.clear();
			tmp_idx_parents[v].second.clear();
			tmp_idx_parents[v].first.shrink_to_fit();
			tmp_idx_parents[v].second.shrink_to_fit();
		}
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	}
	

	void undirected_weighted_sigpath(WGraph& wgraph, double &avg_degree, double &max_degree){
	    cout << "Building Coverage Ordering Based Labels" << endl;
	    if(wgraph.vertices.size() == 0 ){
            cout << "Empty graph "<< endl;
	    }
        if(wgraph.vertices.size()==1){
            numOfVertices = 1;
        }
	    inv.resize(numOfVertices);
	    rank.resize(numOfVertices);

	    vector<NodeID> deg_inv;
	    vector<NodeID> deg_rank;
	    deg_inv.resize(numOfVertices);
	    deg_rank.resize(numOfVertices);
	    double num_degree = 0;
	    double temp_degree = 0;
	    NodeID max_degree_node;
	    //temp variable for each node and degree(the random float is probably to avoid same degree)
	    vector<pair<float, NodeID> > deg(numOfVertices);
	    for (size_t v = 0; v < numOfVertices; ++v) {
	        deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) + float(random()) / RAND_MAX, v);
	        num_degree += wgraph.vertices[v + 1] - wgraph.vertices[v];
	        temp_degree = wgraph.vertices[v + 1] - wgraph.vertices[v];
	        if(temp_degree > max_degree){
	            max_degree = temp_degree;
	            max_degree_node = v;
	        }
	    }
	    //swap(deg[max_degree_node].second,deg[0].second);
	    avg_degree = num_degree/numOfVertices;
	    cout << "Number of Degree: " << avg_degree << endl;
	    cout << "Number of Vertices: " << numOfVertices << endl;
	    cout << "Number of Max Degree: " << max_degree << endl;

	    //sort the nodes into degree from largest
	    sort(deg.rbegin(), deg.rend());



	    for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
	    last_available = 0;
		
	    vector<NodeID> parent_tree(numOfVertices, numOfVertices);
	    vector<NodeID> root_hop(numOfVertices, 0);
	    vector<NodeID> coverage(numOfVertices, 0);
	    vector<NodeID> depth(numOfVertices, 0);
	    vector<NodeID> last_alive(numOfVertices, 0);
	    vector<NodeID> last_hop(numOfVertices, 0);
	    vector<NodeID> descendants;
	    descendants.reserve(numOfVertices);
		
	    vector<NodeID> candidate;
	    candidate.reserve(numOfVertices);
				
	    vector<NodeID> acc_count;
	    acc_count.resize(numOfVertices);
		
	    vector<bool> vis(numOfVertices);
	    queue<NodeID> visited_que;
	    vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
	    benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
	    vector<bool> usd(numOfVertices, false);

	    // Preparing basic structure for pl algorithm.
	    vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
	    vector<pair<vector<NodeID>, vector<EdgeWeight> > >
	    tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
                                         vector<EdgeWeight>(1, INF_WEIGHT)));

		
	    vector<pair<vector<NodeID>, vector<NodeID> > >
	    tmp_idx_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
                                                 vector<NodeID>(1, numOfVertices)));
				
				
	    vector<NodeID> parents(numOfVertices, numOfVertices);
		
	    NodeID chosen = deg_inv[0];
		 
	    NodeID taken = numOfVertices; // /5;
	    NodeID last_iter_cover = 0;
	    NodeID last_iter_hop = 0;
	    unordered_map<int, double> actual_stats;
	    unordered_map<int, double> over_stats;
		
	    vector<NodeID> upwardcoverage(numOfVertices, 0);
		
	    long long currentsum = 0;
		
	    for(NodeID i = 0; i < taken; ++i){
	        inv[i] = chosen;
	        rank[chosen] = i;
	        if(usd[chosen] == true) cout << "shit at " << chosen <<  endl;

	        //find the distances of the current source node to all connected nodes, cover is the amount of nodes it visited
	        NodeID actual_cover = labeling_source_dij(chosen, parent_tree, coverage, descendants, root_hop, last_hop, pqueue, visited_que, distances, vis, dst_r, usd, i, tmp_idx,tmp_idx_parents ,parents, wgraph);
	        //	NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
	        currentsum += actual_cover;
		
	        if(i == numOfVertices - 1 ) break;
			
	        calcover(descendants, parent_tree, coverage, root_hop, last_alive, depth);
			
	        vector<NodeID> upwardsPow;
	        NodeID max_degree = -1;
	        vector<NodeID> sigpath = findSigPath(chosen, coverage, parent_tree, usd,upwardsPow, wgraph, max_degree);


	        if(sigpath.size()<=2){
	            chosen = topcan(deg_inv, usd, last_available);
	        }


	        if(sigpath.size() > 2){
	            NodeID maxDegree = -1;
	            NodeID chosenhop = -1;
	            if( (double)max_degree / (double)sigpath.size() < 2 && sigpath.size() > 3){
	                NodeID max_degree_v = -1;
	                max_degree = -1;
	                //NodeID min_degree = INF_WEIGHT;
	                if(sigpath.size()>2){
	                    NodeID maxtrend = -1;
	                    NodeID chosentrend = 0;
	                    //cout << "Coverage Trends:";
	                    for(NodeID j = 0; j < sigpath.size() - 1; ++j){
	                        //NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
	                        //cout << " " << curDegree;
	                        //NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
	                        //cout << " " <<  j << ":" << abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
	                        NodeID cd = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
	                        if(max_degree < cd){
	                            max_degree = cd;
	                            max_degree_v = sigpath[j];
	                        }
	                        //if(min_degree > cd)
	                        //	min_degree = cd;
	                        if(maxtrend <  cd * abs(coverage[sigpath[j+1]] - coverage[sigpath[j]])){
	                            maxtrend = cd * abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
	                            chosentrend = j;
	                            chosen = sigpath[j];
	                        }
	                    }
	                }
	            }else{
	                for(NodeID j = 0; j < sigpath.size(); ++j){
	                    NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
	                    if(curDegree > maxDegree){
	                        maxDegree = curDegree;
	                        chosen = sigpath[j];
	                        chosenhop = j;
	                    }

	                }
	            }
			
		
		
	            //	if(sigpath.size() < 10 && sigpath.size() >2 &&  max_degree / min_degree > 10)
	            //	chosen = max_degree_v;
	            //cout << endl;
	            //cout << chosentrend << endl;
		 
	            //cout << "iteration:" << i << " - " << currentsum << " - " << usd[chosen] << "-" << chosenhop << " - " << sigpath.size() << endl;
			
	            clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
	        }
	    }
		
	        vector<NodeID> remaining;
	        remaining.reserve(numOfVertices - taken);
		
		
	        for(NodeID i = 0; i < numOfVertices; ++i){
	            if(usd[deg_inv[i]] == false)
	                remaining.push_back(deg_inv[i]);
	        }
		
	        for(NodeID i = taken; i < numOfVertices; ++i){
	            chosen = remaining[i-taken];
	            if(usd[chosen] == true) cout << "shit" << endl;
	            inv[i] = chosen;
	            rank[chosen] = i;
	            labeling_source_dij(chosen, parent_tree, coverage, descendants, root_hop, last_hop, pqueue, visited_que, distances, vis, dst_r, usd, i, tmp_idx,tmp_idx_parents, parents, wgraph);
	            clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
	        }
		
		
	        plabels.index_.resize(numOfVertices);
		
	        double amount = 0;
	        for (NodeID v = 0; v < numOfVertices; ++v)  {
	            NodeID k = tmp_idx[v].first.size();
			
	            amount += k ;
			
	            plabels.index_[v].spt_v.resize(k);
	            plabels.index_[v].spt_d.resize(k);
	            plabels.index_[v].spt_p.resize(k);
	            for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
	            for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
	            for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_p[i] = tmp_idx_parents[v].second[i];
			
			
	            tmp_idx[v].first.clear();
	            tmp_idx[v].second.clear();

	            tmp_idx[v].first.shrink_to_fit();
	            tmp_idx[v].second.shrink_to_fit();
			
	            tmp_idx_parents[v].first.clear();
	            tmp_idx_parents[v].second.clear();
	            tmp_idx_parents[v].first.shrink_to_fit();
	            tmp_idx_parents[v].second.shrink_to_fit();
	        }
	        cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	    }
	



	Coverage_Ordering_Path(WGraph& wgraph, double &average_degree, double&max_degree){
		undirected_weighted_sigpath(wgraph, average_degree, max_degree);
	}


	Coverage_Ordering_Path(Graph& graph, vector<NodeID> border){
		walk_stats(graph, border);
	}
	

	};

#endif
