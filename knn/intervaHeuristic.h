#pragma once
#include "searchnode.h"
#include "searchinstance.h"
#include "successor.h"
#include "mesh.h"
#include "point.h"
#include "cpool.h"
#include "timer.h"
#include "expansion.h"
#include <queue>
#include <vector>
#include <ctime>
#include <unordered_map>

namespace polyanya {

class IntervalHeuristic{
    typedef std::priority_queue<SearchNodePtr, std::vector<SearchNodePtr>,
                                PointerComp<SearchNode> > pq;
    private:
        int K = 1;
        warthog::mem::cpool* node_pool;
        MeshPtr mesh;
        Point start;
        int start_id;
        std::vector<std::vector<std::string>> query_keyword;
        std::vector<Point> goals;
        std::vector<std::unordered_map<std::string,int>> object_keyword;

        // kNN has k final node
        std::vector<SearchNodePtr> final_nodes;
        // poly_id: goal1, goal2, ...
        std::vector<std::vector<int>> end_polygons;
        // gid: polyid
        std::vector<int> goal_polyid; 
        // <i, v>: reached ith goal with cost v
        std::map<int, double> reached;
        pq open_list;

        // Best g value for a specific vertex.
        std::vector<double> root_g_values;
        // Contains the current search id if the root has been reached by
        // the search.
        std::vector<int> root_search_ids;  // also used for root-level pruning
        int search_id;

        warthog::timer timer;

        // Pre-initialised variables to use in search().
        Successor* search_successors;
        SearchNode* search_nodes_to_push;
				bool isZero = false;

        void init() {
            verbose = false;
            search_successors = new Successor [mesh->max_poly_sides + 2];
            search_nodes_to_push = new SearchNode [mesh->max_poly_sides + 2];
            node_pool = new warthog::mem::cpool(sizeof(SearchNode));
            init_root_pruning();
        }

        void init_root_pruning() {
            assert(mesh != nullptr);
            search_id = 0;
            size_t num_vertices = mesh->mesh_vertices.size();
            root_g_values.resize(num_vertices);
            root_search_ids.resize(num_vertices);
            fill(root_search_ids.begin(), root_search_ids.end(), 0);
        }

        void init_search() {
            assert(node_pool);
            node_pool->reclaim();
            search_id++;
            open_list = pq();
            final_nodes = std::vector<SearchNodePtr>();
            reached.clear();
            nodes_generated = 0;
            nodes_pushed = 0;
            nodes_popped = 0;
            nodes_pruned_post_pop = 0;
            successor_calls = 0;
            gen_initial_nodes();
        }
        void set_end_polygon();
        void gen_initial_nodes();
        int succ_to_node(
            SearchNodePtr parent, Successor* successors,
            int num_succ, SearchNodePtr nodes
        );

        void print_node(SearchNodePtr node, std::ostream& outfile);

    public:
        int nodes_generated;        // Nodes stored in memory
        int nodes_pushed;           // Nodes pushed onto open
        int nodes_popped;           // Nodes popped off open
        int nodes_pruned_post_pop;  // Nodes we prune right after popping off
        int successor_calls;        // Times we call get_successors
        int query_counter;
        bool verbose;

        IntervalHeuristic() { }
        IntervalHeuristic(MeshPtr m) : mesh(m) { init(); }
        IntervalHeuristic(int k, MeshPtr m, Point s, std::vector<Point> gs) :
            K(k), mesh(m), start(s), goals(gs) { init(); }
        IntervalHeuristic(IntervalHeuristic const &) = delete;
        void operator=(IntervalHeuristic const &x) = delete;
        ~IntervalHeuristic() {
            if (node_pool) {
                delete node_pool;
            }
            delete[] search_successors;
            delete[] search_nodes_to_push;
        }

        void set_K(int k) { this->K = k; }

        void set_start(Point s, int id) {
          start = s;
          start_id = id;
        }

        void set_goals(std::vector<Point> gs) {
            goals.clear();
            for (const auto it: gs)
              goals.push_back(it);
            final_nodes.clear();
            set_end_polygon();
        }

        void set_goals(std::vector<Point> gs, std::vector<Point> new_gs) {
            goals.clear();
            for (const auto it: gs)
                goals.push_back(it);
            final_nodes.clear();
            set_end_polygon();
            for(const auto it:new_gs){
                goals.push_back(it);
                goal_polyid.push_back(-1);
            }
        }


    inline void insert_end_polygon(int id){
        int poly_id = get_point_location_in_search(goals[id], mesh, verbose).poly1;
        end_polygons[poly_id].push_back(id);
        goal_polyid[id] = poly_id;
    }
    inline void delete_endpolys(int gid) {
            int curpid = goal_polyid[gid];
        // remove gid in current poly
            if (curpid != -1) {
            for (int i=0; i<(int)end_polygons[curpid].size(); i++) {
                if (end_polygons[curpid][i] == gid) {
                    int tmp = end_polygons[curpid].back();
                    end_polygons[curpid][end_polygons[curpid].size()-1] = gid;
                    end_polygons[curpid][i] = tmp;
                    end_polygons[curpid].pop_back();
                    break;
                }
            }
        }
    }

        inline void maintain_endpolys(int gid, int curpid, int newpid) {
          goal_polyid[gid] = newpid;
          // remove gid in current poly
          if (curpid != -1) {
            for (int i=0; i<(int)end_polygons[curpid].size(); i++) {
              if (end_polygons[curpid][i] == gid) {
                int tmp = end_polygons[curpid].back();
                end_polygons[curpid][end_polygons[curpid].size()-1] = gid;
                end_polygons[curpid][i] = tmp;
                end_polygons[curpid].pop_back();
                break;
              }
            }
          }
          else {
            std::cerr << "gid: " << gid << " poly id is -1" << std::endl; 
          }
          // add gid to new poly
          if (newpid != -1) {
            end_polygons[newpid].push_back(gid);
          }
          else {
            std::cerr << "gid: " << gid << " new poly id is -1" << std::endl; 
          }
        }

        inline void update_goals(std::vector<Point> gs, int& update_cnt) {
          update_cnt = 0;
          for (int i=0; i<(int)goals.size(); i++) 
          if (goals[i] != gs[i]) {
            // update poly
            goals[i] = gs[i];
            int pid = goal_polyid[i];
            const PolyContainment res = mesh->poly_contains_point(pid, gs[i]);
            if (res.type == PolyContainment::OUTSIDE) {
              update_cnt++;
              int newpid = get_point_location_in_search(gs[i], mesh, verbose).poly1;
              maintain_endpolys(i, pid, newpid);
            } 
            else {
              continue;
            }
          }
          final_nodes.clear();
        }

        void set_keyword(std::vector<std::vector<std::string>>q_keyword, std::vector<std::vector<std::string>>o_keyword){
            query_keyword = q_keyword;
            object_keyword.resize(o_keyword.size());
            for(int i = 0; i < o_keyword.size(); ++i){
                for(int j = 0; j < o_keyword[i].size(); ++j){
                    object_keyword[i][o_keyword[i][j]] = 1;
                }
            }
            object_keyword.shrink_to_fit();
        }
        int search();

        int keyword_search();

        double get_cost(int k) {
          if (k >= (int)final_nodes.size()) {
            return -1;
          }
          return final_nodes[k]->f;
        }

        double get_search_micro()
        {
            return timer.elapsed_time_micro();
        }

        void get_path_points(std::vector<Point>& out, int k);
        void get_goals(std::vector<Point> &out, int k);
        void print_search_nodes(std::ostream& outfile, int k);
        void deal_final_node(const SearchNodePtr node);
        void gen_final_nodes(const SearchNodePtr node, const Point& rootPoint);
        void gen_keyword_final_nodes(const SearchNodePtr node, const Point& rootPoint);
				void setZero(bool flag) {
					this->isZero = flag;
				}

        double get_gid(int k) {
          if (k >= (int)final_nodes.size()) return -1;
          else return final_nodes[k]->goal_id;
        }

        int get_goal_ord(int gid) {
          for (int i=0; i<K; i++) if (final_nodes[i]->goal_id == gid) return i;
          return -1;
        }


};

}
