#pragma once
#include "searchinstance.h"
#include "point.h"
#include "timer.h"
// #include "RStarTree.h"
// #include "RStarTreeUtil.h"
#include "IRStarTree.h"
#include "IRStarTreeUtil.h"
#include <chrono>
#include <queue>
#include <vector>
#include <ctime>
#include <unordered_map>

using namespace std;

namespace polyanya {

namespace irs = irstar;

class IERPolyanya {
    typedef std::priority_queue<SearchNodePtr, std::vector<SearchNodePtr>,
                                PointerComp<SearchNode> > pq;
    private:
        int K = 1;
        Point start;
        std::vector<Point> goals;
        warthog::timer timer;
        SearchInstance* polyanya;
        std::vector<std::unordered_map<std::string,int>> object_keyword;
        std::vector<std::vector<std::string>> query_keyword;
        int query_id;

        void init() {
            verbose = false;
        }

        void initRtree() {
          assert(irte == nullptr);
          irte = new irs::RStarTree();
          rtEntries.clear();
          gids.clear();
          for (int i=0; i<(int)goals.size(); i++) gids.push_back(i);
          for (int& i: gids) {
            const Point& it = goals[i];
            irs::Mbr mbr(it.x, it.x, it.y, it.y);
            irs::LeafNodeEntry leaf(mbr, (irs::Data_P)(&i), object_keyword[i]);
            rtEntries.push_back(leaf);
          }
          for (auto& it: rtEntries){
              irte->insertData(&it);
          }
        }

        void init_search() {
          num_in_heap = 0;
          search_cost = 0;
          rtree_cost = 0;
          nodes_generated = 0;
          nodes_pushed = 0;
          nodes_popped = 0;
          tot_hit = 0;
        }

    public:
        int nodes_generated;        // Nodes stored in memory
        int nodes_pushed;           // Nodes pushed onto open
        int nodes_popped;           // Nodes popped off open
        int tot_hit;
        double search_cost;
        double rtree_cost;
        int num_in_heap;
        int query_counter;
        bool verbose;
        // rs::RStarTree* rte = nullptr;
        irs::RStarTree* irte = nullptr;
        std::vector<irs::LeafNodeEntry> rtEntries;
        std::vector<int> gids;

        struct object{
            int id;
            Point previous_location;
            Point new_location;
        };
        // object location at each timestamp
        vector<vector<object>> object_timestamp;

        IERPolyanya () { }
        IERPolyanya(SearchInstance* si): polyanya(si) {
          irte = nullptr; 
        };
        IERPolyanya(IERPolyanya const &) = delete;
        void operator=(IERPolyanya const &x) = delete;
        ~IERPolyanya() {
            if (irte)
              delete irte;
        }

        void set_K(int k) { this->K = k; }
        void set_goals(std::vector<Point> gs) {
          goals = std::vector<Point>(gs);
          initRtree();
        }
        void set_start(Point s, int id) {
          start = s;
          query_id = id;
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

        inline void update_goals(const std::vector<Point>& gs, int& update_cnt) {
          update_cnt = 0;
          for (size_t i=0; i<gs.size(); i++)
          if (goals[i] != gs[i]) {
            // update rtree
            update_cnt ++;
            goals[i] = gs[i];
            auto& entry = rtEntries.at(i);
            irs::Mbr cur(goals[i].x, goals[i].x, goals[i].y, goals[i].y);
            irs::Mbr nxt(gs[i].x, gs[i].x, gs[i].y, gs[i].y);
            irte->update(&entry, cur, nxt);
          }
        }

        vector<double> search();

        vector<double> keyword_search();

    void load_object_timestamp(int timestamp, istream &infile);

    void object_update(int timestamp, int object_id, int &updates);
};

}
