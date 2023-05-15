//
// Created by goldi on 9/01/23.
//

#ifndef EBHL_IRTREEEHL_H
#define EBHL_IRTREEEHL_H
#include "IRStarTree.h"
#include "IRStarTreeUtil.h"
#include "ebhl.h"
#include "ebhl_query_v2.h"
#include "unordered_map"


using namespace std;
namespace irs = irstar;

namespace polyanya{
    class IRTreeEhl{
        typedef EBHL* EBHL_Ptr;
        typedef EBHL_query_v2* EBHL_QPtr;
    private:
        int K;
        Point start;
        vector<Point> goals;
        std::vector<std::unordered_map<std::string,int>> object_keyword;
        std::vector<std::vector<std::string>> query_keyword;
        int query_id;
        const EBHL_QPtr ehl_query;

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

            for (auto& it: rtEntries)
                irte->insertData(&it);
        }

    public:
        double search_cost;
        double rtree_cost;
        int num_in_heap;
        irs::RStarTree* irte = nullptr;
        std::vector<irs::LeafNodeEntry> rtEntries;
        std::vector<int> gids;
        int query_counter;
        struct object{
            int id;
            Point previous_location;
            Point new_location;
        };
        // object location at each timestamp
        vector<vector<object>> object_timestamp;

        IRTreeEhl(EBHL_QPtr ebq):ehl_query(ebq){
            irte = nullptr;
        };
        ~IRTreeEhl(){
            delete ehl_query;
            if(irte){
                delete irte;
            }
        }

        inline void set_start(Point s, int id) {
            start = s;
            query_id = id;
        }

        inline void set_K(int k) { this->K = k; }
        inline void set_goals(std::vector<Point> gs) {
            goals = std::vector<Point>(gs);
            initRtree();
        }

        void initRtree(const std::vector<Point>& gs) {
            assert(irte == nullptr);
            irte = new irs::RStarTree();
            for (int i=0; i<gs.size(); i++) {
                auto& it = rtEntries[i];
                // cout << "insertData: "; it.print(".");
                irte->insertData(&it);
            }
        }

        inline void set_goals(const std::vector<Point>& gs,
                              const std::vector<Point>& newg) {
            goals.clear();
            for (auto it: gs) goals.push_back(it);
            for (auto it: newg) goals.push_back(it);
            gids.resize(goals.size());
            for (int i=0; i<goals.size(); i++) gids[i] = i;

            rtEntries.clear();
            for (int i=0; i<goals.size(); i++) {
                const Point& it = goals[i];
                irs::Mbr mbr(it.x, it.x, it.y, it.y);
                irs::LeafNodeEntry leaf(mbr, (irs::Data_P)(&gids[i]), object_keyword[i]);
                rtEntries.push_back(leaf);
            }
            initRtree(gs);
        }
        //add in additional goals to be inserted
        inline void insert_goals(Point gs, int index) {
            goals.push_back(gs);
            gids.push_back(index);
        }

        inline void set_keyword(std::vector<std::vector<std::string>>q_keyword, std::vector<std::vector<std::string>>o_keyword){
            query_keyword = q_keyword;
            object_keyword.resize(o_keyword.size());
            for(int i = 0; i < o_keyword.size(); ++i){
                for(int j = 0; j < o_keyword[i].size(); ++j){
                    object_keyword[i][o_keyword[i][j]] = 1;
                }
            }
            object_keyword.shrink_to_fit();
        }

        inline void add_keyword(std::vector<std::vector<std::string>>o_keyword){
            int size = object_keyword.size();
            object_keyword.resize(o_keyword.size());
            for(int i = size; i < o_keyword.size(); ++i){
                for(int j = 0; j < o_keyword[i].size(); ++j){
                    object_keyword[i][o_keyword[i][j]] = 1;
                }
            }
        }

        inline void delete_object(const int&object_id){
            auto& entry = rtEntries.at(object_id);
            irte->delete_entry(&entry);
        }
//        inline void insert(const int &object_id){
//            const Point& it = goals[object_id];
//            irs::Mbr mbr(it.x, it.x, it.y, it.y);
//            irs::LeafNodeEntry leaf(mbr, (irs::Data_P)(&gids[object_id]), object_keyword[object_id]);
//            rtEntries.push_back(leaf);
//            irte->insertData(&leaf);
//        }

        inline void insert(const int &object_id){
            const Point& it = goals[object_id];
            irs::Mbr mbr(it.x, it.x, it.y, it.y);
            auto& leaf = rtEntries[object_id];
            irte->insertData(&leaf);
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
        inline void load_object_timestamp(int timestamp, istream &infile) {
            if(!infile){
                cout << "Cannot find file!" << endl;
            }
            object_timestamp.resize(timestamp);
            int time;
            object o;
            while (infile >> o.id >> time >> o.previous_location.x >> o.previous_location.y >> o.new_location.x >> o.new_location.y)
            {
                object_timestamp[time].push_back(o);
            }
        }

        inline double get_ehl_distance(Point &q, Point &c){
            ehl_query->set_start_goal(q,c);
            ehl_query->search();
            double distance = ehl_query->get_cost();
            return distance;
        }

        inline vector<double> keyword_search() {
            num_in_heap = 0;
            search_cost = 0;
            rtree_cost = 0;
            warthog::timer rtreetime;
            warthog::timer ehltime;
            priority_queue<double, vector<double>> maxh;
            vector<pair<double, Point>> edists;
            vector<double> odists;
            irs::MinHeap heap;
            irs::Point P(start.x, start.y);
            double D = sqrt(irs::RStarTreeUtil::minDis2(P, irte->root->mbrn));
            heap.push(irs::MinHeapEntry(D, irte->root));
            irs::MinHeapEntry cur(INF, (irs::Entry_P)nullptr);
            while (true) {
                rtreetime.start();
                cur = irs::RStarTreeUtil::iNearestNeighbourKw(heap, P,query_keyword[query_id]);
                //cur = irs::RStarTreeUtil::iNearestNeighbour(heap, P);
                rtreetime.stop();
                rtree_cost += rtreetime.elapsed_time_micro();
                num_in_heap++;
                if (cur.key == INF) // not found
                    break;
                int gid = *((int*)cur.entryPtr->data);
                double de = goals[gid].distance(start);
                if ((int)maxh.size() == K && maxh.top() <= de)
                    break;
                query_counter= 0;
                for(int i = 0; i < query_keyword[query_id].size(); ++i){
                    if(object_keyword[gid].find(query_keyword[query_id][i])!=object_keyword[gid].end()){
                        query_counter++;
                    }
                }
                if(query_counter == query_keyword[query_id].size()) {
                    ehltime.start();
                    double dist = get_ehl_distance(start,goals[gid]);
                    ehltime.stop();
                    search_cost += ehltime.elapsed_time_micro();
                    if (dist == INF) continue;
                    if ((int)maxh.size() < K) {
                        maxh.push(dist);
                    } else if (dist < maxh.top()) {
                        maxh.pop();
                        maxh.push(dist);
                    }
                }
            }
            while (!maxh.empty()) {
                odists.push_back(maxh.top());
                maxh.pop();
            }
            sort(odists.begin(), odists.end());
            return odists;
        }



    };
}
#endif //EBHL_IRTREEEHL_H
