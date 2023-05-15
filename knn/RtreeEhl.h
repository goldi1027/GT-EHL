//
// Created by Jinchun Du on 12/12/2022.
//

#ifndef RTREE_EHL
#define RTREE_EHL

#include "RStarTree.h"
#include "RStarTreeUtil.h"
#include "ebhl_v2.h"
#include "ebhl_query_v2.h"
#include "unordered_map"


using namespace std;
namespace rs = rstar;

namespace polyanya{
    class RtreeEhl{
        typedef EBHL_V2* EBHL_Ptr;
        typedef EBHL_query_v2* EBHL_QPtr;
    private:
        int K;
        Point start;
        vector<Point> goals;
        std::vector<std::unordered_map<std::string,int>> object_keyword;
        std::vector<std::vector<std::string>> query_keyword;
        int query_id;
        const EBHL_QPtr ehl_query;

        inline void initRtree() {
            assert(rte == nullptr);
            rte = new rs::RStarTree();
            rtEntries.clear();
            gids.clear();
            for (int i=0; i<(int)goals.size(); i++) gids.push_back(i);
            for (int& i: gids) {
                const Point& it = goals[i];
                rs::Mbr mbr(it.x, it.x, it.y, it.y);
                rs::LeafNodeEntry leaf(mbr, (rs::Data_P)(&i));
                rtEntries.push_back(leaf);
            }

            for (auto& it: rtEntries)
                rte->insertData(&it);
        }

    public:
        double search_cost;
        double rtree_cost;
        rs::RStarTree* rte = nullptr;
        std::vector<rs::LeafNodeEntry> rtEntries;
        std::vector<int> gids;
        struct object{
            int id;
            Point previous_location;
            Point new_location;
        };
        // object location at each timestamp
        vector<vector<object>> object_timestamp;

        RtreeEhl(EBHL_QPtr ebq):ehl_query(ebq){};
        ~RtreeEhl(){
            delete ehl_query;
        }

        inline void set_start(Point s) {
            start = s;
            rtree_cost = 0;
        }

        inline void set_K(int k) { this->K = k; }
        inline void set_goals(std::vector<Point> gs) {
            goals = std::vector<Point>(gs);
            initRtree();
        }

        inline void set_keyword(std::vector<std::vector<std::string>>q_keyword, std::vector<std::vector<std::string>>o_keyword, int start_id){
            query_id = start_id;
            query_keyword = q_keyword;
            object_keyword.resize(o_keyword.size());
            for(int i = 0; i < o_keyword.size(); ++i){
                for(int j = 0; j < o_keyword[i].size(); ++j){
                    object_keyword[i][o_keyword[i][j]] = 1;
                }
            }
            object_keyword.shrink_to_fit();
        }

        inline void load_object_timestamp(int timestamp, int object_num, istream &infile) {
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

        inline vector<double> search(){
            priority_queue<double, vector<double>> maxh;
            vector<pair<double, Point>> edists;
            vector<double> odists;
            rs::MinHeap heap;
            rs::Point P(start.x, start.y);
            double D = sqrt(rs::RStarTreeUtil::minDis2(P, rte->root->mbrn));
            heap.push(rs::MinHeapEntry(D, rte->root));
            rs::MinHeapEntry cur(INF, (rs::Entry_P)nullptr);
            while (true) {
                auto stime = std::chrono::steady_clock::now();
                cur = rs::RStarTreeUtil::iNearestNeighbour(heap, P);
                auto etime = std::chrono::steady_clock::now();
                rtree_cost += std::chrono::duration_cast<std::chrono::microseconds>(etime- stime).count();
                if (cur.key == INF) // not found
                    break;
                int gid = *((int*)cur.entryPtr->data);
                double de = goals[gid].distance(start);
                if ((int)maxh.size() == K && maxh.top() <= de)
                    break;
                double dist = get_ehl_distance(start,goals[gid]);
                if (dist == INF) continue;
                if ((int)maxh.size() < K) {
                    maxh.push(dist);
                } else if (dist < maxh.top()) {
                    maxh.pop();
                    maxh.push(dist);
                }
            }
            while (!maxh.empty()) {
                odists.push_back(maxh.top());
                maxh.pop();
            }
            sort(odists.begin(), odists.end());
            return odists;
        }

       inline vector<double> keyword_search() {
            priority_queue<double, vector<double>> maxh;
            vector<pair<double, Point>> edists;
            vector<double> odists;
            rs::MinHeap heap;
            rs::Point P(start.x, start.y);
            double D = sqrt(rs::RStarTreeUtil::minDis2(P, rte->root->mbrn));
            heap.push(rs::MinHeapEntry(D, rte->root));
            rs::MinHeapEntry cur(INF, (rs::Entry_P)nullptr);
            while (true) {
                auto stime = std::chrono::steady_clock::now();
                cur = rs::RStarTreeUtil::iNearestNeighbour(heap, P);
                auto etime = std::chrono::steady_clock::now();
                rtree_cost += std::chrono::duration_cast<std::chrono::microseconds>(etime- stime).count();
                if (cur.key == INF) // not found
                    break;
                int gid = *((int*)cur.entryPtr->data);
                double de = goals[gid].distance(start);
                if ((int)maxh.size() == K && maxh.top() <= de)
                    break;
                if(object_keyword[gid][query_keyword[query_id][0]] > 0 && object_keyword[gid][query_keyword[query_id][1]] > 0){
                    double dist = get_ehl_distance(start,goals[gid]);
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

        inline void object_update(int timestamp, int object_id, int&updated){
            if(object_timestamp[timestamp][object_id].previous_location != object_timestamp[timestamp][object_id].new_location){
                updated++;
                rs::LeafNodeEntry& it = rtEntries[object_id];
                rs::Mbr cur(object_timestamp[timestamp][object_id].previous_location.x, object_timestamp[timestamp][object_id].previous_location.y,
                            object_timestamp[timestamp][object_id].previous_location.x, object_timestamp[timestamp][object_id].previous_location.y);
                rs::Mbr after(object_timestamp[timestamp][object_id].new_location.x, object_timestamp[timestamp][object_id].new_location.y,
                              object_timestamp[timestamp][object_id].new_location.x, object_timestamp[timestamp][object_id].new_location.y);
                rte->update(&it, cur,after);
            }
        }


    };
}

#endif //RTREE_EHL
