#include "IERPolyanya.h"
#include "geometry.h"
#include "vertex.h"
#include "point.h"
#include "consts.h"
#include <queue>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <ctime>

using namespace std;

namespace polyanya {
    namespace irs = irstar;
    vector<double> IERPolyanya::search() {
        init_search();
        priority_queue<double, vector<double>> maxh;
        vector<pair<double, Point>> edists;
        vector<double> odists;
        irs::MinHeap heap;
        irs::Point P(start.x, start.y);
        double D = sqrt(irs::RStarTreeUtil::minDis2(P, irte->root->mbrn));
        heap.push(irs::MinHeapEntry(D, irte->root));
        irs::MinHeapEntry cur(INF, (irs::Entry_P)nullptr);
        while (true) {
            auto stime = std::chrono::steady_clock::now();
            cur = irs::RStarTreeUtil::iNearestNeighbour(heap, P);
            auto etime = std::chrono::steady_clock::now();
            rtree_cost += std::chrono::duration_cast<std::chrono::microseconds>(etime- stime).count();
            if (cur.key == INF) // not found
                break;
            int gid = *((int*)cur.entryPtr->data);
            double de = goals[gid].distance(start);
            if ((int)maxh.size() == K && maxh.top() <= de)
                break;
            polyanya->set_start_goal(start, goals[gid]);
            bool found = polyanya->search();
            tot_hit++;
            search_cost += polyanya->get_search_micro();
            nodes_generated += polyanya->nodes_generated;
            if (!found) continue;
            double dist = polyanya->get_cost();
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

    void IERPolyanya::load_object_timestamp(int timestamp, istream& infile){
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

    void IERPolyanya::object_update(int timestamp,int object_id, int&updated){
        if(object_timestamp[timestamp][object_id].previous_location != object_timestamp[timestamp][object_id].new_location){
            updated++;
            irs::LeafNodeEntry& it = rtEntries[object_id];
            irs::Mbr cur(object_timestamp[timestamp][object_id].previous_location.x, object_timestamp[timestamp][object_id].previous_location.y,
                         object_timestamp[timestamp][object_id].previous_location.x, object_timestamp[timestamp][object_id].previous_location.y);
            irs::Mbr after(object_timestamp[timestamp][object_id].new_location.x, object_timestamp[timestamp][object_id].new_location.y,
                           object_timestamp[timestamp][object_id].new_location.x, object_timestamp[timestamp][object_id].new_location.y);
            irte->update(&it, cur,after);
        }
    }
    vector<double> IERPolyanya::keyword_search() {
        warthog::timer rtimer;
        init_search();
        priority_queue<double, vector<double>> maxh;
        vector<pair<double, Point>> edists;
        vector<double> odists;
        irs::MinHeap heap;
        irs::Point P(start.x, start.y);
        double D = sqrt(irs::RStarTreeUtil::minDis2(P, irte->root->mbrn));
        heap.push(irs::MinHeapEntry(D, irte->root));
        irs::MinHeapEntry cur(INF, (irs::Entry_P)nullptr);
        while (true) {
            //auto stime = std::chrono::steady_clock::now();
            rtimer.start();
            cur = irs::RStarTreeUtil::iNearestNeighbourKw(heap, P, query_keyword[query_id]);
            rtimer.stop();
            num_in_heap++;
//            auto etime = std::chrono::steady_clock::now();
//            rtree_cost += std::chrono::duration_cast<std::chrono::microseconds>(etime- stime).count();
            rtree_cost += rtimer.elapsed_time_micro();
            if (cur.key == INF) // not found
                break;
            int gid = *((int*)cur.entryPtr->data);
            double de = goals[gid].distance(start);
            if ((int)maxh.size() == K && maxh.top() <= de)
                break;
            query_counter= 0;
            for(int i = 0; i < query_keyword[query_id].size(); ++i){
                if(object_keyword[gid][query_keyword[query_id][i]] > 0){
                    query_counter++;
                }
//                if(object_keyword[gid].find(query_keyword[query_id][i])!=object_keyword[gid].end()){
//                    query_counter++;
//                }
            }
            if(query_counter == query_keyword[query_id].size()){
                polyanya->set_start_goal(start, goals[gid]);
                bool found = polyanya->search();
                tot_hit++;
                search_cost += polyanya->get_search_micro();
                nodes_generated += polyanya->nodes_generated;
                if (!found) continue;
                double dist = polyanya->get_cost();
                if ((int)maxh.size() < K) {
                    maxh.push(dist);
                } else if (dist < maxh.top()) {
                    maxh.pop();
                    maxh.push(dist);
                }
            }
//            if(object_keyword[gid][query_keyword[query_id][0]] > 0 && object_keyword[gid][query_keyword[query_id][1]] > 0){
//                polyanya->set_start_goal(start, goals[gid]);
//                bool found = polyanya->search();
//                tot_hit++;
//                search_cost += polyanya->get_search_micro();
//                nodes_generated += polyanya->nodes_generated;
//                if (!found) continue;
//                double dist = polyanya->get_cost();
//                if ((int)maxh.size() < K) {
//                    maxh.push(dist);
//                } else if (dist < maxh.top()) {
//                    maxh.pop();
//                    maxh.push(dist);
//                }
//            }
        }
        while (!maxh.empty()) {
            odists.push_back(maxh.top());
            maxh.pop();
        }
        sort(odists.begin(), odists.end());
        return odists;
    }
}
