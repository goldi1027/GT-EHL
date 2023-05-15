//
// Created by Jinchun Du on 12/9/2022.
//
/**
 * Grid tree implementation
 *
 */
#ifndef EBHL_EHL_TREE_H
#define EBHL_EHL_TREE_H
#include <vector>
#include <unordered_map>
#include "ebhl.h"
#include "ebhl_query_v2.h"

using namespace std;
namespace polyanya{
    typedef EBHL* EBHL_Ptr;
    typedef EBHL_query_v2* EBHL_QPtr;

    class EhlKnn{
    private:
        //initialise ehl for ehl knn
        const EBHL_QPtr ehl_query;
    public:
        struct cor{
            Point min_point;
            Point max_point;
        };
        //parent node's children indexes
        struct children{
            int top_right;
            int top_left;
            int bottom_right;
            int bottom_left;
        };

        //object with timestamp
        struct raw_object{
            int id;
            Point previous_location;
            Point new_location;
        };

        //object information
        struct object{
            int id;
            Point cordinate;
            int cellID;
            vector<string> keywords;
            unordered_map<string,int> object_keywords;
        };

        //query data structure
        struct query{
            int id;
            Point cordinate;
            double query_score;
            vector<string> query_keyword;
            vector<int> knn_object;
            vector<int> cell_visited;
        };
        //info for each level
        struct node{
           unsigned nodeID;
           cor cordinate;
           int parent;
           //level also used as cellID for object nodes
           int level;
           children child;
           double distance;
           bool has_object;
           unordered_map<string,int> keyword_list;
           vector<int> cell_object_list;
           vector<int> cell_query_list;
           void print(const node& p) const {
               printf("addr: %p id: %d\n", (void*)&p, p.nodeID);
           }
        };

        struct heapCompare
        {
            bool operator()(const node* x,
                            const node* y) const
            {
                return x->distance > y->distance;
            }
        };

        typedef node* nodePtr;

        typedef std::priority_queue<nodePtr, vector<nodePtr>,
                heapCompare> pq;

        struct objectCompare{
            bool operator() (const node *i, const node*j) {return(i->distance < j->distance);}
        };

        typedef std::priority_queue<nodePtr, vector<nodePtr>,
                objectCompare> object_pq;


        int height, width, grid_size, num_of_column,num_of_rows;
        int k_num;
        int query_counter;
        int timestamp_num;
        vector<node>ehl_tree;
        pq min_heap;
        object_pq max_heap;


        //object data structure
        vector<node> object_list;
        vector<vector<string>> object_keyword_list;

        //object timestamp ist keeps track of all the object's location through all timestamps
        vector<vector<raw_object>> object_timestamp_list;
        vector<int> knn_test;

        //insertion
        vector<int> added_objects;
        vector<int> free_objects;

        //query data structure
        vector<query> query_list;
        vector<node> objects_in_heap;

        EhlKnn(int mapHeight, int mapWidth, int size,  EBHL_QPtr ebq) : ehl_query(ebq) {initialise_tree(mapHeight, mapWidth, size);}


        ~EhlKnn(){
            ehl_tree.clear();
        }

        /**
         * Initialise Grid Tree and begins to subdivide from root node
         * @param mapHeight
         * @param mapWidth
         * @param size
         */

        void initialise_tree(int mapHeight, int mapWidth, int size){
            ehl_tree.clear();
            height = mapHeight;
            width = mapWidth;
            grid_size = size;
            num_of_rows = ceil((double) width / grid_size);
            num_of_column = ceil((double) height / grid_size);


            //initialise root node for tree
            node root_node;
            root_node.cordinate.min_point.x = 0;
            root_node.cordinate.min_point.y = 0;
            root_node.cordinate.max_point.x = width;
            root_node.cordinate.max_point.y = height;
            root_node.level = 0;
            root_node.parent = -1;
            root_node.nodeID = 0;
            root_node.distance = 0;
            root_node.has_object = true;
            ehl_tree.push_back(root_node);
            decompose_node(root_node);
        }



        /**
         * Load keyword information for all objects
         * Stores object keyword in an object list
         * @param infile
         */
        void load_keywords(istream& infile){
            if(!infile){
                cout << "Cannot load object file!" << endl;
            }
            int object_index;
            string keyword;
            object_list.resize(object_timestamp_list[0].size());
            object_keyword_list.resize(object_timestamp_list[0].size());
            while(infile >> object_index >> keyword){
                object_list[object_index].keyword_list[keyword]++;
                object_list[object_index].nodeID = object_index;
                object_list[object_index].cordinate.max_point = object_timestamp_list[0][object_index].previous_location;
                object_list[object_index].parent = -2;
                object_keyword_list[object_index].emplace_back(keyword);
            }
        }

        /**
         * Loading keywords for insertion experiment
         * @param infile
         */
        void load_insertion_keywords(istream& infile){
            if(!infile){
                cout << "Cannot load object file!" << endl;
            }
            for(int i = 0; i < object_list.size(); ++i){
                added_objects.push_back(object_list[i].nodeID);
            }
            int object_index;
            string keyword;
            int size = object_list.size();
            object_list.resize(object_timestamp_list[0].size());
            object_keyword_list.resize(object_timestamp_list[0].size());
            while(infile >> object_index >> keyword){
                object_list[object_index].keyword_list[keyword]++;
                object_list[object_index].nodeID = object_index;
                object_list[object_index].cordinate.max_point = object_timestamp_list[0][object_index].previous_location;
                object_list[object_index].parent = -2;
                object_keyword_list[object_index].emplace_back(keyword);
            }
            for(int i = size; i < object_list.size(); ++i){
                free_objects.push_back(object_list[i].nodeID);
            }
        }


        void load_object(int timestamp, istream& infile){
            if(!infile){
                cout << "Cannot load object file!" << endl;
            }
            object_timestamp_list.resize(timestamp);
            int time;
            raw_object o;
            while (infile >> o.id >> time >> o.previous_location.x >> o.previous_location.y >> o.new_location.x >> o.new_location.y)
            {
                object_timestamp_list[time].push_back(o);
            }
        }

        void load_insertion_object(istream& infile){
            if(!infile){
                cout << "Cannot load object file!" << endl;
            }
            int time;
            raw_object o;
            while (infile >> o.id >> time >> o.previous_location.x>> o.previous_location.y >> o.new_location.x >> o.new_location.y)
            {
                object_timestamp_list[time].push_back(o);
            }
        }

        void load_query(ifstream& infile){
            if(!infile){
                cout << "Cannot load query file!" << endl;
            }
            query current_query;
            int id_counter = 0;
            while(infile >> current_query.cordinate.x >> current_query.cordinate.y){
                current_query.id = id_counter;
                id_counter ++;
                query_list.push_back(current_query);
            }
        }

        void load_query_keywords(ifstream& infile){
            int query_index;
            string keyword;
            while(infile >> query_index >> keyword){
                query_list[query_index].query_keyword.emplace_back(keyword);
            }
        }

        /**
         * Starting from root node (entire game map), continue to break down in 4 children until the child size is smaller than the given grid size
         * @param current_node
         */
        void decompose_node(node &current_node){
            //if either x or y exceeds grid size
            if(current_node.cordinate.max_point.x - current_node.cordinate.min_point.x > grid_size or current_node.cordinate.max_point.y - current_node.cordinate.min_point.y > grid_size){
                double xMid = (current_node.cordinate.max_point.x + current_node.cordinate.min_point.x) / 2.0;
                double yMid = (current_node.cordinate.max_point.y + current_node.cordinate.min_point.y) / 2.0;
                node child_1;
                child_1.cordinate.min_point = current_node.cordinate.min_point;
                child_1.cordinate.max_point.x = xMid;
                child_1.cordinate.max_point.y = yMid;
                child_1.parent = current_node.nodeID;
                child_1.nodeID = ehl_tree.size();
                child_1.level = current_node.level + 1;
                child_1.distance = INF;
                child_1.has_object = false;
                ehl_tree.push_back(child_1);
                decompose_node(child_1);

                node child_2;
                child_2.cordinate.min_point.x = xMid;
                child_2.cordinate.min_point.y = current_node.cordinate.min_point.y;
                child_2.cordinate.max_point.x = current_node.cordinate.max_point.x;
                child_2.cordinate.max_point.y = yMid;
                child_2.parent = current_node.nodeID;
                child_2.nodeID = ehl_tree.size();
                child_2.level = current_node.level + 1;
                child_2.distance = INF;
                child_2.has_object = false;
                ehl_tree.push_back(child_2);
                decompose_node(child_2);

                node child_3;
                child_3.cordinate.min_point.x = current_node.cordinate.min_point.x;
                child_3.cordinate.min_point.y = yMid;
                child_3.cordinate.max_point.x = xMid;
                child_3.cordinate.max_point.y = current_node.cordinate.max_point.y;
                child_3.parent = current_node.nodeID;
                child_3.level = current_node.level + 1;
                child_3.nodeID = ehl_tree.size();
                child_3.distance = INF;
                child_3.has_object = false;
                ehl_tree.push_back(child_3);
                decompose_node(child_3);

                node child_4;
                child_4.cordinate.min_point.x = xMid;
                child_4.cordinate.min_point.y = yMid;
                child_4.cordinate.max_point = current_node.cordinate.max_point;
                child_4.parent = current_node.nodeID;
                child_4.level = current_node.level + 1;
                child_4.nodeID = ehl_tree.size();
                child_4.distance = INF;
                child_4.has_object = false;
                ehl_tree.push_back(child_4);
                decompose_node(child_4);

                ehl_tree[current_node.nodeID].child.top_left = child_3.nodeID;
                ehl_tree[current_node.nodeID].child.top_right = child_4.nodeID;
                ehl_tree[current_node.nodeID].child.bottom_left = child_1.nodeID;
                ehl_tree[current_node.nodeID].child.bottom_right = child_2.nodeID;
            }
            else{
                //if the node is at it's cell level, then it doesn't have any children
                ehl_tree[current_node.nodeID].child.top_left = -1;
                ehl_tree[current_node.nodeID].child.top_right = -1;
                ehl_tree[current_node.nodeID].child.bottom_left = -1;
                ehl_tree[current_node.nodeID].child.bottom_right = -1;
            }

        }

        /**
         * Add keyword for each object in the list to the Grid Tree
         */
        void initialise_object_keywords(){
            for(int i = 0; i < object_list.size(); ++i){
                traverse_keywords(i, ehl_tree[0]);
            }
        }

        /**
         * Traverse down the Grid Tree to add in the keyword from the root node to all of it's branches right till the base cell
         * @param index
         * @param current_node
         */
        void traverse_keywords(const int &index, const node &current_node) {
            //add keywords to the current node's keyword list
            ehl_tree[current_node.nodeID].has_object = true;
            for(int i = 0; i < object_keyword_list[index].size(); ++i){
                ehl_tree[current_node.nodeID].keyword_list[object_keyword_list[index][i]]++;
            }
            //if the current node has children
            if (current_node.child.top_left != -1) {
                Point center;
                center.x = (current_node.cordinate.max_point.x + current_node.cordinate.min_point.x) / 2;
                center.y = (current_node.cordinate.max_point.y + current_node.cordinate.min_point.y) / 2;
                //determine which node the object should be placed in according to the coordinates
                if (object_list[index].cordinate.max_point.x >= center.x) {
                    //Indicates right side
                    if (object_list[index].cordinate.max_point.y >= center.y) {
                        return traverse_keywords(index, ehl_tree[current_node.child.top_right]);
                    } else {
                        return traverse_keywords(index, ehl_tree[current_node.child.bottom_right]);
                    }
                } else {
                    if (object_list[index].cordinate.max_point.y >= center.y) {
                        return traverse_keywords(index, ehl_tree[current_node.child.top_left]);
                    } else {
                        return traverse_keywords(index, ehl_tree[current_node.child.bottom_left]);
                    }
                }
            }
            //add the object index to the current base cell's object list
            ehl_tree[current_node.nodeID].cell_object_list.push_back(index);
            object_list[index].level = current_node.nodeID;
        }


        int keyword_traverse_test(const object &p, const node &current_node) {
            if (current_node.child.top_left != -1) {
                Point center;
                center.x = (current_node.cordinate.max_point.x + current_node.cordinate.min_point.x) / 2;
                center.y = (current_node.cordinate.max_point.y + current_node.cordinate.min_point.y) / 2;

                if (p.cordinate.x >= center.x) {
                    //Indicates right side
                    if (p.cordinate.y >= center.y) {
                        return keyword_traverse_test(p, ehl_tree[current_node.child.top_right]);
                    } else {
                        return keyword_traverse_test(p, ehl_tree[current_node.child.bottom_right]);
                    }
                } else {
                    if (p.cordinate.y >= center.y) {
                        return keyword_traverse_test(p, ehl_tree[current_node.child.top_left]);
                    } else {
                        return keyword_traverse_test(p, ehl_tree[current_node.child.bottom_left]);
                    }
                }
            }
            return current_node.nodeID;
        }
        /**
         * Inserting object in Grid Tree without keyword information
         * @param p
         * @param current_node
         * @return
         */
        const node& traverse(const object &p, const node &current_node) {
            if (current_node.child.top_left != -1) {
                ehl_tree[current_node.nodeID].has_object = true;
                Point center;
                center.x = (current_node.cordinate.max_point.x + current_node.cordinate.min_point.x) / 2;
                center.y = (current_node.cordinate.max_point.y + current_node.cordinate.min_point.y) / 2;

                if (p.cordinate.x >= center.x) {
                    //Indicates right side
                    if (p.cordinate.y >= center.y) {
                        return traverse(p, ehl_tree[current_node.child.top_right]);
                    } else {
                        return traverse(p, ehl_tree[current_node.child.bottom_right]);
                    }
                } else {
                    if (p.cordinate.y >= center.y) {
                        return traverse(p, ehl_tree[current_node.child.top_left]);
                    } else {
                        return traverse(p, ehl_tree[current_node.child.bottom_left]);
                    }
                }
            }
            //ehl_tree[current_node.nodeID].keyword_list[to_string(p.id)]++;
            ehl_tree[current_node.nodeID].has_object = true;
            ehl_tree[current_node.nodeID].cell_object_list.push_back(p.id);
            return current_node;
        }


        /**
         * Calculate min distance of a given point and the current rectangle
         * @param Point p
         * @param current_rectangle
         * @return
         */
        double get_min_distance(Point &p, node &current_rectangle){
            double dist = 0.0;
            double min_i = 0.0; // minimum distance in dimension i

            min_i = min(current_rectangle.cordinate.max_point.x - p.x, p.x - current_rectangle.cordinate.min_point.x);
            if(min_i >= 0){
                dist += 0;
            }
            else{
                dist += pow(double(min_i), 2.0);
            }

            min_i = min(current_rectangle.cordinate.max_point.y - p.y, p.y - current_rectangle.cordinate.min_point.y);
            if(min_i >= 0){
                dist += 0;
            }
            else{
                dist += pow(double(min_i), 2.0);
            }
            dist = sqrt(dist);
            return dist;
        }

        /**
         * Retrieve shortest distance through EHL for two given points
         * @param q
         * @param c
         * @return
         */
        double get_ehl_distance(Point &q, Point &c){
            ehl_query->set_start_goal(q,c);
            ehl_query->search();
            double distance = ehl_query->get_cost();
            return distance;
        }

        void test(vector<double> &time){
            warthog::timer timer =  warthog::timer ();
            int query_id = 0;
            for(int i = 0; i < timestamp_num; ++i){
                initialise_tree(height,width,grid_size);
                for(int j = 0; j < object_list.size(); ++j){
                    object_list[j].cordinate.max_point = object_timestamp_list[i][j].new_location;
                }
                initialise_object_keywords();
                int counter;
                for(int j = query_id; j < query_id + 100; ++j){
                    timer.start();
                    keyword_knn(j);
                    timer.stop();
                    time[j] = timer.elapsed_time_micro();
                    reverse(query_list[j].knn_object.begin(),query_list[j].knn_object.end());
                    counter = j;
                }
                query_id = counter + 1;
            }
        }
        /**
         * Keyword experiment ran for all queries of a given map
         * @param time
         * @param distances
         * @param total_candidates
         * @param update_time
         * @param updated
         */
        void keyword_knn_experiment(vector<double> &time, vector<double> &distances, vector<int> &total_candidates, vector<double> &update_time, vector<int> &updated){
            warthog::timer timer =  warthog::timer ();
            warthog::timer update_timer =  warthog::timer ();
            int query_id = 0;
            //iterates over all the timestamp
            for(int i = 0; i < timestamp_num; ++i){
                //Updates for each timestamp after 0
                if(i > 0){
                    int number_of_changes = 0;
                    update_timer.start();
                    for(int j = 0; j < object_timestamp_list[i].size(); ++j){
                        bool changed = keyword_tree_update(j,object_timestamp_list[i][j].new_location);
                        if(changed){
                            number_of_changes++;
                        }
                    }
                    update_timer.stop();
                    update_time[i] = update_timer.elapsed_time_micro();
                    updated[i] = number_of_changes;
                }
                int counter;
                //query id increments by 100 for each timestamp
                for(int j = query_id; j < query_id + 100; ++j){
                    timer.start();
                    keyword_knn(j);
                    timer.stop();
                    time[j] = timer.elapsed_time_micro();
                    reverse(query_list[j].knn_object.begin(),query_list[j].knn_object.end());
                    counter = j;
                }
                query_id = counter + 1;
            }
        }

        /**
         * Keyword insertion experiment
         * @param insertion_time
         * @param deletion_time
         * @param probability
         */
        void keyword_insertion_experiment(vector<double> &insertion_time, vector<double> &deletion_time, double probability){
            warthog::timer insert_timer =  warthog::timer ();
            warthog::timer delete_timer =  warthog::timer ();
            //calculate number of objects to be inserted and deleted
            int changed = ceil((probability * added_objects.size()) / 2);
            for(int time = 0; time < timestamp_num; ++time){
                vector<int> newly_inserted;
                //find enough objects to insert
                while(newly_inserted.size() < changed){
                    int index = rand() % free_objects.size();
                    newly_inserted.push_back(free_objects[index]);
                    //delete from free object lists
                    free_objects.erase(free_objects.begin()+index);
                }
                //find enough objects to delete
                vector<int> deleted_objects;
                while(deleted_objects.size() < changed){
                    int index = rand() % added_objects.size();
                    deleted_objects.push_back(added_objects[index]);
                    //delete from object already present in tree list
                    added_objects.erase(added_objects.begin()+index);
                }
                insert_timer.start();
                //inserts
                for(int i = 0; i < newly_inserted.size(); ++i){
                    traverse_keywords(newly_inserted[i], ehl_tree[0]);
                    added_objects.push_back(newly_inserted[i]);
                }
                insert_timer.stop();
                insertion_time.push_back(insert_timer.elapsed_time_micro());

                delete_timer.start();
                //deletes
                for(int i = 0; i < deleted_objects.size(); ++i){
                    object_deletion(deleted_objects[i]);
                    free_objects.push_back(deleted_objects[i]);
                }
                delete_timer.stop();
                deletion_time.push_back(delete_timer.elapsed_time_micro());
                //cout << insert_timer.elapsed_time_micro() << " " << delete_timer.elapsed_time_micro() << endl;
            }

        }

        /**
         * Finds the distance of all objects for a query to test correctness of GT-EHL
         * Without keyword information - static knn
         * @param distances
         * @param query_id
         * @param total_candidates
         * @return
         */
        bool knn_tree_test(double &distances,int query_id, int& total_candidates){
            vector<pair<double,int>> object_distance_list;
            for(int i = 0; i < object_timestamp_list[0].size(); ++i){
                double distance = get_ehl_distance(query_list[query_id].cordinate, object_list[i].cordinate.max_point);
                object_distance_list.push_back(make_pair(distance,i));
            }
            sort(object_distance_list.begin(), object_distance_list.end());
            for(int i = 0; i < k_num; ++i){
                knn_test.push_back(object_distance_list[i].second);
                //cout << object_distance_list[i].first << " " << object_distance_list[i].second << endl;
            }
            total_candidates = object_timestamp_list[0].size();
            distances = object_distance_list[0].first;
            return true;
        }

        /**
         * Finds the distance of all objects that matches the keyword for a given query
         * @param distances
         * @param query_id
         * @param total_candidates
         * @return
         */
        bool keyword_knn_test(double &distances,int query_id, int& total_candidates){
            vector<pair<double,int>> object_distance_list;
            for(int i = 0; i < object_list.size(); ++i){
                int query_counter = 0;
                for(int j = 0; j < query_list[query_id].query_keyword.size(); ++j){
                    if(find(object_keyword_list[i].begin(), object_keyword_list[i].end(), query_list[query_id].query_keyword[j]) != object_keyword_list[i].end()){
                        query_counter++;
                    }
                }
                if(query_counter == query_list[query_id].query_keyword.size()){
                    double distance = get_ehl_distance(query_list[query_id].cordinate, object_list[i].cordinate.max_point);
                    object_distance_list.push_back(make_pair(distance,i));
                }
            }
            sort(object_distance_list.begin(), object_distance_list.end());
            total_candidates = object_distance_list.size();
            distances = object_distance_list[0].first;
            if(object_distance_list.size() > k_num){
                for(int i = 0; i < k_num; ++i){
                    knn_test.push_back(object_distance_list[i].second);
//                    if(object_distance_list[i].first == INF){
//                        distances = INF;
//                    }
                }
            }
            else{
                for(int i = 0; i < object_distance_list.size(); ++i){
                    knn_test.push_back(object_distance_list[i].second);
                }
            }
            return true;
        }

        /**
         * Finds knn of objects without textual information
         * Uses Euclidean distance as lowerbound for objects
         * @param query_id
         * @param k
         * @param ehl_time
         * @param heap_push
         * @return
         */
        bool basic_knn_euclidean(int &query_id, int& k, int &ehl_time, int &heap_push){
            clear_queue();
            double query_score = INF;
            vector<node> object_in_heap;
            nodePtr root = &ehl_tree[0];
            min_heap.push(root);
            warthog::timer timer2 =  warthog::timer ();
            timer2.start();
            while(!min_heap.empty()){
                nodePtr current_node = min_heap.top(); min_heap.pop();
                if(current_node->distance >= query_score){
                    for(int i = 0; i < k; ++i){
                        nodePtr current_object = max_heap.top(); max_heap.pop();
                        query_list[query_id].knn_object.push_back(current_object->nodeID);
                    }
                    return true;
                }
                if(current_node->parent == -2){
                    if(max_heap.size() < k){
                        current_node->distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        ehl_time ++;
                        max_heap.push(current_node);
                        if(max_heap.size() == k){
                            query_score = max_heap.top()->distance;
                        }
                    }
                    else if(current_node->distance <= query_score){
                        current_node->distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        ehl_time ++;
                        if(current_node->distance < query_score){
                            max_heap.push(current_node);
                            max_heap.pop();
                            query_score = max_heap.top()->distance;
                        }
                    }
                }
                else if(current_node->child.top_left == -1){
                    for(int i = 0; i < current_node->cell_object_list.size(); ++i){
                        node object;
                        //euclidean distance calculated here
                        object.distance = query_list[query_id].cordinate.distance(object_list[current_node->nodeID].cordinate.max_point);
                        object.nodeID = current_node->cell_object_list[i];
                        object.parent = -2;
                        object.child.top_left = -5;
                        object_in_heap.push_back(object);
                        nodePtr object_ptr = &object;
                        min_heap.push(object_ptr);
                        heap_push++;
                    }
                }
                else{
                    //intermediate nodes
                    if(current_node->parent != -2){
                        //check whether the children has objects inside it or not
                        if(ehl_tree[current_node->child.top_left].has_object){
                            ehl_tree[current_node->child.top_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_left]);
                            if(ehl_tree[current_node->child.top_left].distance <= query_score){
                                nodePtr node_to_push = &ehl_tree[current_node->child.top_left];
                                min_heap.push(node_to_push);
                                heap_push++;
                            }
                        }
                        if(ehl_tree[current_node->child.top_right].has_object){
                            ehl_tree[current_node->child.top_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_right]);
                            if(ehl_tree[current_node->child.top_right].distance <= query_score){
                                nodePtr node_to_push = &ehl_tree[current_node->child.top_right];
                                min_heap.push(node_to_push);
                                heap_push++;
                            }
                        }
                        if(ehl_tree[current_node->child.bottom_left].has_object){
                            ehl_tree[current_node->child.bottom_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_left]);
                            if(ehl_tree[current_node->child.bottom_left].distance <= query_score){
                                nodePtr node_to_push = &ehl_tree[current_node->child.bottom_left];
                                min_heap.push(node_to_push);
                                heap_push++;
                            }
                        }
                        if(ehl_tree[current_node->child.bottom_right].has_object){
                            ehl_tree[current_node->child.bottom_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_right]);
                            if(ehl_tree[current_node->child.bottom_right].distance <= query_score){
                                nodePtr node_to_push = &ehl_tree[current_node->child.bottom_right];
                                min_heap.push(node_to_push);
                                heap_push++;
                            }
                        }
                    }
                }

            }
            for(int i = 0; i < k; ++i){
                nodePtr current_object = max_heap.top(); max_heap.pop();
                query_list[query_id].knn_object.push_back(current_object->nodeID);
            }
            return true;
        }

        /**
         * Deleting object in Grid Tree
         * @param object_id
         */
        //TODO could optimise the insertion and deletion of using object pointers instead of node id and object id
        void object_deletion(const int &object_id){
            int cell_id = object_list[object_id].level;
            ehl_tree[cell_id].cell_object_list.erase(remove(ehl_tree[cell_id].cell_object_list.begin(),ehl_tree[cell_id].cell_object_list.end(),object_id));
            for(int i = 0; i < object_keyword_list[object_id].size(); ++i){
                if(ehl_tree[cell_id].keyword_list[object_keyword_list[object_id][i]] > 1){
                    ehl_tree[cell_id].keyword_list[object_keyword_list[object_id][i]]--;
                }
                else{
                    ehl_tree[cell_id].keyword_list.erase(object_keyword_list[object_id][i]);
                }
            }
            //delete all entries of keyword until it reaches root.
            int parent_index = ehl_tree[cell_id].parent;
            while(parent_index > 0){
                for(int i = 0; i < object_keyword_list[object_id].size(); ++i){
                    if(ehl_tree[parent_index].keyword_list[object_keyword_list[object_id][i]] > 1){
                        ehl_tree[parent_index].keyword_list[object_keyword_list[object_id][i]]--;
                    }
                    else{
                        ehl_tree[parent_index].keyword_list.erase(object_keyword_list[object_id][i]);
                    }
                }
                parent_index = ehl_tree[parent_index].parent;
            }
            object_list[object_id].level = -1;
        }

        /**
         * Update Grid Tree when an object moves in a timestamp
         * @param object_id
         * @param new_location
         * @return
         */
        bool keyword_tree_update(const int &object_id, Point &new_location){
            int cell_id = object_list[object_id].level;
            object_list[object_id].cordinate.max_point = new_location;
            //if the object moved away from the previous grid cell
            if(get_min_distance(new_location, ehl_tree[cell_id]) != 0){
                //remove object id from previous location
                ehl_tree[cell_id].cell_object_list.erase(remove(ehl_tree[cell_id].cell_object_list.begin(),ehl_tree[cell_id].cell_object_list.end(),object_id));
                for(int i = 0; i < object_keyword_list[object_id].size(); ++i){
                    if(ehl_tree[cell_id].keyword_list[object_keyword_list[object_id][i]] > 1){
                        ehl_tree[cell_id].keyword_list[object_keyword_list[object_id][i]]--;
                    }
                    else{
                        ehl_tree[cell_id].keyword_list.erase(object_keyword_list[object_id][i]);
                    }
                }

                //delete keyword entries in previous location's parents
                int parent_index = ehl_tree[cell_id].parent;
                while(get_min_distance(new_location, ehl_tree[parent_index]) != 0){
                    for(int i = 0; i < object_keyword_list[object_id].size(); ++i){
                        if(ehl_tree[parent_index].keyword_list[object_keyword_list[object_id][i]] > 1){
                            ehl_tree[parent_index].keyword_list[object_keyword_list[object_id][i]]--;
                        }
                        else{
                            ehl_tree[parent_index].keyword_list.erase(object_keyword_list[object_id][i]);
                        }

                    }
                    parent_index = ehl_tree[parent_index].parent;
                }
                //object_list[object_id].cordinate.max_point = new_location;
                traverse_keywords(object_id, ehl_tree[parent_index]);
                return true;
            }
            return false;
        }
        bool knn(int &query_id){
            clear_queue();
            double query_score = INF;
            nodePtr root = &(ehl_tree[0]);
            min_heap.push(root);
            while(!min_heap.empty()){
                nodePtr current_node = min_heap.top(); min_heap.pop();
                //current_node->print(*current_node);
                if(current_node->distance >= query_score){
                    for(int i = 0; i < k_num; ++i){
                        nodePtr current_object = max_heap.top(); max_heap.pop();
                        query_list[query_id].knn_object.push_back(current_object->nodeID);
                    }
                    return true;
                }
                //node is object
                if(current_node->parent == -2){
//                    cout << "Object node pop out ";
//                    current_node->print(*current_node);
                    if(max_heap.size() < k_num){
                        object_list[current_node->nodeID].distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        nodePtr object_ptr = &(object_list[current_node->nodeID]);
                        max_heap.push(object_ptr);
                        if(max_heap.size() == k_num){
                            query_score = max_heap.top()->distance;
                        }
                    }
                    else if(current_node->distance <= query_score){
                        current_node->distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        if(current_node->distance < query_score){
                            max_heap.push(current_node);
                            max_heap.pop();
                            query_score = max_heap.top()->distance;
                        }
                    }
                }
                //node is cell
                else if(current_node->child.top_left == -1){
                    for(int i = 0; i < current_node->cell_object_list.size(); ++i){
                            //euclidean distance calculated here
                            object_list[current_node->cell_object_list[i]].distance = query_list[query_id].cordinate.distance(object_list[current_node->cell_object_list[i]].cordinate.max_point);
                            nodePtr object_ptr = &(object_list[current_node->cell_object_list[i]]);
                            //object.print(objects_in_heap[pointer_index]);
                            min_heap.push(object_ptr);
                    }
                }
                else{
                    //intermediate nodes
                    if(current_node->parent != -2){
                        //check whether the children has objects inside it or not
                        if(ehl_tree[current_node->child.top_left].has_object){
                            ehl_tree[current_node->child.top_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_left]);
                            if(ehl_tree[current_node->child.top_left].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.top_left]);
                                min_heap.push(node_to_push);
                            }
                        }
                        if(ehl_tree[current_node->child.top_right].has_object){
                            ehl_tree[current_node->child.top_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_right]);
                            if(ehl_tree[current_node->child.top_right].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.top_right]);
                                min_heap.push(node_to_push);
                            }
                        }
                        if(ehl_tree[current_node->child.bottom_left].has_object){
                            ehl_tree[current_node->child.bottom_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_left]);
                            if(ehl_tree[current_node->child.bottom_left].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.bottom_left]);
                                min_heap.push(node_to_push);
                            }
                        }
                        if(ehl_tree[current_node->child.bottom_right].has_object){
                            ehl_tree[current_node->child.bottom_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_right]);
                            if(ehl_tree[current_node->child.bottom_right].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.bottom_right]);
                                min_heap.push(node_to_push);
                            }
                        }
                    }
                }
            }
            int max_heap_size = max_heap.size();
            for(int i = 0; i < max_heap_size; ++i){
                nodePtr current_object = max_heap.top(); max_heap.pop();
                //current_object->print(*current_object);
                query_list[query_id].knn_object.push_back(current_object->nodeID);
            }
            return true;
        }

        bool keyword_knn(int &query_id){
            clear_queue();
            double query_score = INF;
            nodePtr root = &(ehl_tree[0]);
            min_heap.push(root);
            while(!min_heap.empty()){
                nodePtr current_node = min_heap.top(); min_heap.pop();
                //current_node->print(*current_node);
                if(current_node->distance >= query_score){
                    for(int i = 0; i < k_num; ++i){
                        nodePtr current_object = max_heap.top(); max_heap.pop();
                        query_list[query_id].knn_object.push_back(current_object->nodeID);
                    }
                    return true;
                }
                //node is object
                if(current_node->parent == -2){
//                    cout << "Object node pop out ";
//                    current_node->print(*current_node);
                    if(max_heap.size() < k_num){
                        object_list[current_node->nodeID].distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        nodePtr object_ptr = &(object_list[current_node->nodeID]);
                        max_heap.push(object_ptr);
                        if(max_heap.size() == k_num){
                            query_score = max_heap.top()->distance;
                        }
                    }
                    else if(current_node->distance <= query_score){
                        current_node->distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        if(current_node->distance < query_score){
                            max_heap.push(current_node);
                            max_heap.pop();
                            query_score = max_heap.top()->distance;
                        }
                    }
                }
                    //node is cell
                else if(current_node->child.top_left == -1){
                    for(int i = 0; i < current_node->cell_object_list.size(); ++i){
                        query_counter = 0;
                        for(int j = 0; j < query_list[query_id].query_keyword.size(); ++j){
                            if(object_list[current_node->cell_object_list[i]].keyword_list.find(query_list[query_id].query_keyword[j])!= object_list[current_node->cell_object_list[i]].keyword_list.end()){
                                query_counter++;
                            }
                        }
                        if(query_counter == query_list[query_id].query_keyword.size()){
                            object_list[current_node->cell_object_list[i]].distance = query_list[query_id].cordinate.distance(object_list[current_node->cell_object_list[i]].cordinate.max_point);
                            nodePtr object_ptr = &(object_list[current_node->cell_object_list[i]]);
                            //object.print(objects_in_heap[pointer_index]);
                            min_heap.push(object_ptr);
                        }
                    }
                }
                else{
                    //intermediate nodes
                    if(current_node->parent != -2){
                        //check whether the children has objects inside it or not
                        query_counter = 0;
                        for(int i = 0; i < query_list[query_id].query_keyword.size(); ++i){
                            if(ehl_tree[current_node->child.top_left].keyword_list.find(query_list[query_id].query_keyword[i])!= ehl_tree[current_node->child.top_left].keyword_list.end()){
                                query_counter++;
                            }
                        }
                        if(query_counter == query_list[query_id].query_keyword.size()){
                            ehl_tree[current_node->child.top_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_left]);
                            if(ehl_tree[current_node->child.top_left].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.top_left]);
                                min_heap.push(node_to_push);
                            }
                        }

                        query_counter = 0;
                        for(int i = 0; i < query_list[query_id].query_keyword.size(); ++i){
                            if(ehl_tree[current_node->child.top_right].keyword_list.find(query_list[query_id].query_keyword[i])!= ehl_tree[current_node->child.top_right].keyword_list.end()){
                                query_counter++;
                            }
                        }
                        if(query_counter == query_list[query_id].query_keyword.size()){
                            ehl_tree[current_node->child.top_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_right]);
                            if(ehl_tree[current_node->child.top_right].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.top_right]);
                                min_heap.push(node_to_push);
                            }
                        }

                        query_counter = 0;
                        for(int i = 0; i < query_list[query_id].query_keyword.size(); ++i){
                            if(ehl_tree[current_node->child.bottom_right].keyword_list.find(query_list[query_id].query_keyword[i])!= ehl_tree[current_node->child.bottom_right].keyword_list.end()){
                                query_counter++;
                            }
                        }
                        if(query_counter == query_list[query_id].query_keyword.size()){
                            ehl_tree[current_node->child.bottom_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_right]);
                            if(ehl_tree[current_node->child.bottom_right].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.bottom_right]);
                                min_heap.push(node_to_push);
                            }
                        }

                        query_counter = 0;
                        for(int i = 0; i < query_list[query_id].query_keyword.size(); ++i){
                            if(ehl_tree[current_node->child.bottom_left].keyword_list.find(query_list[query_id].query_keyword[i])!= ehl_tree[current_node->child.bottom_left].keyword_list.end()){
                                query_counter++;
                            }
                        }
                        if(query_counter == query_list[query_id].query_keyword.size()){
                            ehl_tree[current_node->child.bottom_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_left]);
                            if(ehl_tree[current_node->child.bottom_left].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.bottom_left]);
                                min_heap.push(node_to_push);
                            }
                        }
                    }
                }
            }
            int max_heap_size = max_heap.size();
            for(int i = 0; i < max_heap_size; ++i){
                nodePtr current_object = max_heap.top(); max_heap.pop();
                //current_object->print(*current_object);
                query_list[query_id].knn_object.push_back(current_object->nodeID);
            }
            return true;
        }

        bool keyword_static(int &query_id){
            clear_queue();
            double query_score = INF;
            nodePtr root = &(ehl_tree[0]);
            min_heap.push(root);
            while(!min_heap.empty()){
                nodePtr current_node = min_heap.top(); min_heap.pop();
                //current_node->print(*current_node);
                if(current_node->distance >= query_score){
                    for(int i = 0; i < k_num; ++i){
                        nodePtr current_object = max_heap.top(); max_heap.pop();
                        query_list[query_id].knn_object.push_back(current_object->nodeID);
                    }
                    return true;
                }
                //node is object
                if(current_node->parent == -2){
//                    cout << "Object node pop out ";
//                    current_node->print(*current_node);
                    if(max_heap.size() < k_num){
                        object_list[current_node->nodeID].distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        nodePtr object_ptr = &(object_list[current_node->nodeID]);
                        max_heap.push(object_ptr);
                        if(max_heap.size() == k_num){
                            query_score = max_heap.top()->distance;
                        }
                    }
                    else if(current_node->distance <= query_score){
                        current_node->distance = get_ehl_distance(query_list[query_id].cordinate,object_list[current_node->nodeID].cordinate.max_point);
                        if(current_node->distance < query_score){
                            max_heap.push(current_node);
                            max_heap.pop();
                            query_score = max_heap.top()->distance;
                        }
                    }
                }
                    //node is cell
                else if(current_node->child.top_left == -1){
                    for(int i = 0; i < current_node->cell_object_list.size(); ++i){
                            object_list[current_node->cell_object_list[i]].distance = query_list[query_id].cordinate.distance(object_list[current_node->cell_object_list[i]].cordinate.max_point);
                            nodePtr object_ptr = &(object_list[current_node->cell_object_list[i]]);
                            //object.print(objects_in_heap[pointer_index]);
                            min_heap.push(object_ptr);
                    }
                }
                else{
                    //intermediate nodes
                    if(current_node->parent != -2){
                        //check whether the children has objects inside it or not
                        if(ehl_tree[current_node->child.top_left].has_object){
                            ehl_tree[current_node->child.top_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_left]);
                            if(ehl_tree[current_node->child.top_left].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.top_left]);
                                min_heap.push(node_to_push);
                            }
                        }
                        if(ehl_tree[current_node->child.top_right].has_object){
                            ehl_tree[current_node->child.top_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.top_right]);
                            if(ehl_tree[current_node->child.top_right].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.top_right]);
                                min_heap.push(node_to_push);
                            }
                        }
                        if(ehl_tree[current_node->child.bottom_right].has_object){
                            ehl_tree[current_node->child.bottom_right].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_right]);
                            if(ehl_tree[current_node->child.bottom_right].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.bottom_right]);
                                min_heap.push(node_to_push);
                            }
                        }
                        if(ehl_tree[current_node->child.bottom_left].has_object){
                            ehl_tree[current_node->child.bottom_left].distance = get_min_distance(query_list[query_id].cordinate, ehl_tree[current_node->child.bottom_left]);
                            if(ehl_tree[current_node->child.bottom_left].distance <= query_score){
                                nodePtr node_to_push = &(ehl_tree[current_node->child.bottom_left]);
                                min_heap.push(node_to_push);
                            }
                        }
                    }
                }
            }
            int max_heap_size = max_heap.size();
            for(int i = 0; i < max_heap_size; ++i){
                nodePtr current_object = max_heap.top(); max_heap.pop();
                //current_object->print(*current_object);
                query_list[query_id].knn_object.push_back(current_object->nodeID);
            }
            return true;
        }

        void clear_queue(){
            pq empty;
            swap(min_heap, empty);

            object_pq object_empty;
            swap(max_heap, object_empty);
        }

        void grid_density_stats(vector<unordered_map<int,int>> &grid_density){
            cout << ehl_tree.size()<< endl;
            for(int i = 0; i < timestamp_num; ++i){
                //Updates for each timestamp after 0
                if(i > 0){
                    for(int j = 0; j < object_timestamp_list[i].size(); ++j){
                        bool changed = keyword_tree_update(j,object_timestamp_list[i][j].new_location);
                    }
                }
                for(int j = 0; j < ehl_tree.size(); ++j){
                    if(ehl_tree[j].child.top_left == -1){
                        int size = ehl_tree[j].cell_object_list.size();
                        grid_density[i][size]++;
                    }
                }
            }
        }
    };
}
#endif //EBHL_EHL_TREE_H
