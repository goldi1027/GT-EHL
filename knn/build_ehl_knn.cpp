//
// Created by Jinchun Du on 26/9/2022.
//
/*
 * Constructs and runs object search
 * Takes in argument of game map, k, object dnesity and GT grid size
 * Finds and returns the k nearest matching keyword for all the queries in that given game map
 */

#include <ehl_tree.h>
#include <visibleAreaSearchInstance.h>
#include "ebhl.h"
#include "ebhl_query_v2.h"
#include "scenario.h"

using namespace std;
namespace pl = polyanya;
using namespace pl;
pl::MeshPtr mp;
typedef EBHL* EBHL_Ptr;
pl:: EBHL_query_v2* ebhlQueryV2;

/**
 * Find and retrieves k nearest object matching textual information for all queries in a given game map for uniform distribution
 * @param dir
 * @param map
 * @param k
 * @param object_number
 * @param grid_size
 */
void build_ehl_knn_keywords(string dir, string map, int k, int object_number, int grid_size){
    cout << "Knn for keywords.. " << map << " k: " << k << " #objects: " << object_number << " grid: " << grid_size << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());

    int map_width, map_height;
    mp->get_grid_width_height(map_width,map_height);
    //x = width, y = height

    vector<pl::Point> turning_point;
    //corresponding vertices number in mesh_vertices
    vector<int> turning_vertices;
    //vertice location in polygon?
    vector<pl::PointLocation> turning_vertices_location;
    vector<pl::Point> obstacle_middle_point;
    int id = 0;
    int poly = 0;
    for( pl::Vertex& v : mp->mesh_vertices){
        if(v.is_turning_vertex && !v.is_ambig) {
            for (int polygon: v.polygons) {
                if (polygon != -1) {
                    //                    p.polygons.insert(polygon);
                    poly = polygon;
                }
            }
            // there is some issue here, old implementation assume that these vertex always inside one polygon only. it doesnt
            // use the correct type actually, I fix it here manually assign a random adjacent polygon
            pl::PointLocation location = {pl::PointLocation::ON_CORNER_VERTEX_UNAMBIG, poly, -1, id, -1};
            turning_vertices_location.push_back(location);
            turning_vertices.push_back(id);
            turning_point.push_back(mp->mesh_vertices[id].p);
            pl::Vertex o1 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[0]];
            pl::Vertex o2 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[1]];
            pl::Point m  = {(o1.p.x+ o2.p.x)/2,
                            (o1.p.y+ o2.p.y)/2
            };
            obstacle_middle_point.push_back(m);
        }
        id ++;
    }
    turning_point.shrink_to_fit();
    obstacle_middle_point.shrink_to_fit();

    //Initialises EHL with grid size 1
    auto ebhl = new pl::EBHL(1,mp->get_map_height(),mp->get_map_width());

    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_1";
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    //initialises EHl query for distance retrieval later
    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //initialise GT-EHL
    auto ehl_knn = new pl::EhlKnn(map_height, map_width, grid_size, ebhlQueryV2);
    ehl_knn->k_num = k;
    ehl_knn->timestamp_num = 50;

    //loads object location for all timestamps
    string knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_knn.csv";
    ifstream knn_file(knn_path);
    ehl_knn->load_object(ehl_knn->timestamp_num,knn_file);

    //load object keywords
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    ehl_knn->load_keywords(object_keyword_file);

    ehl_knn->initialise_object_keywords();

    //load query locations
    string query_path = "dataset/knn_dataset/"+dir+"/location/"+map +"_query.csv";
    ifstream queryfile(query_path);
    ehl_knn->load_query(queryfile);
    ehl_knn->query_list.shrink_to_fit();

    //load query keywords
    string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number)  +"_query_keyword.csv";
    ifstream query_keyword_file(query_keyword_path);
    ehl_knn->load_query_keywords(query_keyword_file);


    cout << "Begin keyword knn computation" << endl;
    vector<vector<int>> knn_list;
    vector<double> time;
    vector<double> update_time;
    update_time.resize(ehl_knn->timestamp_num);
    knn_list.resize(ehl_knn->query_list.size());
    time.resize(ehl_knn->query_list.size());
    //fill(time.begin(), time.end(),0);

    vector<double> distances(ehl_knn->query_list.size());
    vector<int> total_candidates(ehl_knn->query_list.size());
    vector<int> updated(ehl_knn->query_list.size());
    warthog::timer timer =  warthog::timer ();

    //runs experiment
    ehl_knn->keyword_knn_experiment(time,distances,total_candidates,update_time, updated);


    //output experiment results
    string dataset_output = "dataset/result/knn/" + dir + "/"  + map + "_" + to_string(k) + "_" + to_string(object_number) + "_keyword_" + to_string(grid_size) +".csv";
    std::ofstream q_file(dataset_output);
    cout << dataset_output << endl;
    q_file<<"Query ID,Total Time\n";
    for (int i = 0; i < knn_list.size(); ++i) {
        q_file << std::fixed << setprecision(8)  << i  << " " << (double)time[i]  <<"\n";

    }
    q_file.close();

    string update_output = "dataset/result/knn/" + dir  + "/" + map + "_" + to_string(k) + "_" + to_string(object_number)  + "_rect4_update_" + to_string(grid_size) +".csv";
    std::ofstream u_file(update_output);
    cout << update_output << endl;
    u_file<<"map,Timestamp,Total Time,#Updates\n";
    for (int i = 0; i < update_time.size(); ++i) {
        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i]  << " " << updated[i] << "\n";

    }
    u_file.close();

    delete ebhl;
    delete ehl_knn;
    delete ebhlQueryV2;
}

/**
 * Find and retrieves k nearest object matching textual information for all queries in a given game map for different distribution (For example, clustered distribution)
 * @param dir
 * @param map
 * @param k
 * @param object_number
 * @param grid_size
 * @param rectangle (denotes the clustered distribution, one rectangle meaning there is one cluster of objects in the given game map)
 */
void build_ehl_knn_keywords_rectangles(string dir, string map, int k, int object_number, int grid_size, int rectangle){
    cout << "Knn for keywords.. " << map << " k: " << k << " #objects: " << object_number << " grid: " << grid_size << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());

    int map_width, map_height;
    mp->get_grid_width_height(map_width,map_height);
    //x = width, y = height

    vector<pl::Point> turning_point;
    //corresponding vertices number in mesh_vertices
    vector<int> turning_vertices;
    //vertice location in polygon?
    vector<pl::PointLocation> turning_vertices_location;
    vector<pl::Point> obstacle_middle_point;
    int id = 0;
    int poly = 0;
    for( pl::Vertex& v : mp->mesh_vertices){
        if(v.is_turning_vertex && !v.is_ambig) {
            for (int polygon: v.polygons) {
                if (polygon != -1) {
                    //                    p.polygons.insert(polygon);
                    poly = polygon;
                }
            }
            // there is some issue here, old implementation assume that these vertex always inside one polygon only. it doesnt
            // use the correct type actually, I fix it here manually assign a random adjacent polygon
            pl::PointLocation location = {pl::PointLocation::ON_CORNER_VERTEX_UNAMBIG, poly, -1, id, -1};
            turning_vertices_location.push_back(location);
            turning_vertices.push_back(id);
            turning_point.push_back(mp->mesh_vertices[id].p);
            pl::Vertex o1 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[0]];
            pl::Vertex o2 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[1]];
            pl::Point m  = {(o1.p.x+ o2.p.x)/2,
                            (o1.p.y+ o2.p.y)/2
            };
            obstacle_middle_point.push_back(m);
        }
        id ++;
    }
    turning_point.shrink_to_fit();
    obstacle_middle_point.shrink_to_fit();
    //initialises EHL with grid size 1
    auto ebhl = new pl::EBHL(1,mp->get_map_height(),mp->get_map_width());

    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_1";
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    //initialises EHL query
    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //initialise ehl_knn
    auto ehl_knn = new pl::EhlKnn(map_height, map_width, grid_size, ebhlQueryV2);
    ehl_knn->k_num = k;
    ehl_knn->timestamp_num = 50;

    //loads object location for all timestamp in a given distribution (rectangle)
    string knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_" + to_string(rectangle) + "_knn.csv";
    ifstream knn_file(knn_path);
    ehl_knn->load_object(ehl_knn->timestamp_num,knn_file);

    //load object keywords
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    ehl_knn->load_keywords(object_keyword_file);

    ehl_knn->initialise_object_keywords();


    //load query location
    string query_path = "dataset/knn_dataset/"+dir+"/location/"+map +"_query.csv";
    ifstream queryfile(query_path);
    ehl_knn->load_query(queryfile);
    ehl_knn->query_list.shrink_to_fit();

    //load query keywords
    string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number)  +"_query_keyword.csv";
    ifstream query_keyword_file(query_keyword_path);
    ehl_knn->load_query_keywords(query_keyword_file);


    cout << "Begin keyword knn computation" << endl;
    vector<vector<int>> knn_list;
    vector<double> time;
    vector<double> update_time;
    update_time.resize(ehl_knn->timestamp_num);
    knn_list.resize(ehl_knn->query_list.size());
    time.resize(ehl_knn->query_list.size());
    //fill(time.begin(), time.end(),0);

    vector<double> distances(ehl_knn->query_list.size());
    vector<int> total_candidates(ehl_knn->query_list.size());
    vector<int> updated(ehl_knn->query_list.size());
    warthog::timer timer =  warthog::timer ();

    //runs experiment
    ehl_knn->keyword_knn_experiment(time,distances,total_candidates,update_time, updated);


    //output experiment results
    string dataset_output = "dataset/result/knn/" + dir + "/rect/"  + dir + "_" + to_string(k) + "_" + to_string(object_number) + "_rect" + to_string(rectangle) + "_keyword_" + to_string(grid_size) +".csv";
    std::ofstream q_file(dataset_output,std::ofstream::out | std::ofstream::app);
    cout << dataset_output << endl;
   // q_file<<"Query ID,Total Time\n";
    for (int i = 0; i < knn_list.size(); ++i) {
        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i]  <<"\n";

    }
    q_file.close();

    string update_output = "dataset/result/knn/" + dir  + "/rect/" + dir + "_" + to_string(k) + "_" + to_string(object_number)  + "_rect" + to_string(rectangle) + "_update_" + to_string(grid_size) +".csv";
    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
    cout << update_output << endl;
    //u_file<<"map,Timestamp,Total Time,#Updates\n";
    for (int i = 0; i < update_time.size(); ++i) {
        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i]  << " " << updated[i] << "\n";

    }
    u_file.close();

    delete ebhl;
    delete ebhlQueryV2;
    delete ehl_knn;
}


void build_ehl_knn_insertion(string dir, string map, int k, int object_number, int grid_size, double probability){
    cout << "Knn for keywords.. " << map << " k: " << k << " #objects: " << object_number << " grid: " << grid_size << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());

    int map_width, map_height;
    mp->get_grid_width_height(map_width,map_height);
    //x = width, y = height

    vector<pl::Point> turning_point;
    //corresponding vertices number in mesh_vertices
    vector<int> turning_vertices;
    //vertice location in polygon?
    vector<pl::PointLocation> turning_vertices_location;
    vector<pl::Point> obstacle_middle_point;
    int id = 0;
    int poly = 0;
    for( pl::Vertex& v : mp->mesh_vertices){
        if(v.is_turning_vertex && !v.is_ambig) {
            for (int polygon: v.polygons) {
                if (polygon != -1) {
                    //                    p.polygons.insert(polygon);
                    poly = polygon;
                }
            }
            // there is some issue here, old implementation assume that these vertex always inside one polygon only. it doesnt
            // use the correct type actually, I fix it here manually assign a random adjacent polygon
            pl::PointLocation location = {pl::PointLocation::ON_CORNER_VERTEX_UNAMBIG, poly, -1, id, -1};
            turning_vertices_location.push_back(location);
            turning_vertices.push_back(id);
            turning_point.push_back(mp->mesh_vertices[id].p);
            pl::Vertex o1 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[0]];
            pl::Vertex o2 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[1]];
            pl::Point m  = {(o1.p.x+ o2.p.x)/2,
                            (o1.p.y+ o2.p.y)/2
            };
            obstacle_middle_point.push_back(m);
        }
        id ++;
    }
    turning_point.shrink_to_fit();
    obstacle_middle_point.shrink_to_fit();
    auto ebhl = new pl::EBHL(1,mp->get_map_height(),mp->get_map_width());

    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_1";
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //initialise ehl_knn
    auto ehl_knn = new pl::EhlKnn(map_height, map_width, grid_size, ebhlQueryV2);
    ehl_knn->k_num = k;
    ehl_knn->timestamp_num = 50;

    cout << "loading knn...." << endl;
    string knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_knn.csv";
    ifstream knn_file(knn_path);
    ehl_knn->load_object(ehl_knn->timestamp_num,knn_file);

    //load object keywords
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    ehl_knn->load_keywords(object_keyword_file);

    ehl_knn->initialise_object_keywords();

    //insertion objects
    string knn_insert_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_knn_insertion.csv";
    ifstream knn_insert_file(knn_insert_path);
    ehl_knn->load_insertion_object(knn_insert_file);

    string object_insertion_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword_insertion.csv";
    ifstream object_insertion_keyword_file(object_insertion_keyword_path);
    ehl_knn->load_insertion_keywords(object_insertion_keyword_file);



    //read query
    string query_path = "dataset/knn_dataset/"+dir+"/location/"+map +"_query.csv";
    ifstream queryfile(query_path);
    ehl_knn->load_query(queryfile);
    ehl_knn->query_list.shrink_to_fit();

    //!
    //load query objects
    string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_query_keyword.csv";
    //string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_" + to_string(query) +"_query_keyword.csv";
    ifstream query_keyword_file(query_keyword_path);
    ehl_knn->load_query_keywords(query_keyword_file);


    cout << "Begin keyword knn computation" << endl;
    vector<vector<int>> knn_list;
    vector<double> time;
    vector<double> update_time;
    update_time.resize(ehl_knn->timestamp_num);
    knn_list.resize(ehl_knn->query_list.size());
    time.resize(ehl_knn->query_list.size());
    //fill(time.begin(), time.end(),0);

    vector<double> distances(ehl_knn->query_list.size());
    vector<int> total_candidates(ehl_knn->query_list.size());
    vector<int> updated(ehl_knn->query_list.size());
    warthog::timer timer =  warthog::timer ();

    vector<double> insert_time;
    vector<double> delete_time;
    ehl_knn->keyword_insertion_experiment(insert_time, delete_time, probability);


    string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number) + "_"+ to_string(probability) + "_insertion_" + to_string(grid_size) +".csv";
    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
    cout << update_output << endl;
    //u_file<<"Timestamp,Total Time,#Updates\n";
    for (int i = 0; i < insert_time.size(); ++i) {
        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)insert_time[i]  << " " << (double)delete_time[i] << "\n";

    }
    u_file.close();

    delete ebhl;
}

void keywords_density_experiment(string dir, string map, int k, int object_number, int grid_size){
    cout << "Keyword experiment for grid density.. " << map << " k: " << k << " #objects: " << object_number << " grid: " << grid_size << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());

    int map_width, map_height;
    mp->get_grid_width_height(map_width,map_height);
    //x = width, y = height

    vector<pl::Point> turning_point;
    //corresponding vertices number in mesh_vertices
    vector<int> turning_vertices;
    //vertice location in polygon?
    vector<pl::PointLocation> turning_vertices_location;
    vector<pl::Point> obstacle_middle_point;
    int id = 0;
    int poly = 0;
    for( pl::Vertex& v : mp->mesh_vertices){
        if(v.is_turning_vertex && !v.is_ambig) {
            for (int polygon: v.polygons) {
                if (polygon != -1) {
                    //                    p.polygons.insert(polygon);
                    poly = polygon;
                }
            }
            // there is some issue here, old implementation assume that these vertex always inside one polygon only. it doesnt
            // use the correct type actually, I fix it here manually assign a random adjacent polygon
            pl::PointLocation location = {pl::PointLocation::ON_CORNER_VERTEX_UNAMBIG, poly, -1, id, -1};
            turning_vertices_location.push_back(location);
            turning_vertices.push_back(id);
            turning_point.push_back(mp->mesh_vertices[id].p);
            pl::Vertex o1 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[0]];
            pl::Vertex o2 = mp->mesh_vertices[mp->mesh_vertices[id].obstacle_edge[1]];
            pl::Point m  = {(o1.p.x+ o2.p.x)/2,
                            (o1.p.y+ o2.p.y)/2
            };
            obstacle_middle_point.push_back(m);
        }
        id ++;
    }
    turning_point.shrink_to_fit();
    obstacle_middle_point.shrink_to_fit();
    auto ebhl = new pl::EBHL(1,mp->get_map_height(),mp->get_map_width());

    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_1";
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //initialise ehl_knn
    auto ehl_knn = new pl::EhlKnn(map_height, map_width, grid_size, ebhlQueryV2);
    ehl_knn->k_num = k;
    ehl_knn->timestamp_num = 50;

    cout << "loading knn...." << endl;
    string knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_knn.csv";
    ifstream knn_file(knn_path);
    ehl_knn->load_object(ehl_knn->timestamp_num,knn_file);

    //load object keywords
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    ehl_knn->load_keywords(object_keyword_file);

    ehl_knn->initialise_object_keywords();



    //read query
    string query_path = "dataset/knn_dataset/"+dir+"/location/"+map +"_query.csv";
    ifstream queryfile(query_path);
    ehl_knn->load_query(queryfile);
    ehl_knn->query_list.shrink_to_fit();


    //load query objects
    //string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_query_keyword.csv";
    string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number)  +"_query_keyword.csv";
    ifstream query_keyword_file(query_keyword_path);
    ehl_knn->load_query_keywords(query_keyword_file);


    cout << "Begin keyword knn computation" << endl;
    vector<vector<int>> knn_list;
    vector<double> time;
    vector<double> update_time;
    update_time.resize(ehl_knn->timestamp_num);
    knn_list.resize(ehl_knn->query_list.size());
    time.resize(ehl_knn->query_list.size());
    //fill(time.begin(), time.end(),0);

    vector<double> distances(ehl_knn->query_list.size());
    vector<int> total_candidates(ehl_knn->query_list.size());
    vector<int> updated(ehl_knn->query_list.size());
    warthog::timer timer =  warthog::timer ();

    vector<unordered_map<int,int>> grid_density(ehl_knn->timestamp_num);
    ehl_knn->grid_density_stats(grid_density);

    for(int i = 0; i < grid_density.size(); ++i){
        for(auto j: grid_density[i]){
            cout << i << " " << j.first << " " << j.second << endl;
        }
    }

    delete ebhl;
}

int main(int argc, char*argv[]){
    try{
        string dir;
        string map;
        string k;
        string object_number;
        string grid_size;
        string rectangle;
        if(argc != 7){
            cerr << argv[0] << "Directory | Map | Grid Size" << endl;
            return 1;
        }
        else{
            dir = argv[1];
            map = argv[2];
            k = argv[3];
            object_number = argv[4];
            grid_size = argv[5];
            rectangle = argv[6];
        }
        //build_ehl_knn(dir, map, stoi(k), stoi(object_number), stoi(grid_size));

        //build_ehl_knn_keywords(dir, map, stoi(k), stoi(object_number), stoi(grid_size));
        build_ehl_knn_keywords_rectangles(dir, map, stoi(k), stoi(object_number), stoi(grid_size), stoi(rectangle));
        //keywords_density_experiment(dir, map, stoi(k), stoi(object_number), stoi(grid_size));
        //build_ehl_knn_insertion(dir, map, stoi(k), stoi(object_number), stoi(grid_size), stod(query));
        //build_ehl_knn_int_keywords(dir, map, stoi(k), stoi(object_number), stoi(grid_size));
    }catch(exception&err){
        cerr << "Stopped on exception lol: " << err.what() << endl;
    }
};