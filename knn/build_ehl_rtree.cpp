//
// Created by Jinchun Du on 8/12/2022.
//

#include "IRTreeEhl.h"
#include "ebhl.h"
#include "ebhl_query_v2.h"
#include <IERPolyanya.h>
#include <intervaHeuristic.h>

using namespace std;
namespace pl = polyanya;
using namespace pl;
pl::MeshPtr mp;
pl::IRTreeEhl* IE;
typedef EBHL* EBHL_Ptr;
pl:: EBHL_query_v2* ebhlQueryV2;
vector<Point> starts;
pl::SearchInstance* si;
pl::IntervalHeuristic* ki;
struct raw_object {
    int id;
    pl::Point previous_location;
    pl::Point new_location;
};
pl::IERPolyanya* ffp;
vector<vector<raw_object>>object_timestamp_list;
vector<vector<string>> query_keywords;
vector<vector<string>> object_keywords;
vector<Point> pts;
vector<Point> added_pts;
vector<int>object_list;
vector<int>free_objects;
vector<int> added_objects;
vector<vector<Point>> bulk_pts;
//mesh polygons
vector<vector<pl::Point>> polys;

void load_poly_mesh(string &polys_path, string &mesh_path) {
    //cin >> mesh_path >> polys_path >> obs_path >> pts_path;

    ifstream meshfile(mesh_path);
    ifstream polysfile(polys_path);

    mp = new pl::Mesh(meshfile);
    si = new pl::SearchInstance(mp);
    ki = new pl::IntervalHeuristic(mp);
    ffp = new pl::IERPolyanya(si);

}

void read_query_data(istream& infile){
    if(!infile){
        cout << "Cannot read query data!" << endl;
    }
    polyanya::Point current_point;
    while(infile >> current_point.x >> current_point.y){
        starts.push_back(current_point);
    }
    starts.shrink_to_fit();
}

void read_ehl_data(int timestamp,  istream& infile){
    if(!infile){
        cout << "Cannot read object data!!" << endl;
    }
    object_timestamp_list.resize(timestamp);
    int time;
    raw_object o;
    while (infile >> o.id >> time >> o.previous_location.x >> o.previous_location.y >> o.new_location.x >> o.new_location.y)
    {
        object_timestamp_list[time].push_back(o);
    }
    bulk_pts.resize(timestamp);
    pts.resize(object_timestamp_list[0].size());
    for(int i = 0; i < object_timestamp_list[0].size(); ++i){
        pts[i] = object_timestamp_list[0][i].previous_location;
    }
    for(int i = 0; i < object_timestamp_list.size(); ++i){
        for(int j = 0; j < object_timestamp_list[i].size(); ++j){
            bulk_pts[i].push_back(object_timestamp_list[i][j].new_location);
        }
    }
}

void load_insertion_object(istream& infile){
    if(!infile){
        cout << "Cannot load object file!" << endl;
    }
    int time;
    raw_object o;
    Point current_ob;
    while (infile >> o.id >> time >> o.previous_location.x>> o.previous_location.y >> o.new_location.x >> o.new_location.y)
    {
        object_timestamp_list[time].push_back(o);
    }
    int size = object_timestamp_list[0].size() /2;
    for(int i = size; i < object_timestamp_list[0].size(); ++i){
        added_pts.push_back(object_timestamp_list[0][i].previous_location);
    }
}

void load_query_keywords(ifstream& infile){
    if(!infile){
        cout << "Cannot read query keyword data!!" << endl;
    }
    int query_index;
    string keyword;
    while(infile >> query_index >> keyword){
        query_keywords[query_index].emplace_back(keyword);
    }
}

void load_object_keywords(ifstream& infile){
    if(!infile){
        cout << "Cannot read object keyword data!!" << endl;
    }
    int object_index;
    string keyword;
    while(infile >> object_index >> keyword){
        object_keywords[object_index].emplace_back(keyword);
    }
}


void load_insertion_keywords(istream& infile){
    if(!infile){
        cout << "Cannot load object file!" << endl;
    }
    for(int i = 0;i < pts.size(); ++i){
        object_list.push_back(i);
        added_objects.push_back(i);
    }

    int object_index;
    string keyword;
    int size = object_list.size();
    object_list.resize(object_timestamp_list[0].size());
    object_keywords.resize(object_timestamp_list[0].size());
    while(infile >> object_index >> keyword){
        object_keywords[object_index].emplace_back(keyword);
    }
    for(int i = size; i < object_list.size(); ++i){
        free_objects.push_back(i);
    }
}


void build_irtree_ehl(string dir, string map, int k, int object_number){
    cout << "Generating rtree knn for " << map << endl;
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

    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_" + to_string(1);
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_" + to_string(1);
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //read query file with keywords
    string scenario_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_query.csv";
    ifstream scenariofile(scenario_path);
    read_query_data(scenariofile);

    query_keywords.resize(starts.size());
    //string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_" + to_string(query) + "_query_keyword.csv";
    string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_query_keyword.csv";
    ifstream query_keyword_file(query_keyword_path);
    load_query_keywords(query_keyword_file);

    //object cordinate info
    pts.clear();
    string ehl_knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) +"_knn.csv";
    ifstream ehl_knn_file(ehl_knn_path);
    read_ehl_data(50,ehl_knn_file);
    ifstream rtree_knn_file(ehl_knn_path);

    //object keyword
    object_keywords.resize(object_timestamp_list[0].size());
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    load_object_keywords(object_keyword_file);

    auto IE = new pl::IRTreeEhl(ebhlQueryV2);
    IE->load_object_timestamp(50,rtree_knn_file);
    IE->set_keyword(query_keywords,object_keywords);
    IE->set_goals(pts);
    IE->set_K(k);

    vector<double> time;
    time.resize(starts.size());

    vector<double> rtime;
    rtime.resize(starts.size());

    vector<double> update_time;

    vector<double> distance;
    distance.resize(starts.size());

    vector<int> num_of_updates;

    vector<double> rtree_time;
    rtree_time.resize(starts.size());

    vector<int> num_in_heap(starts.size());

    warthog::timer timer;
    warthog::timer update_timer;

    //setup of IER-Pol
    si = new pl::SearchInstance(mp);
    ffp = new pl::IERPolyanya(si);
    ifstream ier_knn_file(ehl_knn_path);
    ffp->load_object_timestamp(50,ier_knn_file);
    ffp->set_keyword(query_keywords, object_keywords);
    ffp->set_goals(pts);
    ffp->set_K(k);

    warthog::timer ffp_timer;
    warthog::timer ffp_update_timer;
    vector<double>ffp_update_time(starts.size());
    vector<int> ffp_num_updates(starts.size());
    vector<double>ffp_time(starts.size());
    vector<double> ffp_rtree_time(starts.size());
    vector<double> test_pol_time(starts.size());


    ki = new pl::IntervalHeuristic(mp);
    ki->set_keyword(query_keywords,object_keywords);
    ki->set_goals(pts);
    ki->set_K(k);

    warthog::timer ih_timer;
    warthog::timer ih_update_timer;
    vector<double>ih_update_time(starts.size());
    vector<int> ih_num_updates(starts.size());
    vector<double>ih_time(starts.size());


    int query_id = 0;
    for(int i = 0; i < object_timestamp_list.size(); ++i){
        if(i > 0){
            int updates = 0;
            update_timer.start();
            IE->update_goals(bulk_pts[i],updates);
            update_timer.stop();
            update_time.push_back(update_timer.elapsed_time_micro());
            num_of_updates.push_back(updates);

            int ffp_updates = 0;
            ffp_update_timer.start();
            ffp->update_goals(bulk_pts[i],ffp_updates);
            ffp_update_timer.stop();
            ffp_update_time.push_back(ffp_update_timer.elapsed_time_micro());

            ffp_num_updates.push_back(ffp_updates);

            //interval heuristic update
            int ki_updates = 0;
            ih_update_timer.start();
            ki->update_goals(bulk_pts[i],ki_updates);
            ih_update_timer.stop();
            ih_update_time.push_back(ih_update_timer.elapsed_time_micro());
            ih_num_updates.push_back(ki_updates);
        }
        int counter;
        for(int j = query_id; j < query_id + 100; ++j){
            IE->set_start(starts[j],j);
            timer.start();
            vector<double>dist = IE->keyword_search();
            timer.stop();
            if(dist.empty()){
                distance[j] = INF;
            }
            else{
                distance[j] = dist[0];
            }
            rtree_time[j] = IE->rtree_cost;
            time[j] = timer.elapsed_time_micro();
            rtime[j] = IE->search_cost;

            ki->set_start(starts[j], j);
            ih_timer.start();
            int actual = ki->keyword_search();
            ih_timer.stop();
            ih_time[j] = ih_timer.elapsed_time_micro();

            //ierpolyanya
            ffp->set_start(starts[j], j);
            ffp_timer.start();
            vector<double> odists2 = ffp->keyword_search();
            ffp_timer.stop();
            ffp_time[j] = ffp_timer.elapsed_time_micro();
            ffp_rtree_time[j] = ffp->rtree_cost;
            test_pol_time[j] = ffp->search_cost;

            for(int i = 0; i < actual; ++i){
                int gid = ki->get_gid(i);
                double dist_ki = ki->get_cost(i);
                cout << dist_ki << " " << odists2[i] << " " << dist[i] << endl;
            }

            counter = j;
        }
        query_id = counter + 1;
    }

    //IER polyanya
//    si = new pl::SearchInstance(mp);
//    ffp = new pl::IERPolyanya(si);
//    ifstream ier_knn_file(ehl_knn_path);
//    ffp->load_object_timestamp(50,ier_knn_file);
//    ffp->set_keyword(query_keywords, object_keywords);
//    ffp->set_goals(pts);
//    ffp->set_K(k);
//
//    warthog::timer ffp_timer;
//    warthog::timer ffp_update_timer;
//    vector<double>ffp_update_time(starts.size());
//    vector<int> ffp_num_updates(starts.size());
//    vector<double>ffp_time(starts.size());
//    vector<double> ffp_rtree_time(starts.size());
//    vector<double> test_pol_time(starts.size());
//
//    query_id = 0;
//    for(int i = 0; i < object_timestamp_list.size(); ++i){
//        //cout << "Timestamp: " << i << endl;
//        //All timestamp except timestamp 0 needs update
//        if(i > 0){
//            int ffp_updates = 0;
//            ffp_update_timer.start();
//            ffp->update_goals(bulk_pts[i],ffp_updates);
//            ffp_update_timer.stop();
//            ffp_update_time.push_back(ffp_update_timer.elapsed_time_micro());
//            ffp_num_updates.push_back(ffp_updates);
//
//        }
//        //run query for each timestamp
//        //run 100 unique queries for each timestamp
//        //counter keeps track of next timestamp start
//        int counter;
//        for(int j = query_id; j < query_id + 100; ++j){
//            //ierpolyanya
//            ffp->set_start(starts[j], j);
//            ffp_timer.start();
//            vector<double> odists2 = ffp->keyword_search();
//            ffp_timer.stop();
//            ffp_time[j] = ffp_timer.elapsed_time_micro();
//            ffp_rtree_time[j] = ffp->rtree_cost;
//            test_pol_time[j] = ffp->search_cost;
//            counter = j;
//        }
//        query_id = counter + 1;
//    }

    string dataset_output = "dataset/result/knn/"  + dir + "/rect/" + dir + "_" + to_string(k) + "_" + to_string(object_number) +  "_comp.csv";
    std::ofstream q_file(dataset_output,std::ofstream::out | std::ofstream::app);
    cout << dataset_output << endl;
    //q_file<<"Map,Query ID,Total Time,EHL Time,Rtree Time,IER Time,IER Rtree Time\n";
    for (int i = 0; i <time.size(); ++i) {
        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i] << " " << rtime[i] << " " << rtree_time[i] << " " << ffp_time[i] << " " << ffp_rtree_time[i] <<
        " " << ih_time[i] << "\n";
    }
    q_file.close();


    string update_output = "dataset/result/knn/" + dir + "/rect/" + dir + "_" + to_string(k) + "_" + to_string(object_number)  + "_update_comp.csv";
    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
    cout << update_output << endl;
    //u_file<<"Map,Timestamp,Total Time,#Updates\n";
    for (int i = 0; i < update_time.size(); ++i) {
        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i] << " " << ffp_update_time[i] << " " << ih_update_time[i] << " " << num_of_updates[i] << "\n";
    }
    u_file.close();


//    string dataset_output = "dataset/result/knn/" + dir + "/"  + map + "_" + to_string(k) + "_" + to_string(object_number) +  "_rtree.csv";
//    std::ofstream q_file(dataset_output);
//    cout << dataset_output << endl;
//    q_file<<"Map,Query ID,Total Time,EHL Time\n";
//    for (int i = 0; i <time.size(); ++i) {
//        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i] << " " << rtime[i] <<"\n";
//    }
//    q_file.close();
//
//
//    string update_output = "dataset/result/knn/" + dir  + "/" + map + "_" + to_string(k) + "_" + to_string(object_number)  + "_update_rtree.csv";
//    std::ofstream u_file(update_output);
//    cout << update_output << endl;
//    u_file<<"Map,Timestamp,Total Time,#Updates\n";
//    for (int i = 0; i < update_time.size(); ++i) {
//        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i] << " " << num_of_updates[i] << "\n";
//    }
//    u_file.close();

//    string dataset_output = "dataset/result/knn/" + dir + "/"  + map + "_" + to_string(k) + "_" + to_string(object_number) +  "_4096_rtree.csv";
//    std::ofstream q_file(dataset_output);
//    cout << dataset_output << endl;
//    q_file<<"Query ID,Total Time,Distance\n";
//    for (int i = 0; i <time.size(); ++i) {
//        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i]  << " " << distance[i] <<"\n";
//    }
//    q_file.close();



    //string dataset_output = "dataset/result/knn/" + dir + "/"  + dir + "_" + to_string(k) + "_" + to_string(object_number) + "_q" + to_string(query) + "_rtree.csv";
//    string dataset_output = "dataset/result/knn/" + dir + "/"  + dir + "_" + to_string(k) + "_" + to_string(object_number) +  "_512_rtree.csv";
//    std::ofstream q_file(dataset_output,std::ofstream::out | std::ofstream::app);
//    cout << dataset_output << endl;
//    //q_file<<"Query ID,Total Time,Rtree Time\n";
//    for (int i = 0; i <time.size(); ++i) {
//        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i]  << " " << distance[i] <<"\n";
//    }
//    q_file.close();

    //string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number) +"_q" + to_string(probability)  + "_update_rtree.csv";
//    string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number)  + "_512_update_rtree.csv";
//    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
//    cout << update_output << endl;
//    //u_file<<"Timestamp,Total Time,#Updates\n";
//    for (int i = 0; i < update_time.size(); ++i) {
//        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i]  << "\n";
//    }
//    u_file.close();

}

void build_irtree_rect(string dir, string map, int k, int object_number, int rect){
    cout << "Generating rtree knn for " << map << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    string poly_path = "dataset/polygons/" + dir + "/" + map + ".map.poly";
    //load_poly_mesh(poly_path, mesh_path);
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

    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_" + to_string(1);
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_" + to_string(1);
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //read query file with keywords
    string scenario_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_query.csv";
    ifstream scenariofile(scenario_path);
    read_query_data(scenariofile);

    query_keywords.resize(starts.size());
    //string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_" + to_string(query) + "_query_keyword.csv";
    string query_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_query_keyword.csv";
    ifstream query_keyword_file(query_keyword_path);
    load_query_keywords(query_keyword_file);

    //object cordinate info
    pts.clear();
    string ehl_knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_" + to_string(rect)+ "_knn.csv";
    ifstream ehl_knn_file(ehl_knn_path);
    read_ehl_data(50,ehl_knn_file);
    ifstream rtree_knn_file(ehl_knn_path);

    //object keyword
    object_keywords.resize(object_timestamp_list[0].size());
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    load_object_keywords(object_keyword_file);

    auto IE = new pl::IRTreeEhl(ebhlQueryV2);
    IE->load_object_timestamp(50,rtree_knn_file);
    IE->set_keyword(query_keywords,object_keywords);
    IE->set_goals(pts);
    IE->set_K(k);

    vector<double> time;
    time.resize(starts.size());

    vector<double> rtime;
    rtime.resize(starts.size());

    vector<double> update_time;

    vector<double> distance;
    distance.resize(starts.size());

    vector<int> num_of_updates;

    vector<double> rtree_time;
    rtree_time.resize(starts.size());

    vector<int> num_in_heap(starts.size());

    warthog::timer timer;
    warthog::timer update_timer;

    //setup of IER-Pol
    si = new pl::SearchInstance(mp);
    ffp = new pl::IERPolyanya(si);
    ifstream ier_knn_file(ehl_knn_path);
    ffp->load_object_timestamp(50,ier_knn_file);
    ffp->set_keyword(query_keywords, object_keywords);
    ffp->set_goals(pts);
    ffp->set_K(k);

    warthog::timer ffp_timer;
    warthog::timer ffp_update_timer;
    vector<double>ffp_update_time;
    vector<int> ffp_num_updates;
    vector<double>ffp_time(starts.size());
    vector<double> ffp_rtree_time(starts.size());
    vector<double> test_pol_time(starts.size());


//    ki = new pl::IntervalHeuristic(mp);
//    ki->set_keyword(query_keywords,object_keywords);
//    ki->set_goals(pts);
//    ki->set_K(k);

    warthog::timer ih_timer;
    warthog::timer ih_update_timer;
    vector<double>ih_update_time;
    vector<int> ih_num_updates;
    vector<double>ih_time(starts.size());

    cout << "start search...." << endl;

    int query_id = 0;
    for(int i = 0; i < object_timestamp_list.size(); ++i){
        if(i > 0){
            int updates = 0;
            update_timer.start();
            IE->update_goals(bulk_pts[i],updates);
            update_timer.stop();
            update_time.push_back(update_timer.elapsed_time_micro());
            num_of_updates.push_back(updates);

            int ffp_updates = 0;
            ffp_update_timer.start();
            ffp->update_goals(bulk_pts[i],ffp_updates);
            ffp_update_timer.stop();
            ffp_update_time.push_back(ffp_update_timer.elapsed_time_micro());

            ffp_num_updates.push_back(ffp_updates);

            //interval heuristic update
//            int ki_updates = 0;
//            ih_update_timer.start();
//            ki->update_goals(bulk_pts[i],ki_updates);
//            ih_update_timer.stop();
//            ih_update_time.push_back(ih_update_timer.elapsed_time_micro());
//            ih_num_updates.push_back(ki_updates);
        }
        int counter;
        for(int j = query_id; j < query_id + 100; ++j){
            IE->set_start(starts[j],j);
            timer.start();
            vector<double>dist = IE->keyword_search();
            timer.stop();
            rtree_time[j] = IE->rtree_cost;
            time[j] = timer.elapsed_time_micro();
            rtime[j] = IE->search_cost;


            //ierpolyanya
            ffp->set_start(starts[j], j);
            ffp_timer.start();
            vector<double> odists2 = ffp->keyword_search();
            ffp_timer.stop();
            ffp_time[j] = ffp_timer.elapsed_time_micro();
            ffp_rtree_time[j] = ffp->rtree_cost;
            test_pol_time[j] = ffp->search_cost;


//            ki->set_start(starts[j], j);
//            ih_timer.start();
//            int actual = ki->keyword_search();
//            //int actual = 5;
//            ih_timer.stop();
//            ih_time[j] = ih_timer.elapsed_time_micro();


            //double dist_ki = ki->get_cost(0);

            if(!dist.empty()){
                if(fabs(odists2[0] - dist[0]) > EPSILON*10) {
                    cout << "WTF " << j  << " "  << odists2[0] << " " << dist[0] << endl;
                }
            }

            counter = j;

        }
        query_id = counter + 1;
    }

    string dataset_output = "dataset/result/knn/"  + dir + "/rect/" + dir + "_" + to_string(k) + "_" + to_string(object_number) + "_" + to_string(rect) +  "_comp.csv";
    std::ofstream q_file(dataset_output,std::ofstream::out | std::ofstream::app);
    cout << dataset_output << endl;
    //q_file<<"Map,Query ID,Total Time,EHL Time,Rtree Time,IER Time,IER Rtree Time\n";
    for (int i = 0; i <time.size(); ++i) {
        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i] << " " << rtime[i] << " " << rtree_time[i] << " " << ffp_time[i] << " " << ffp_rtree_time[i] << "\n";
    }
    q_file.close();


    string update_output = "dataset/result/knn/" + dir + "/rect/" + dir + "_" + to_string(k) + "_" + to_string(object_number) + "_" + to_string(rect) + "_update_comp.csv";
    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
    cout << update_output << endl;
    //u_file<<"Map,Timestamp,Total Time,#Updates\n";
    for (int i = 0; i < update_time.size(); ++i) {
        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i] << " " << ffp_update_time[i] << " " << num_of_updates[i] << "\n";
    }
    u_file.close();

//    string dataset_output = "dataset/result/knn/" + dir + "/"  + map + "_" + to_string(k) + "_" + to_string(object_number) +  "_rtree.csv";
//    std::ofstream q_file(dataset_output);
//    cout << dataset_output << endl;
//    q_file<<"Map,Query ID,Total Time,EHL Time\n";
//    for (int i = 0; i <time.size(); ++i) {
//        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i] << " " << rtime[i] <<"\n";
//    }
//    q_file.close();
//
//
//    string update_output = "dataset/result/knn/" + dir  + "/" + map + "_" + to_string(k) + "_" + to_string(object_number)  + "_update_rtree.csv";
//    std::ofstream u_file(update_output);
//    cout << update_output << endl;
//    u_file<<"Map,Timestamp,Total Time,#Updates\n";
//    for (int i = 0; i < update_time.size(); ++i) {
//        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i] << " " << num_of_updates[i] << "\n";
//    }
//    u_file.close();

//    string dataset_output = "dataset/result/knn/" + dir + "/"  + map + "_" + to_string(k) + "_" + to_string(object_number) +  "_4096_rtree.csv";
//    std::ofstream q_file(dataset_output);
//    cout << dataset_output << endl;
//    q_file<<"Query ID,Total Time,Distance\n";
//    for (int i = 0; i <time.size(); ++i) {
//        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i]  << " " << distance[i] <<"\n";
//    }
//    q_file.close();



    //string dataset_output = "dataset/result/knn/" + dir + "/"  + dir + "_" + to_string(k) + "_" + to_string(object_number) + "_q" + to_string(query) + "_rtree.csv";
//    string dataset_output = "dataset/result/knn/" + dir + "/"  + dir + "_" + to_string(k) + "_" + to_string(object_number) +  "_512_rtree.csv";
//    std::ofstream q_file(dataset_output,std::ofstream::out | std::ofstream::app);
//    cout << dataset_output << endl;
//    //q_file<<"Query ID,Total Time,Rtree Time\n";
//    for (int i = 0; i <time.size(); ++i) {
//        q_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)time[i]  << " " << distance[i] <<"\n";
//    }
//    q_file.close();

    //string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number) +"_q" + to_string(probability)  + "_update_rtree.csv";
//    string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number)  + "_512_update_rtree.csv";
//    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
//    cout << update_output << endl;
//    //u_file<<"Timestamp,Total Time,#Updates\n";
//    for (int i = 0; i < update_time.size(); ++i) {
//        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)update_time[i]  << "\n";
//    }
//    u_file.close();

}

void build_irtree_insert(string dir, string map, int k, int object_number, double probability){
    //object cordinate info
    pts.clear();
    //string ehl_knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_p" + to_string(probability) +"_knn.csv";
    string ehl_knn_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) +"_knn.csv";
    ifstream ehl_knn_file(ehl_knn_path);
    read_ehl_data(50,ehl_knn_file);
    ifstream rtree_knn_file(ehl_knn_path);

    //object keyword
    object_keywords.resize(object_timestamp_list[0].size());
    string object_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword.csv";
    ifstream object_keyword_file(object_keyword_path);
    load_object_keywords(object_keyword_file);

    //insertion objects
    string knn_insert_path = "dataset/knn_dataset/" + dir + "/location/" + map + "_" + to_string(object_number) + "_knn_insertion.csv";
    ifstream knn_insert_file(knn_insert_path);
    load_insertion_object(knn_insert_file);

    string object_insertion_keyword_path = "dataset/knn_dataset/" + dir + "/keyword/" + map + "_" + to_string(object_number) + "_object_keyword_insertion.csv";
    ifstream object_insertion_keyword_file(object_insertion_keyword_path);
    load_insertion_keywords(object_insertion_keyword_file);

    auto IE = new pl::IRTreeEhl(ebhlQueryV2);
    IE->load_object_timestamp(50,rtree_knn_file);
    IE->set_keyword(query_keywords,object_keywords);
    IE->set_goals(pts,added_pts);

    //add in additional object information for insertion and deletion
    IE->add_keyword(object_keywords);
    vector<double> insertion_time;

    vector<double> deletion_time;




    warthog::timer insert_timer =  warthog::timer ();
    warthog::timer delete_timer =  warthog::timer ();
    //calculate number of objects to be inserted and deleted
    int changed = ceil((probability * added_objects.size()) / 2);
    for(int time = 0; time < object_timestamp_list.size(); ++time) {
        vector<int> newly_inserted;
        //find enough objects to insert
        while (newly_inserted.size() < changed) {
            int index = rand() % free_objects.size();
            newly_inserted.push_back(free_objects[index]);
            //delete from free object lists
            free_objects.erase(free_objects.begin() + index);
        }
        //find enough objects to delete
        vector<int> deleted_objects;
        while (deleted_objects.size() < changed) {
            int index = rand() % added_objects.size();
            deleted_objects.push_back(added_objects[index]);
            //delete from object already present in tree list
            added_objects.erase(added_objects.begin() + index);
        }
        insert_timer.start();
        //inserts
        for (int i = 0; i < newly_inserted.size(); ++i) {
            IE->insert(newly_inserted[i]);
            added_objects.push_back(newly_inserted[i]);
        }
        insert_timer.stop();
        insertion_time.push_back(insert_timer.elapsed_time_micro());


        delete_timer.start();
        //deletes
        for (int i = 0; i < deleted_objects.size(); ++i) {
            IE->delete_object(deleted_objects[i]);
            free_objects.push_back(deleted_objects[i]);
        }
        delete_timer.stop();
        deletion_time.push_back(delete_timer.elapsed_time_micro());
    }


    string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number) +"_" + to_string(probability)  + "_insertion_rtree.csv";
    //string update_output = "dataset/result/knn/" + dir  + "/" + dir + "_" + to_string(k) + "_" + to_string(object_number)  + "_512_update_rtree.csv";
    std::ofstream u_file(update_output,std::ofstream::out | std::ofstream::app);
    cout << update_output << endl;
    //u_file<<"Timestamp,Total Time,#Updates\n";
    for (int i = 0; i < insertion_time.size(); ++i) {
        u_file << std::fixed << setprecision(8)  << map << " " << i  << " " << (double)insertion_time[i]  << " " << (double)deletion_time[i] << "\n";
    }
    u_file.close();

}
int main(int argc, char*argv[]){
    try{
        string dir;
        string map;
        string k;
        string object_number;
        string rectangle;

        if(argc != 6){
            cerr << argv[0] << "Directory | Map | Grid Size" << endl;
            return 1;
        }
        else{
            dir = argv[1];
            map = argv[2];
            k = argv[3];
            object_number = argv[4];
            rectangle = argv[5];

        }

        //build_rstartree_ehl(dir,map,stoi(k),stoi(object_number));

        //build_irtree_ehl(dir,map,stoi(k),stoi(object_number));

        build_irtree_rect(dir,map,stoi(k),stoi(object_number),stoi(rectangle));

        //build_irtree_insert(dir,map,stoi(k),stoi(object_number),stod(probability));
    }catch(exception&err){
        cerr << "Stopped on exception lol: " << err.what() << endl;
    }
};