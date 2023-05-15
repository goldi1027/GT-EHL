/*
 * Computes EHL query for a given scenario
 */
#include <stdio.h>
#include <iostream>
#include "point.h"
#include "consts.h"
#include "scenario.h"
#include "edge.h"
#include "coverage_ordering_path.h"
#include "ebhl.h"
#include "ebhl_query_v2.h"
#include "searchinstance.h"

using namespace std;
namespace pl = polyanya;
pl::MeshPtr mp;
pl::SearchInstance* si;
pl:: EBHL_query_v2* ebhlQueryV2;
int runtimes = 5;


//main query ran for EHL
int final_path_adjacent( string dir_name, string map_name,int grid_size, vector<double> &p_time, vector<double> &p_distance, vector<int> &path_length){
    std::cout<<"Current version: Final Experiment" << std::endl;
    vector<polyanya::Scenario> out;
    string mesh_path = "dataset/merged-mesh/"+dir_name+"/"+map_name +"-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    //test with polyanya
    si = new pl::SearchInstance(mp);

    string grid_path = "dataset/grid/"+dir_name+"/"+map_name +".map";

    // mark obstacle edge;
    mp->pre_compute_obstacle_edge_on_vertex();
    //wrap it in class later?
    mp->mark_turning_point(grid_path.c_str());
    //turning point is the actual point which would be used in search
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
    auto ebhl = new pl::EBHL(grid_size,mp->get_map_height(),mp->get_map_width());


    string grid_label_path = "dataset/ehl/" + dir_name + "/" + map_name + ".nt_adj_"+ to_string(grid_size);
    ebhl->load_adjacent_list(grid_label_path.c_str());

    string triangles_path = "dataset/ehl/" + dir_name + "/" + map_name + ".nt_triangles_"+ to_string(grid_size);
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir_name + "/" + map_name + ".label";
    string order = "dataset/hub_label/" + dir_name + "/" + map_name + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);


    string scenario_path = "dataset/scenarios/"+dir_name+"/"+map_name +".map.scen";
    ifstream scenariofile(scenario_path);
    polyanya::load_scenarios(scenariofile, out);
    warthog::timer timer =  warthog::timer ();

    p_distance.resize(out.size());
    path_length.resize(out.size());
    p_time.resize(out.size());
    fill(p_time.begin(), p_time.end(), 0);

    double unique_label_size;
    double error_count = 0;
    vector<double> pl_distance(out.size());
    fill(pl_distance.begin(),pl_distance.end(),0);

    int pl_query = 0;
    for(const polyanya::Scenario& s : out){
        si->set_start_goal(s.start,s.goal);
        si->search();
        double poly_search_cost =  si->get_cost();
        pl_distance[pl_query] = poly_search_cost;
        pl_query++;
    }
    for ( int i = 0; i < runtimes; i++) {
        int query_number = 0;
        for (const polyanya::Scenario& s : out) {
            timer.start();
            ebhlQueryV2->set_start_goal(s.start,s.goal);
            ebhlQueryV2->search_path();
            double distance = ebhlQueryV2->get_cost();
            auto current_path = ebhlQueryV2->get_path();
            timer.stop();

            p_time[query_number] = p_time[query_number] + timer.elapsed_time_micro();
            p_distance[query_number] = distance;
            path_length[query_number] = current_path.size();


            double epsilon = 0.01f;
            if (fabs(distance - pl_distance[query_number])>EPSILON*10) {
                // dont want to handle ambiguous case

                if(!ebhlQueryV2->is_start_or_target_ambiguous()){
//                                        vector<pl::Point> path;
//                                        si->get_path_points(path);
//                                        for(const auto& p : path){
//                                            std::cout<<p<<std::endl;
//                                        }
                    std::cout<<std::endl;
                    for(int i = 0; i < current_path.size(); ++i){
                        cout << current_path[i] << endl;
                    }
                    cout << "Error: "<< query_number << s.start << s.goal << pl_distance[query_number] << " " << distance << endl;
                    ++error_count;
                }
            }
            query_number++;
        }


    }
    double total = 0 ;
    for(double d_time  :p_time){
        total = total + d_time;
    }

    std::cout<<"ebhl: average time: "<<total/out.size()/runtimes <<" micro secs" <<std::endl;
    cout << "Error count: " << error_count/runtimes << " out of " << out.size() << endl;


    delete mp;
    delete ebhl;
    delete ebhlQueryV2;
    delete si;
    return 0;
}

void  benchmark_map(  string dir_name, string map_name, int grid_size){


    vector<double> p_time;
    vector<double> p_distance;
    vector<int> path_length;


   final_path_adjacent(dir_name, map_name, grid_size, p_time, p_distance, path_length);

// The below code is to generate the output of all queries for each map with information of runtime, path length, and distance per query. As it is for a single map only, for an entire benchmark, merging these files is necessary.
//    string query_output = "dataset/result/query/" + dir_name + "/final/" + to_string(grid_size) + "/" + map_name  +  ".csv";
//    std::ofstream q_file(query_output);
//    q_file<<"Map,Query Number,Time,Distance,Path Length\n";
//    for (int i = 0; i < p_distance.size(); ++i) {
//        q_file << std::fixed << setprecision(8) << map_name << "," << i << "," << p_time[i]/runtimes << "," << p_distance[i] << "," << path_length[i] <<"\n";
//    }
//    q_file.close();
}


int main(int argc, char* argv[]) {
    if (argc > 1) {
        string dir_name = string (argv[1]);
        string map_name = string(argv[2]);
        string grid_size = string(argv[3]);
        benchmark_map(dir_name, map_name, stoi(grid_size));

    }
}
