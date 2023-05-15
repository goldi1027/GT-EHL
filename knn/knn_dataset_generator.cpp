//
// Created by Jinchun Du on 22/9/2022.
// Generates dataset for EHL knn
//

/*
 * Generate k start and target objects
 * Generate a portion of the objects to be moving at each timestamp
 * Object info: Object ID | timestamp | previous location | new location
 * Generates different number of objects for each map
 * Generates 1000 queries with different keywords targeted at different object size
 */
#include <iosfwd>
#include <iostream>
#include <visibleAreaSearchInstance.h>
#include "ebhl.h"
#include "ebhl_query_v2.h"
#include "point.h"
#include <unordered_map>

using namespace std;
namespace pl = polyanya;
using namespace pl;
typedef EBHL* EBHL_Ptr;
pl:: EBHL_query_v2* ebhlQueryV2;
pl::MeshPtr mp;

struct object_info{
    int id;
    int timestamp;
    Point previous_location;
    Point new_location;
};

//Rectangle for cluster distribution
struct rectangle{
    Point min_point;
    Point max_point;
};

int map_width;
int map_height;
vector<vector<string>> keyword_list;

/**
 * Read in input keyword file and store it in a keyword list
 * @param infile
 */
void read_dataset(istream& infile){
    if(!infile){
        cout << "Cannot read data!" << endl;
    }
    string line;
    while(getline(infile,line)){
        vector<string> temp_keyword;
        istringstream ss(line);
        string word;
        while (ss >> word){
            temp_keyword.emplace_back(word);
        }
        keyword_list.push_back(temp_keyword);
    }
}

/**
 * Generates a random point in the traversable area of a given map
 * @return Random point
 */
Point generate_random_point(){
    Point obj_start;
    bool start_bool = true;
    //ensures the random points generated are in mesh
    while(start_bool){
        int x_start = rand() % map_width;
        int y_start = rand() % map_height;
        obj_start.x = double(x_start);
        obj_start.y = double(y_start);
        const PointLocation& start_pl = mp->get_point_location(obj_start);
        if(start_pl.type == 1){
            start_bool = false;
        }
    }
    return obj_start;
}

/**
 * Generate a random point inside a given rectangle dimensions
 * @param rectangles
 * @return Random point
 */

Point generate_random_point_in_rectangle(vector<rectangle> rectangles){
    Point obj_start;
    bool start_bool = true;
    int rectangle_index = rand() % rectangles.size();
    while(start_bool){
        //ensure the random point is within the randomly selected rectangle bounding box
        int x_range = int(rectangles[rectangle_index].max_point.x) - int(rectangles[rectangle_index].min_point.x) + 1;
        int y_range = int(rectangles[rectangle_index].max_point.y) - int(rectangles[rectangle_index].min_point.y) + 1;
        int x_start = rand() % x_range + int(rectangles[rectangle_index].min_point.x);
        int y_start = rand() % y_range + int(rectangles[rectangle_index].min_point.y);
        obj_start.x = double(x_start);
        obj_start.y = double(y_start);
        const PointLocation& start_pl = mp->get_point_location(obj_start);
        if(start_pl.type == 1){
            start_bool = false;
        }
    }
    if(obj_start.x > rectangles[rectangle_index].max_point.x && obj_start.x < rectangles[rectangle_index].min_point.x
    && obj_start.y > rectangles[rectangle_index].max_point.y && obj_start.y < rectangles[rectangle_index].min_point.y){
        cout << "Something wrong " << obj_start << " " << rectangles[rectangle_index].min_point << " " << rectangles[rectangle_index].max_point << endl;
    }

    return obj_start;
}

/**
 * Generates the trajectory path of an object
 * Finds a list of possible targets for each object, then the shortest path for each is found
 * This is to ensure the object has enough trajectory path to move in the timestamp given
 * The start and target could only be in the given rectangle's dimensions
 * @param  Rectangles
 * @param object path trajectory
 */

//find the trajectory of an object in different rectangles as source and targets
void calculate_location_rectangle(vector<Point> &intermediate_points, vector<rectangle> rectangles){
    vector<Point> target_list;
    //Finds a sufficient amount of targets for each object
    //The target should not be the same as the previous one (or they will not be able to move)
    //40 is just a number which will definitely generate sufficient trajectory for an object to move in 50 timestamps
    for(int i = 0; i < 40; ++i){
        Point random_point = generate_random_point_in_rectangle(rectangles);
        if(i > 0 && random_point!= target_list[i-1]){
            target_list.push_back(random_point);
        }
        else if(i == 0){
            target_list.push_back(random_point);
        }
    }

    //Shortest path of the object's targets are found
    vector<Point> shortest_path;
    for(int i = 1; i < target_list.size(); ++i){
        ebhlQueryV2->set_start_goal(target_list[i-1],target_list[i]);
        ebhlQueryV2->search_path();
        auto path = ebhlQueryV2->get_path();
        for(int j = 0; j < path.size()-1; ++j){
            shortest_path.push_back(path[j]);
        }
    }
    //Generates the trajectory from the shortest path
    //The trajectory is set to move at 1 Euclidean unit per timestamp
    //If the path ends but distance is less than 1, it just moves to the end of the current shortest path for this timestamp
    for(int i = 1; i < shortest_path.size(); ++i){
        double distance = shortest_path[i-1].distance(shortest_path[i]);
        double increment_x = (shortest_path[i].x - shortest_path[i-1].x) / distance;
        double increment_y = (shortest_path[i].y - shortest_path[i-1].y) / distance;
        //push start point into final path list
        Point start_point = shortest_path[i-1];
        bool flag = true;
        while(flag && distance > 0){
            intermediate_points.push_back(start_point);
            Point current_point;
            current_point.x = start_point.x + increment_x;
            current_point.y = start_point.y + increment_y;
            start_point = current_point;
            //Checks whether the move is less than the distance of the shortest path or not
            if(increment_x > 0){
                if(current_point.x >= shortest_path[i].x){
                    flag = false;
                }
            }else{
                if(current_point.x <= shortest_path[i].x){
                    flag = false;
                }
            }
        }

    }

}

void calculate_location(vector<Point> &intermediate_points){
    vector<Point> target_list;
    for(int i = 0; i < 40; ++i){
        Point random_point = generate_random_point();
        if(i > 0 && random_point!= target_list[i-1]){
            target_list.push_back(random_point);
        }
        else if(i == 0){
            target_list.push_back(random_point);
        }
    }
    vector<Point> shortest_path;
    for(int i = 1; i < target_list.size(); ++i){
        ebhlQueryV2->set_start_goal(target_list[i-1],target_list[i]);
        ebhlQueryV2->search_path();
        auto path = ebhlQueryV2->get_path();
        for(int j = 0; j < path.size()-1; ++j){
            shortest_path.push_back(path[j]);
        }
    }
    for(int i = 1; i < shortest_path.size(); ++i){
        double distance = shortest_path[i-1].distance(shortest_path[i]);
        double increment_x = (shortest_path[i].x - shortest_path[i-1].x) / distance;
        double increment_y = (shortest_path[i].y - shortest_path[i-1].y) / distance;
        //push start point into final path list
        Point start_point = shortest_path[i-1];
        bool flag = true;
        while(flag && distance > 0){
            intermediate_points.push_back(start_point);
            Point current_point;
            current_point.x = start_point.x + increment_x;
            current_point.y = start_point.y + increment_y;
            start_point = current_point;
            if(increment_x > 0){
                if(current_point.x >= shortest_path[i].x){
                    flag = false;
                }
            }else{
                if(current_point.x <= shortest_path[i].x){
                    flag = false;
                }
            }
        }

    }

}
/**
 * Generate a rectangle box that the at least one point is in the traversable area
 * @param rectangle_height
 * @param rectangle_width
 * @return rectangle
 */
Point generate_rectangle_start_point(int rectangle_height, int rectangle_width){
    Point rectangle_start;
    bool start_bool = true;
    //ensures the random points generated are in mesh
    while(start_bool){
        int x_start = rand() % map_width;
        int y_start = rand() % map_height;
        rectangle_start.x = double(x_start);
        rectangle_start.y = double(y_start);
        const PointLocation& start_pl = mp->get_point_location(rectangle_start);
        //checks whether extending the rectangle will exceed the map's dimensions
        if(start_pl.type == 1 && (rectangle_start.x + rectangle_width < map_width) && (rectangle_start.y + rectangle_height < map_height)){
            start_bool = false;
        }
    }
    return rectangle_start;
}

/**
 * Generate rectangle boxes of the traversable areas of a given map according to the ratio given
 *
 * @param ratio
 * @param num_rects
 * @return vector of rectangle boxes generated
 *
 */
vector<rectangle> generate_rectangles(int ratio, int num_rects){
    //rectangle dimensions are dependent on the ratio of a given map's dimensions
    int rectangle_height = ceil(map_height/ratio);
    int rectangle_width = ceil(map_width/ratio);
    vector<rectangle> rectangles;

    int current_count = 0;
    while(current_count < num_rects){
        Point rectangle_start = generate_rectangle_start_point(rectangle_height, rectangle_width);
        rectangle current_rectangle;
        current_rectangle.min_point = rectangle_start;
        current_rectangle.max_point.x = rectangle_start.x + rectangle_width;
        current_rectangle.max_point.y = rectangle_start.y + rectangle_height;
        rectangles.push_back(current_rectangle);
        current_count++;
    }
    return rectangles;
}

void generate_random_point(vector<Point> &query_points){
    Point random_point;
    bool target_bool = true;
    while(target_bool){
        int x_target = rand() % map_width;
        int y_target = rand() % map_height;
        random_point.x = double(x_target);
        random_point.y = double(y_target);
        const PointLocation& target_pl = mp->get_point_location(random_point);
        if(target_pl.type == 1){
            target_bool = false;
        }
    }
    query_points.push_back(random_point);
}
/**
 * Generates knn dataset with given map and number of objects(k)
 * @param dir
 * @param map
 * @param k
 */
void generate_knn_dataset(string dir, string map, string k){
    cout << "Generating knn dataset for " << map << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());
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


    //load grid size 1 EHL
    auto ebhl = new pl::EBHL(1,mp->get_map_height(),mp->get_map_width());
    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_"+ "1";
    ebhl->load_adjacent_list(grid_label_path.c_str());
    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_"+ "1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    //load knn keyword dataset for current directory
    string dataset_path = "dataset/knn_dataset/" + dir + "_cleaned.txt";
    ifstream dataset_file(dataset_path);
    read_dataset(dataset_file);

    //get number of traversable grids
    int traversable_grids = 0;
    mp->get_traversable_grids(grid_path.c_str(), traversable_grids);
    cout << traversable_grids << endl;
//    cout << (double)k/(double)100 << endl;
//    int object_size = ((double)k/100) * traversable_grids;
    cout <<k << " " << traversable_grids << endl;
    double k_num = stod(k);
    //number of objects is determined by the ratio of traversable grids
    int object_size = ceil(k_num * (double)traversable_grids);

    cout << object_size << endl;
    //generates the position of each object at each timestamp
    //stores it in a vector of defined object structure - Object ID | Previous Location | New Location
    //index defines timestamp
    vector<vector<Point>> object_timestamp_paths;
    object_timestamp_paths.resize(object_size);
    int timestamp_count = 50;

    for(int i = 0; i < object_size; ++i){
        vector<Point> intermediate_paths;
        calculate_location(intermediate_paths);
        object_timestamp_paths[i] = intermediate_paths;
    }

    //generates a start and target for each object first
    //shortest path is then calculated
    //if the path is less than the timestamp, iteratively generates new target and shortest path

    vector<vector<object_info>> final_object_dataset;
    final_object_dataset.resize(timestamp_count);
    for(int i = 0; i < object_timestamp_paths.size(); ++i){
        //probability of the object moving at that timestamp
        int previous_location_index = 0;
        for(int j = 0; j < timestamp_count; ++j){
            //probability
            bool moving = (rand() % 10) < 7;
            //bool moving = true;
            object_info current_object;
            current_object.id = i;
            current_object.timestamp = j;
            //initial location
            if(j == 0){
                current_object.previous_location = object_timestamp_paths[i][0];
                current_object.new_location = object_timestamp_paths[i][0];
            }else if(moving){
                current_object.previous_location = object_timestamp_paths[i][previous_location_index];
                current_object.new_location = object_timestamp_paths[i][previous_location_index + 1];
                previous_location_index++;
            }else{
                current_object.previous_location = object_timestamp_paths[i][previous_location_index];
                current_object.new_location = object_timestamp_paths[i][previous_location_index];
            }
            final_object_dataset[j].push_back(current_object);
        }
    }



    string dataset_output = "dataset/knn_dataset/" + dir + "/location/"  + map + "_" + to_string(int(k_num*100)) + "_knn.csv";
    std::ofstream q_file(dataset_output);
    cout << dataset_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < final_object_dataset.size(); ++i) {
        for(int j = 0; j < final_object_dataset[i].size(); ++j)
        q_file << std::fixed << setprecision(8)  << final_object_dataset[i][j].id << " " << final_object_dataset[i][j].timestamp << " " << final_object_dataset[i][j].previous_location.x << " " << final_object_dataset[i][j].previous_location.y
        << " " << final_object_dataset[i][j].new_location.x << " " << final_object_dataset[i][j].new_location.y <<"\n";
    }
    q_file.close();

    vector<vector<string>> object_keywords(object_size);

    //used for generating query keywords
    vector<vector<string>> query_object_index;

    for(int i = 0; i < object_size; ++i){
        //determine which keyword subset it should chose from
        int keyword_index = rand() % keyword_list.size();
        int number_of_keywords;
        //randomly assign a number of keywords from the subset
        //ensure if there's only one object, it has at least 2 keywords to fulfill query criteria
        if(object_size < 2){
            number_of_keywords = rand() % (keyword_list[keyword_index].size() - 2)+ 2;
        }
        else{
            number_of_keywords = rand() % (keyword_list[keyword_index].size() - 1) + 1;
        }
        unordered_map<string,int> keyword;
        vector<string> temp_keyword_set = keyword_list[keyword_index];
        int counter = temp_keyword_set.size();
        //keep generating keywords until it has met the keyword number
        while(counter!=number_of_keywords){
            int index = rand() % counter;
            temp_keyword_set.erase(temp_keyword_set.begin()+index);
            counter--;
        }
        for(int j = 0; j < temp_keyword_set.size(); ++j){
            object_keywords[i].emplace_back(temp_keyword_set[j]);
        }
        //a query should contain at least 2 keywords by default, so it is put into the pool for query keywords to be chosen from
        if(number_of_keywords > 1){
            query_object_index.push_back(object_keywords[i]);
        }

    }


    //output for object keyword
    string keyword_output = "dataset/knn_dataset/" + dir + "/keyword/"  + map + "_" + to_string(int(k_num*100)) + "_object_keyword.csv";
    std::ofstream k_file(keyword_output);
    cout << keyword_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < object_keywords.size(); ++i) {
        for(int j = 0; j < object_keywords[i].size(); ++j)
            k_file << fixed << setprecision(8) << i << " " << object_keywords[i][j] << "\n";
    }
    k_file.close();

    //generate query keywords
    int num_of_queries = 5000;
    vector<vector<string>> query_keywords(num_of_queries);
    for(int i = 0; i < num_of_queries; ++i) {
        unordered_map<string, int> keywords;
        //randomly chose a object keyword set for the query keywords to be chosen from
        int index = rand() % query_object_index.size();
        vector<string> temp_keyword_set = query_object_index[index];
        int counter = temp_keyword_set.size();
        while(counter!=2){
            int index = rand() % counter;
            temp_keyword_set.erase(temp_keyword_set.begin()+index);
            counter--;
        }
        for(int j = 0; j < temp_keyword_set.size(); ++j){
            query_keywords[i].emplace_back(temp_keyword_set[j]);
        }
    }

    string query_output = "dataset/knn_dataset/" + dir + "/keyword/"  + map  + "_" + to_string(int(k_num*100)) + "_query_keyword.csv";
    std::ofstream query_file(query_output);
    cout << query_output << endl;
    for (int i = 0; i < query_keywords.size(); ++i) {
        for(int j = 0; j < query_keywords[i].size(); ++j){
            query_file << fixed << setprecision(8) << i << " " << query_keywords[i][j] << endl;
        }
    }
    query_file.close();


}


void generate_knn_dataset_insertion(string dir, string map, string k){
    cout << "Generating knn dataset for " << map << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());
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
    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_"+ "1";
    ebhl->load_adjacent_list(grid_label_path.c_str());
    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_"+ "1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    string dataset_path = "dataset/knn_dataset/" + dir + "_cleaned.txt";
    ifstream dataset_file(dataset_path);
    read_dataset(dataset_file);

    int traversable_grids = 0;
    mp->get_traversable_grids(grid_path.c_str(), traversable_grids);
    cout << traversable_grids << endl;
//    cout << (double)k/(double)100 << endl;
//    int object_size = ((double)k/100) * traversable_grids;
    cout <<k << " " << traversable_grids << endl;
    double k_num = stod(k);
    int object_size = ceil(k_num * (double)traversable_grids) * 2;

    cout << object_size << endl;
    //generates the position of each object at each timestamp
    //stores it in a vector of defined object structure - Object ID | Previous Location | New Location
    //index defines timestamp
    vector<vector<Point>> object_timestamp_paths;
    object_timestamp_paths.resize(object_size);
    int timestamp_count = 50;

    for(int i = 0; i < object_size; ++i){
        vector<Point> intermediate_paths;
        calculate_location(intermediate_paths);
        object_timestamp_paths[i] = intermediate_paths;
    }

    //generates a start and target for each object first
    //shortest path is then calculated
    //if the path is less than the timestamp, iteratively generates new target and shortest path

    vector<vector<object_info>> final_object_dataset;
    final_object_dataset.resize(timestamp_count);
    for(int i = 0; i < object_timestamp_paths.size(); ++i){
        //probability of the object moving at that timestamp
        int previous_location_index = 0;
        for(int j = 0; j < timestamp_count; ++j){
            //probability
            bool moving = (rand() % 10) < 7;
            //bool moving = true;
            object_info current_object;
            current_object.id = i;
            current_object.timestamp = j;
            //initial location
            if(j == 0){
                current_object.previous_location = object_timestamp_paths[i][0];
                current_object.new_location = object_timestamp_paths[i][0];
            }else if(moving){
                current_object.previous_location = object_timestamp_paths[i][previous_location_index];
                current_object.new_location = object_timestamp_paths[i][previous_location_index + 1];
                previous_location_index++;
            }else{
                current_object.previous_location = object_timestamp_paths[i][previous_location_index];
                current_object.new_location = object_timestamp_paths[i][previous_location_index];
            }
            final_object_dataset[j].push_back(current_object);
        }
    }



    string dataset_output = "dataset/knn_dataset/" + dir + "/location/"  + map + "_" + to_string(int(k_num*100)) + "_knn.csv";
    std::ofstream q_file(dataset_output);
    cout << dataset_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < final_object_dataset.size(); ++i) {
        for(int j = 0; j < final_object_dataset[i].size()/2; ++j)
            q_file << std::fixed << setprecision(8)  << final_object_dataset[i][j].id << " " << final_object_dataset[i][j].timestamp << " " << final_object_dataset[i][j].previous_location.x << " " << final_object_dataset[i][j].previous_location.y
                   << " " << final_object_dataset[i][j].new_location.x << " " << final_object_dataset[i][j].new_location.y <<"\n";
    }
    q_file.close();

    string insertion_k_output = "dataset/knn_dataset/" + dir + "/location/"  + map + "_" + to_string(int(k_num*100)) + "_knn_insertion.csv";
    std::ofstream d_file(insertion_k_output);
    cout << insertion_k_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < final_object_dataset.size(); ++i) {
        for(int j = final_object_dataset[i].size()/2; j < final_object_dataset[i].size(); ++j)
            d_file << std::fixed << setprecision(8)  << final_object_dataset[i][j].id << " " << final_object_dataset[i][j].timestamp << " " << final_object_dataset[i][j].previous_location.x << " " << final_object_dataset[i][j].previous_location.y
                   << " " << final_object_dataset[i][j].new_location.x << " " << final_object_dataset[i][j].new_location.y <<"\n";
    }
    d_file.close();


//    vector<string> keyword_list {"strength","dexterity","constitution","intelligence","wisdom","charisma","acrobatics",
//                                 "animal","arcana","athletics","deception","history","intimidation","insight","investigation",
//                                 "medicine","nature","perception","performance","persuasion","religion","sleight","stealth","survival"};


    vector<vector<string>> object_keywords(object_size);

    //used for generating query keywords
    vector<vector<string>> query_object_index;

    for(int i = 0; i < object_size; ++i){
        int keyword_index = rand() % keyword_list.size();
        int number_of_keywords;
        //ensure if there's only one object, it has at least 2 keywords to fulfill query criteria
        if(object_size < 2){
            number_of_keywords = rand() % (keyword_list[keyword_index].size() - 2)+ 2;
        }
        else{
            number_of_keywords = rand() % (keyword_list[keyword_index].size() - 1) + 1;
        }
        unordered_map<string,int> keyword;
        vector<string> temp_keyword_set = keyword_list[keyword_index];
        int counter = temp_keyword_set.size();
        while(counter!=number_of_keywords){
            int index = rand() % counter;
            temp_keyword_set.erase(temp_keyword_set.begin()+index);
            counter--;
        }
        for(int j = 0; j < temp_keyword_set.size(); ++j){
            object_keywords[i].emplace_back(temp_keyword_set[j]);
        }
        if(number_of_keywords > 1){
            query_object_index.push_back(object_keywords[i]);
        }

    }



    string keyword_output = "dataset/knn_dataset/" + dir + "/keyword/"  + map + "_" + to_string(int(k_num*100)) + "_object_keyword.csv";
    std::ofstream k_file(keyword_output);
    cout << keyword_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < object_keywords.size()/2; ++i) {
        for(int j = 0; j < object_keywords[i].size(); ++j)
            k_file << fixed << setprecision(8) << i << " " << object_keywords[i][j] << "\n";
    }
    k_file.close();


    string insertion_output = "dataset/knn_dataset/" + dir + "/keyword/"  + map + "_" + to_string(int(k_num*100)) + "_object_keyword_insertion.csv";
    std::ofstream i_file(insertion_output);
    cout << insertion_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = object_keywords.size()/2; i < object_keywords.size(); ++i) {
        for(int j = 0; j < object_keywords[i].size(); ++j)
            i_file << fixed << setprecision(8) << i << " " << object_keywords[i][j] << "\n";
    }
    i_file.close();
}
//additional parameter of ratio of how big the rectangle should be, number of rectangles stands how many rectangles should be generated
//10 stands for 10%
void generate_rectangle_object_sets(string dir, string map, string k, string ratio, string num_rects){
    cout << "Generating rectangle object sets for " << map << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());
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
    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_"+ "1";
    ebhl->load_adjacent_list(grid_label_path.c_str());
    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_"+ "1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    string dataset_path = "dataset/knn_dataset/" + dir + "_cleaned.txt";
    ifstream dataset_file(dataset_path);
    read_dataset(dataset_file);

    int traversable_grids = 0;
    mp->get_traversable_grids(grid_path.c_str(), traversable_grids);
    cout << traversable_grids << endl;
    cout <<k << " " << traversable_grids << endl;
    double k_num = stod(k);
    int object_size = ceil(k_num * (double)traversable_grids);

    cout << object_size << endl;
    //generates the position of each object at each timestamp
    //stores it in a vector of defined object structure - Object ID | Previous Location | New Location
    //index defines timestamp
    vector<vector<Point>> object_timestamp_paths;
    object_timestamp_paths.resize(object_size);
    int timestamp_count = 50;


    //generate rectangles
    vector<rectangle> rectangles;
    rectangles = generate_rectangles(stoi(ratio),stoi(num_rects));

    for(int i = 0; i < rectangles.size(); ++i){
        cout << rectangles[i].min_point << " " << rectangles[i].max_point << endl;
    }

    //generate start and target within the rectangles boxes

    for(int i = 0; i < object_size; ++i){
        vector<Point> intermediate_paths;
        calculate_location_rectangle(intermediate_paths,rectangles);
        object_timestamp_paths[i] = intermediate_paths;
    }

    vector<vector<object_info>> final_object_dataset;
    final_object_dataset.resize(timestamp_count);
    for(int i = 0; i < object_timestamp_paths.size(); ++i){
        //probability of the object moving at that timestamp
        int previous_location_index = 0;
        for(int j = 0; j < timestamp_count; ++j){
            //probability
            bool moving = (rand() % 10) < 7;
            //bool moving = true;
            object_info current_object;
            current_object.id = i;
            current_object.timestamp = j;
            //initial location
            if(j == 0){
                current_object.previous_location = object_timestamp_paths[i][0];
                current_object.new_location = object_timestamp_paths[i][0];
            }else if(moving){
                current_object.previous_location = object_timestamp_paths[i][previous_location_index];
                current_object.new_location = object_timestamp_paths[i][previous_location_index + 1];
                previous_location_index++;
            }else{
                current_object.previous_location = object_timestamp_paths[i][previous_location_index];
                current_object.new_location = object_timestamp_paths[i][previous_location_index];
            }
            final_object_dataset[j].push_back(current_object);
        }
    }



    string dataset_output = "dataset/knn_dataset/" + dir + "/location/"  + map + "_" + to_string(int(k_num*100)) +"_" + num_rects + "_knn.csv";
    std::ofstream q_file(dataset_output);
    cout << dataset_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < final_object_dataset.size(); ++i) {
        for(int j = 0; j < final_object_dataset[i].size(); ++j)
            q_file << std::fixed << setprecision(8)  << final_object_dataset[i][j].id << " " << final_object_dataset[i][j].timestamp << " " << final_object_dataset[i][j].previous_location.x << " " << final_object_dataset[i][j].previous_location.y
                   << " " << final_object_dataset[i][j].new_location.x << " " << final_object_dataset[i][j].new_location.y <<"\n";
    }
    q_file.close();

}
//Each map has two location files: object knn location and query location
//Each map for each object size has two files: object keyword and query keyword
void generate_query_set(string dir, string map, int k){
    cout << "Generating knn dataset for " << map << endl;
    string mesh_path = "dataset/merged-mesh/" + dir + "/" + map + "-merged.mesh";
    ifstream meshfile(mesh_path);
    mp = new pl::Mesh(meshfile);
    mp->pre_compute_obstacle_edge_on_vertex();
    string grid_path = "dataset/grid/" + dir + "/" + map + ".map";
    mp->mark_turning_point(grid_path.c_str());

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
    string grid_label_path = "dataset/ebhl/" + dir + "/" + map + ".nt_adj_"+ "1";
    ebhl->load_adjacent_list(grid_label_path.c_str());
    string triangles_path = "dataset/ebhl/" + dir + "/" + map + ".nt_triangles_"+ "1";
    ebhl->load_non_taut_triangles(triangles_path.c_str());

    ebhlQueryV2 = new pl::EBHL_query_v2(mp, ebhl, turning_point, obstacle_middle_point);
    string label = "dataset/hub_label/" + dir + "/" + map + ".label";
    string order = "dataset/hub_label/" + dir + "/" + map + ".order";
    ebhlQueryV2 ->load_hub_label_on_convext_vertex(label, order);

    vector<Point> query_points;
    for(int i = 0; i < k; ++i){
        generate_random_point(query_points);
    }



    string dataset_output = "dataset/knn_dataset/" + dir + "/location/"  + map  + "_query.csv";
    std::ofstream q_file(dataset_output);
    cout << dataset_output << endl;
    //q_file<<"ID,Timestamp,Previous Location,New Location\n";
    for (int i = 0; i < query_points.size(); ++i) {
        q_file << std::fixed << setprecision(8)  << query_points[i].x << " " << query_points[i].y << "\n";
    }
    q_file.close();

//    vector<string> keyword_list {"strength","dexterity","constitution","intelligence","wisdom","charisma","acrobatics",
//                                 "animal","arcana","athletics","deception","history","intimidation","insight","investigation",
//                                 "medicine","nature","perception","performance","persuasion","religion","sleight","stealth","survival"};
//
//    vector<vector<string>> query_keywords(k);
//
//    for(int i = 0; i < k; ++i) {
//        unordered_map<string, int> keywords;
//        while (keywords.size() < 2) {
//            int index = rand() % keyword_list.size();
//            if (keywords.find(keyword_list[index]) == keywords.end()) {
//                keywords[keyword_list[index]] = 1;
//            }
//        }
//        for(auto & it: keywords){
//            query_keywords[i].emplace_back(it.first);
//        }
//    }
//
//    string query_output = "dataset/knn_dataset/" + dir + "/"  + map  + "_query_keyword.csv";
//    std::ofstream k_file(query_output);
//    cout << query_output << endl;
//    for (int i = 0; i < query_keywords.size(); ++i) {
//        for(int j = 0; j < query_keywords[i].size(); ++j){
//            k_file << fixed << setprecision(8) << i << " " << query_keywords[i][j] << endl;
//        }
//    }
//    k_file.close();


}

int main(int argc, char*argv[]){
    try{
        string dir;
        string map;
        string k;

        string ratio;
        string num_rects;
        if(argc != 6){
            cerr << argv[0] << "Directory | Map | K number" << endl;
            return 1;
        }
        else{
            dir = argv[1];
            map = argv[2];
            k = argv[3];

            ratio = argv[4];
            num_rects = argv[5];
        }
        //generate_knn_dataset(dir, map, k);

        generate_rectangle_object_sets(dir,map,k,ratio,num_rects);
        //generate_knn_dataset_insertion(dir, map, k);

        //generates a fixed number of query number
        //generate_query_set(dir, map, 5000);
    }catch(exception&err){
        cerr << "Stopped on exception lol: " << err.what() << endl;
    }
};
