#include "intervaHeuristic.h"
#include "geometry.h"
#include "searchnode.h"
#include "successor.h"
#include "vertex.h"
#include "mesh.h"
#include "point.h"
#include "consts.h"
#include <queue>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <ctime>

namespace polyanya {

int IntervalHeuristic::succ_to_node(
    SearchNodePtr parent, Successor* successors, int num_succ,
    SearchNodePtr nodes) {
  // copy from searchinstance.cpp
  assert(mesh != nullptr);
  const Polygon& polygon = mesh->mesh_polygons[parent->next_polygon];
  const std::vector<int>& V = polygon.vertices;
  const std::vector<int>& P = polygon.polygons;

  double right_g = -1, left_g = -1;
  int out = 0;
  for (int i = 0; i < num_succ; i++) {
    const Successor& succ = successors[i];
    const int next_polygon = P[succ.poly_left_ind];
    if (next_polygon == -1) continue;
    // no end_polygon in knn
    if (mesh->mesh_polygons[next_polygon].is_one_way && end_polygons[next_polygon].empty()) continue;
    const int left_vertex = V[succ.poly_left_ind];
    const int right_vertex = succ.poly_left_ind? V[succ.poly_left_ind - 1]: V.back();

    // We implicitly set h to be zero and let search() update it.
    const auto process = [&](const int root, const double g) {
      if (root != -1) {
        assert(root >= 0 && root < (int) root_g_values.size());
        if (root_search_ids[root] != search_id) {
          // first time reaching root
          root_search_ids[root] = search_id;
          root_g_values[root] = g;
        } else {
          // visited
          if (root_g_values[root] + EPSILON < g) return;
          else root_g_values[root] = g;
        }
      }
      nodes[out++] = {nullptr, root, succ.left, succ.right, left_vertex, right_vertex, next_polygon, g, g};
    };

    const Point& parent_root = (parent->root == -1? start: mesh->mesh_vertices[parent->root].p);
    #define get_g(new_root) parent->g + parent_root.distance(new_root)

    switch (succ.type) {
      case Successor::RIGHT_NON_OBSERVABLE:
        if (right_g == -1) right_g = get_g(parent->right);
        process(parent->right_vertex, right_g);
        break;

      case Successor::OBSERVABLE:
        process(parent->root, parent->g);
        break;

      case Successor::LEFT_NON_OBSERVABLE:
        if (left_g == -1) left_g = get_g(parent->left);
        process(parent->left_vertex, left_g);
        break;

      default:
        assert(false);
        break;
    }
  }
  #undef get_g
  return out;
}

void IntervalHeuristic::set_end_polygon() {
  end_polygons.resize(mesh->mesh_polygons.size());
  goal_polyid.resize(goals.size());
  fill(goal_polyid.begin(), goal_polyid.end(), -1);
  for (int i=0; i<(int)mesh->mesh_polygons.size(); i++) end_polygons[i].clear();
  for (int i=0; i<(int)goals.size(); i++) {
    int poly_id = get_point_location_in_search(goals[i], mesh, verbose).poly1;
    if (poly_id == -1) continue;
    assert(poly_id < (int)end_polygons.size());
    end_polygons[poly_id].push_back(i);
    goal_polyid[i] = poly_id;
  }
}


void IntervalHeuristic::gen_initial_nodes() {
  // revised from SearchInstance::gen_initial_nodes
  // modify:
  // 1. h value for get_lazy() is 0
  // 2. no end_polygon in knn search
  const PointLocation pl = get_point_location_in_search(start, mesh, verbose);
  #define get_lazy(next, left, right) new (node_pool->allocate()) SearchNode \
    {nullptr, -1, start, start, left, right, next, 0, 0}
  #define v(vertex) mesh->mesh_vertices[vertex]

  const auto push_lazy = [&](SearchNodePtr lazy) {
    const int poly = lazy->next_polygon;
    if (poly == -1) return;

    if (!end_polygons[poly].empty()) {
      for (int gid: end_polygons[poly]) {
        const Point& goal = goals[gid];
        SearchNodePtr final_node = get_lazy(lazy->next_polygon, lazy->left_vertex, lazy->right_vertex);
        final_node->f += start.distance(goal);
        final_node->set_reached();
        final_node->set_goal_id(gid);

        #ifndef NDEBUG
        if (verbose) {
          std::cerr << "generating init node: ";
          print_node(final_node, std::cerr);
          std::cerr << std::endl;
        }
        #endif
        open_list.push(final_node);
        nodes_generated++;
        nodes_pushed++;
      }
    }

    const std::vector<int>& vertices = mesh->mesh_polygons[poly].vertices;
    Successor* successors = new Successor [vertices.size()];
    int last_vertex = vertices.back();
    int num_succ = 0;
    for (int i = 0; i < (int) vertices.size(); i++) {
      const int vertex = vertices[i];
      if (vertex == lazy->right_vertex ||
          last_vertex == lazy->left_vertex) {
        last_vertex = vertex;
        continue;
      }
      successors[num_succ++] = {Successor::OBSERVABLE, v(vertex).p, v(last_vertex).p, i};
      last_vertex = vertex;
    }
    SearchNode* nodes = new SearchNode [num_succ];
    const int num_nodes = succ_to_node(lazy, successors, num_succ, nodes);

    delete[] successors;
    for (int i = 0; i < num_nodes; i++) {
      SearchNodePtr nxt = new (node_pool->allocate()) SearchNode(nodes[i]);
      const Point& nxt_root = (nxt->root == -1? start: mesh->mesh_vertices[nxt->root].p);
			if (!this->isZero)
				nxt->f += get_interval_heuristic(nxt_root, nxt->left, nxt->right);
      nxt->parent = lazy;
      #ifndef NDEBUG
      if (verbose) {
        std::cerr << "generating init node: ";
        print_node(nxt, std::cerr);
        std::cerr << std::endl;
      }
      #endif
      open_list.push(nxt);
      nodes_pushed++;

      if (!end_polygons[nxt->next_polygon].empty()) {
        gen_final_nodes(nxt, nxt_root);
      }
    }
    delete[] nodes;
    nodes_generated += num_nodes;
    nodes_pushed += num_nodes;
  };
  switch (pl.type) {
    case PointLocation::NOT_ON_MESH:
      break;
    case PointLocation::ON_CORNER_VERTEX_AMBIG:
      if (pl.poly1 == -1) break;
    case PointLocation::ON_CORNER_VERTEX_UNAMBIG:
    case PointLocation::IN_POLYGON:
      {
        SearchNodePtr lazy = get_lazy(pl.poly1, -1, -1);
        push_lazy(lazy);
        nodes_generated++;
      }
      break;
    case PointLocation::ON_MESH_BORDER:
      {
        SearchNodePtr lazy = get_lazy(pl.poly1, -1, -1);
        push_lazy(lazy);
        nodes_generated++;
      }
      break;
    case PointLocation::ON_EDGE:
      {
        SearchNodePtr lazy1 = get_lazy(pl.poly2, pl.vertex1, pl.vertex2);
        SearchNodePtr lazy2 = get_lazy(pl.poly1, pl.vertex2, pl.vertex1);
        push_lazy(lazy1);
        nodes_generated++;
        push_lazy(lazy2);
        nodes_generated++;
      }
      break;
    case PointLocation::ON_NON_CORNER_VERTEX:
      {
        for (int& poly: v(pl.vertex1).polygons) {
          SearchNodePtr lazy = get_lazy(poly, pl.vertex1, pl.vertex1);
          push_lazy(lazy);
          nodes_generated++;
        }
      }
      break;

    default:
      assert(false);
      break;
  }
  #undef v
  #undef get_lazy
}

#define root_to_point(root) ((root) == -1 ? start : mesh->mesh_vertices[root].p)

int IntervalHeuristic::search() {
  init_search();
  timer.start();
  if (mesh == nullptr) {
    timer.stop();
    return 0;
  }

  while (!open_list.empty()) {
    SearchNodePtr node = open_list.top(); open_list.pop();

    #ifndef NDEBUG
    if (verbose) {
      std::cerr << "popped off: ";
      print_node(node, std::cerr);
      std::cerr << std::endl;
    }
    #endif

    nodes_popped++;
    if (node->reached) {
      deal_final_node(node);
      if ((int)final_nodes.size() == K) break;
      continue;
    }
    const int root = node->root;
    if (root != -1) {
      assert(root < (int) root_g_values.size());
      if (root_search_ids[root] == search_id) {
        // We've been here before!
        // Check whether we've done better.
        if (root_g_values[root] + EPSILON < node->g) {
          nodes_pruned_post_pop++;

          #ifndef NDEBUG
          if (verbose) std::cerr << "node is dominated!" << std::endl;
          #endif
          continue;
        }
      }
    }
    int num_nodes = 1;
    search_nodes_to_push[0] = *node;
    // Intermediate pruning
    do {
      SearchNode cur = search_nodes_to_push[0];
      int num_succ = get_successors(cur, start, *mesh, search_successors);
      successor_calls++;
      num_nodes = succ_to_node(&cur, search_successors, num_succ, search_nodes_to_push);
      //break;
      if (num_nodes == 1) { // we should continue
        // Did we turn?
        if (cur.g != search_nodes_to_push[0].g) {
          // Turned. Set the parent of this, and set the current
          // node pointer to this after allocating space for it.
          search_nodes_to_push[0].parent = node;
          node = new (node_pool->allocate()) SearchNode(search_nodes_to_push[0]);
          nodes_generated++;
        }
        if (!end_polygons[search_nodes_to_push[0].next_polygon].empty()) {
          const SearchNodePtr nxt = new (node_pool->allocate()) SearchNode(search_nodes_to_push[0]);
          nxt->parent = node;
          const Point& nxt_root = nxt->root == -1? start: mesh->mesh_vertices[nxt->root].p;
          gen_final_nodes(nxt, nxt_root);
        }
        #ifndef NDEBUG
        if (verbose) {
          std::cerr << "\tintermediate: ";
          print_node(&search_nodes_to_push[0], std::cerr);
          std::cerr << std::endl;
        }
        #endif
      }
    } while (num_nodes == 1); // if num_nodes == 0, we still break

    for (int i = 0; i < num_nodes; i++) {
      // update h value before we push
      const SearchNodePtr nxt = new (node_pool->allocate()) SearchNode(search_nodes_to_push[i]);
      const Point& nxt_root = (nxt->root == -1 ? start: mesh->mesh_vertices[nxt->root].p);
			if (!this->isZero)
				nxt->f += get_interval_heuristic(nxt_root, nxt->left, nxt->right);
      nxt->parent = node;
      #ifndef NDEBUG
      if (verbose) {
        std::cerr << "\tpushing: ";
        print_node(nxt, std::cerr);
        std::cerr << std::endl;
      }
      #endif
      open_list.push(nxt);
      nodes_pushed++;
      nodes_generated++;

      // when nxt can be final_node
      int nxt_poly = nxt->next_polygon;
      if (!end_polygons[nxt_poly].empty()) {
        gen_final_nodes(nxt, nxt_root);
      }
    }
  }
  timer.stop();
  return (int)final_nodes.size();
}

int IntervalHeuristic::keyword_search() {
        init_search();
        timer.start();
        if (mesh == nullptr) {
            timer.stop();
            return 0;
        }

        while (!open_list.empty()) {
            SearchNodePtr node = open_list.top(); open_list.pop();

#ifndef NDEBUG
            if (verbose) {
                std::cerr << "popped off: ";
                print_node(node, std::cerr);
                std::cerr << std::endl;
            }
#endif

            nodes_popped++;
            if (node->reached) {
//                if((object_keyword[node->goal_id][query_keyword[start_id][0]] > 0) && (object_keyword[node->goal_id][query_keyword[start_id][1]] > 0)){
//                    deal_final_node(node);
//                    if ((int)final_nodes.size() == K) break;
//                    continue;
//                }
                query_counter = 0;
                for(int i = 0; i < query_keyword[start_id].size(); ++i){
                    if(object_keyword[node->goal_id].find(query_keyword[start_id][i]) != object_keyword[node->goal_id].end()){
                        query_counter++;
                    }
                }
                if(query_counter == query_keyword[start_id].size()){
                    deal_final_node(node);
                    if ((int)final_nodes.size() == K) break;
                    continue;
                }
            }
            const int root = node->root;
            if (root != -1) {
                assert(root < (int) root_g_values.size());
                if (root_search_ids[root] == search_id) {
                    // We've been here before!
                    // Check whether we've done better.
                    if (root_g_values[root] + EPSILON < node->g) {
                        nodes_pruned_post_pop++;

#ifndef NDEBUG
                        if (verbose) std::cerr << "node is dominated!" << std::endl;
#endif
                        continue;
                    }
                }
            }
            int num_nodes = 1;
            search_nodes_to_push[0] = *node;
            // Intermediate pruning
            do {
                SearchNode cur = search_nodes_to_push[0];
                int num_succ = get_successors(cur, start, *mesh, search_successors);
                successor_calls++;
                num_nodes = succ_to_node(&cur, search_successors, num_succ, search_nodes_to_push);
                //break;
                if (num_nodes == 1) { // we should continue
                    // Did we turn?
                    if (cur.g != search_nodes_to_push[0].g) {
                        // Turned. Set the parent of this, and set the current
                        // node pointer to this after allocating space for it.
                        search_nodes_to_push[0].parent = node;
                        node = new (node_pool->allocate()) SearchNode(search_nodes_to_push[0]);
                        nodes_generated++;
                    }
                    if (!end_polygons[search_nodes_to_push[0].next_polygon].empty()) {
                        const SearchNodePtr nxt = new (node_pool->allocate()) SearchNode(search_nodes_to_push[0]);
                        nxt->parent = node;
                        const Point& nxt_root = nxt->root == -1? start: mesh->mesh_vertices[nxt->root].p;
                        gen_keyword_final_nodes(nxt, nxt_root);
                    }
#ifndef NDEBUG
                    if (verbose) {
                        std::cerr << "\tintermediate: ";
                        print_node(&search_nodes_to_push[0], std::cerr);
                        std::cerr << std::endl;
                    }
#endif
                }
            } while (num_nodes == 1); // if num_nodes == 0, we still break

            for (int i = 0; i < num_nodes; i++) {
                // update h value before we push
                const SearchNodePtr nxt = new (node_pool->allocate()) SearchNode(search_nodes_to_push[i]);
                const Point& nxt_root = (nxt->root == -1 ? start: mesh->mesh_vertices[nxt->root].p);
                if (!this->isZero)
                    nxt->f += get_interval_heuristic(nxt_root, nxt->left, nxt->right);
                nxt->parent = node;
#ifndef NDEBUG
                if (verbose) {
                    std::cerr << "\tpushing: ";
                    print_node(nxt, std::cerr);
                    std::cerr << std::endl;
                }
#endif
                open_list.push(nxt);
                nodes_pushed++;
                nodes_generated++;

                // when nxt can be final_node
                int nxt_poly = nxt->next_polygon;
                if (!end_polygons[nxt_poly].empty()) {
                    gen_keyword_final_nodes(nxt, nxt_root);
                }
            }
        }
        timer.stop();
        return (int)final_nodes.size();
    }

void IntervalHeuristic::print_node(SearchNodePtr node, std::ostream& outfile) {
  outfile << "root=" << root_to_point(node->root) << "; left=" << node->left
          << "; right=" << node->right << "; f=" << node->f << ", g="
          << node->g;
}

void IntervalHeuristic::get_path_points(std::vector<Point>& out, int k) {
  if (k >= (int)goals.size()) return;
  assert((int)final_nodes.size() <= K);
  assert(final_nodes[k]->goal_id != -1);
  assert(final_nodes[k]->reached == true);
  out.clear();
  out.push_back(goals[final_nodes[k]->goal_id]);
  SearchNodePtr cur = final_nodes[k];

  while (cur != nullptr) {
    if (root_to_point(cur->root) != out.back()) out.push_back(root_to_point(cur->root));
    cur = cur->parent;
  }
  std::reverse(out.begin(), out.end());
}

void IntervalHeuristic::get_goals(std::vector<Point> &out, int k){
    for(int i = 0; i < k; ++i){
        out.push_back(goals[i]);
    }
}

void IntervalHeuristic::print_search_nodes(std::ostream& outfile, int k) {
  if (k > (int)final_nodes.size()) return;
  SearchNodePtr cur = final_nodes[k];
  while (cur != nullptr) {
    print_node(cur, outfile);
    outfile << std::endl;
    mesh->print_polygon(outfile, cur->next_polygon);
    outfile << std::endl;
    cur = cur->parent;
  }
}

void IntervalHeuristic::deal_final_node(const SearchNodePtr node) {

  const Point& goal = goals[node->goal_id];
  const int final_root = [&]() {
      const Point& root = root_to_point(node->root);
      const Point root_goal = goal - root;
      // If root-left-goal is not CW, use left.
      if (root_goal * (node->left - root) < -EPSILON) {
          return node->left_vertex;
      }
      // If root-right-goal is not CCW, use right.
      if ((node->right - root) * root_goal < -EPSILON)
      {
          return node->right_vertex;
      }
      // Use the normal root.
      return node->root;
  }();

  assert(node->goal_id != -1);
  assert(reached.find(node->goal_id) == reached.end() || reached[node->goal_id] < node->f + EPSILON);

  if (reached.find(node->goal_id) == reached.end()) {
    int end_polygon = node->next_polygon;
    const SearchNodePtr true_final =
      new (node_pool->allocate()) SearchNode
      {node, final_root, goal, goal, -1, -1, end_polygon, node->f, node->g};
    true_final->set_reached();
    true_final->set_goal_id(node->goal_id);
    reached[node->goal_id] = node->f;
    final_nodes.push_back(true_final);
    nodes_generated++;
  }
}

void IntervalHeuristic::gen_final_nodes(const SearchNodePtr node, const Point& rootPoint) {
    assert(node->next_polygon != -1);
    for (int gid: end_polygons[node->next_polygon]) {
          const Point& goal = goals[gid];
          SearchNodePtr final_node = new (node_pool->allocate()) SearchNode(*node);
          final_node->set_reached();
          final_node->set_goal_id(gid);
          final_node->f = final_node->g + get_h_value(rootPoint, goal, node->left, node->right);
#ifndef NDEBUG
          if (verbose) {
              std::cerr << "\tpushing: ";
              print_node(final_node, std::cerr);
              std::cerr << std::endl;
          }
#endif
          open_list.push(final_node);
          nodes_generated++;
          nodes_pushed++;
      }
}

void IntervalHeuristic::gen_keyword_final_nodes(const SearchNodePtr node, const Point& rootPoint) {
        assert(node->next_polygon != -1);
        for (int gid: end_polygons[node->next_polygon]) {
            query_counter = 0;
            for(int i = 0; i < query_keyword[start_id].size(); ++i){
                if(object_keyword[gid].find(query_keyword[start_id][i]) != object_keyword[gid].end()){
                       query_counter++;
                    }
            }
            if(query_counter == query_keyword[start_id].size()){
                const Point& goal = goals[gid];
                SearchNodePtr final_node = new (node_pool->allocate()) SearchNode(*node);
                final_node->set_reached();
                final_node->set_goal_id(gid);
                final_node->f = final_node->g + get_h_value(rootPoint, goal, node->left, node->right);
#ifndef NDEBUG
                if (verbose) {
                    std::cerr << "\tpushing: ";
                    print_node(final_node, std::cerr);
                    std::cerr << std::endl;
                }
#endif
                open_list.push(final_node);
                nodes_generated++;
                nodes_pushed++;
            }
//            if((object_keyword[gid][query_keyword[start_id][0]] > 0) && (object_keyword[gid][query_keyword[start_id][1]] > 0)){
//                const Point& goal = goals[gid];
//                SearchNodePtr final_node = new (node_pool->allocate()) SearchNode(*node);
//                final_node->set_reached();
//                final_node->set_goal_id(gid);
//                final_node->f = final_node->g + get_h_value(rootPoint, goal, node->left, node->right);
//#ifndef NDEBUG
//                if (verbose) {
//                    std::cerr << "\tpushing: ";
//                    print_node(final_node, std::cerr);
//                    std::cerr << std::endl;
//                }
//#endif
//                open_list.push(final_node);
//                nodes_generated++;
//                nodes_pushed++;
//            }

        }
    }

#undef root_to_point

}
