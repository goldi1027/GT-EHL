#pragma once
#include <algorithm>
#include <queue>
#include <unordered_map>
#include "irUtil2D.h"

namespace irstar {

#define MAX_CHILD (pageSize - sizeof(RTreeNode)) / sizeof(RTreeNode)

const size_t DIM = 2;
const size_t LEAF_LEVEL = 0;
const size_t PAGESIZE = 4096;

const double SPLIT_FACTOR = 0.4;
const double REINSERT_FACTOR = 0.3;
const size_t NEAR_MINIMUM_OVERLAP_FACTOR = 32;

/*******************************************************************************
 * Mbr - minimum bounding rectangle
 ******************************************************************************/

typedef struct Mbr
{
    Coord coord[DIM][2];

    Mbr();
    Mbr(Coord c[DIM]);
    Mbr(Coord c[DIM][2]);
    Mbr(Coord xmin, Coord xmax, Coord ymin, Coord ymax);

    void print();

    static Mbr getMbr(Mbr &mbr1, Mbr &mbr2);
    static double getOverlap(Mbr &mbr1, Mbr &mbr2);
    static double getEnlarge(Mbr &base, Mbr &add);

    double getMargin() const;
    double getArea() const;
    inline double getCenter(size_t dim);

    void enlarge(Mbr &add);
    void init();
}Mbr;

struct kwcnt{
  unordered_map<string, int> cnt;

  kwcnt() {
    cnt = unordered_map<string, int>();
  }
  kwcnt(unordered_map<string, int> kw) {
    cnt = unordered_map<string, int>(kw);
  }

  void reset() {
    cnt.clear();
  }

  kwcnt operator-(const kwcnt& other) const {
    kwcnt res = kwcnt{this->cnt};
    for (auto& it: other.cnt) {
      assert(res.cnt.find(it.first) != res.cnt.end());
      assert(res.cnt[it.first] >= it.second);
      res.cnt[it.first] -= it.second;
    }
    for (auto it=res.cnt.begin(); it != res.cnt.end(); ) {
      if (it->second == 0) it = res.cnt.erase(it);
      else it++;
    }
    return res;
  }

  void remove(const kwcnt& other) {
    for (auto it: other.cnt) {
      if (cnt.find(it.first) != cnt.end() && cnt[it.first] > 0) {
        assert(cnt[it.first] >= it.second);
        cnt[it.first] -= it.second;
      }
    }
    for (auto it=cnt.begin(); it != cnt.end();) {
      if (it->second == 0) it = cnt.erase(it);
      else it++;
    }
  }

  void update(const kwcnt& other) {
    for (auto it: other.cnt) {
      if (cnt.find(it.first) == cnt.end()) cnt[it.first] = it.second;
      else cnt[it.first] += it.second;
    }
  }

  bool equal(const kwcnt& other) {
    if (cnt.size() != other.cnt.size()) return false;
    for (const auto& it: other.cnt) {
      assert(it.second > 0);
      if (cnt.find(it.first) == cnt.end()) return false;
      if (cnt[it.first] != it.second) return false;
    }
    return true;
  }

  bool equal(const vector<string>& other) {
    if (cnt.size() != other.size()) return false;
    for (const auto& it: other) {
      if (cnt.find(it) == cnt.end()) return false;
      assert(cnt.find(it)->second > 0);
    }
    return true;
  }

  bool has(const vector<string>& kws) const {
    for (auto &it: kws) {
      auto k = cnt.find(it);
      assert(k == cnt.end() || k->second > 0);
      if (k == cnt.end()) return false;
    }
    return true;
  }

  bool check() {
    for (auto& it: cnt) if (it.second <= 0) return false;
    return true;
  }
};

/*******************************************************************************
 * RTreeNode
 ******************************************************************************/
struct RTreeNode;
struct LeafNodeEntry;
typedef void* Data_P;
typedef RTreeNode* Node_P;
typedef vector<Node_P> Node_P_V;
typedef list<Node_P> Node_P_L;
typedef LeafNodeEntry* Entry_P;
typedef vector<Entry_P> Entry_P_V;

struct LeafNodeEntry
{
    Mbr mbre;
    Data_P data;
    kwcnt kw;
    double value;
    Node_P parent = nullptr;

    LeafNodeEntry(Mbr mbre, Data_P data, unordered_map<string, int> kw);
    LeafNodeEntry(Mbr mbre, Data_P data);

    void print(string ident);
};

struct RTreeNode
{
    size_t level; // Level of the node. All leaf nodes are at level 0.
    size_t aggregate;
    Mbr mbrn;
    Node_P_V* children;
    Entry_P_V* entries;
    Node_P parent;
    double value;
    kwcnt kw;

    RTreeNode(size_t level);
    ~RTreeNode();

    bool operator < (const RTreeNode& node) const;

    void print(string ident);

    size_t size();
    void insert(Node_P childPtr);
    void insert(Entry_P entryPtr);
    void take(RTreeNode &source, size_t from); // Copy and remove elements from source.
    void resize(size_t from);
    void adjustMbr();
    void enlargeMbr(Mbr &mbr);
    void decreaseAggregate(size_t amount, kwcnt cnt);
    void increaseAggregate(size_t amount, kwcnt cnt);
    void prepareEnlargementSort(Mbr &add);
    void prepareDisSort();

    inline size_t entryNum() {
      size_t res = 0;
      if (level) {
        for (const auto it: *children) {
          res += it->entryNum();
        }
      }
      else {
        res += entries->size();
      }
      return res;
    }

    inline Data_P checkDuplicate() {
      queue<Node_P> q;
      q.push(this);
      set<Data_P> vis;
      Data_P res = nullptr;
      while (!q.empty()) {
        Node_P c = q.front(); q.pop();
        assert(c != nullptr);
        if (c->level) {
          for (const auto& it: *c->children) 
            q.push(it);
        }
        else {
          for (const auto& it: *c->entries) {
            if (vis.find(it->data) == vis.end())
              vis.insert(it->data);
            else {
              res = it->data;

              cout << "duplicate: ";
              it->mbre.print();
              if (c->parent != nullptr)
                c->parent->print("xxx");
              else
                c->print("xxx");
            }
          }
        }
      }
      return res;
    }

    inline bool validateAgg() {
      assert(this->kw.check());
      if (this->level) { // a node
        for (const auto it: *children) {
          if (!it->validateAgg()) return false;
        }
        kwcnt cnt;
        size_t num = 0;
        for (const auto it: *this->children) {
          cnt.update(it->kw);
          num += it->aggregate;
        }
        if (!cnt.equal(this->kw)) {
          cerr << "keywords cnt not agree with subtrees" << endl;
          assert(false);
          return false;
        }
        if (num != this->aggregate) {
          cerr << "aggregate not agree with subtrees" << endl;
          assert(false);
          return false;
        }
        return true;
      }
      else { // a leaf node
        kwcnt cnt;
        size_t num = this->entries->size();
        for (const auto it: *this->entries) {
          cnt.update(it->kw);
        }
        if (!cnt.equal(this->kw)) {
          cerr << "keywords cnt not agree with subtrees" << endl;
          assert(false);
          return false;
        }
        if (num != this->aggregate) {
          cerr << "aggregate not agree with subtrees" << endl;
          assert(false);
          return false;
        }
        return true;
      }
    }

};

/*******************************************************************************
 * Sorting
 ******************************************************************************/

typedef struct AxisSort
{
    size_t axis;

    AxisSort(size_t axis);

    bool operator() (Node_P child1, Node_P child2);
    bool operator() (Entry_P entry1, Entry_P entry2);
}AxisSort;

/*******************************************************************************
 * RStarTree
 ******************************************************************************/

typedef struct RStarTree
{
    static size_t pageSize; // Page size in bytes.
    static size_t maxChild; // M -- node capacity -- maximum number of entries
    static size_t minChild; // m -- minimum number of entries in a node

    RTreeNode* root;

    RStarTree();
    ~RStarTree();

    void print();

    /// Insertion
    void insertData(Entry_P entryPtr);

    void update(Entry_P entryPtr, Mbr cur, Mbr after);

    void delete_entry(Entry_P entryPtr);

    private:
        /// Insertion
        void insert(Node_P childPtr, Entry_P entryPtr, size_t desiredLevel, size_t &overflowLevel); // with OverflowTreatment

        /// Insertion - ChooseSubtree
        Node_P chooseSubtree(Mbr &mbr, size_t desiredLevel);

        /// Insertion - OverflowTreatment - ReInsert
        void reInsert(RTreeNode &node, size_t &overflowLevel);

        /// Insertion - OverflowTreatment - Split
        void split(RTreeNode &node);

        size_t chooseSplitAxis(RTreeNode &node);
        double computeS(RTreeNode &node);
        size_t chooseSplitIndex(RTreeNode &node, size_t axis);
}RStarTree;

}


