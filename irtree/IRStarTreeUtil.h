#pragma once
#include "IRStarTree.h"
#include "KwData.h"

namespace irstar {
/*******************************************************************************
 * MinHeapEntry
 ******************************************************************************/
typedef struct MinHeapEntry
{
    MinHeapEntry(double key, Node_P nodePtr)
        : key(key), nodePtr(nodePtr), entryPtr(NULL)
        { }

    MinHeapEntry(double key, Entry_P entryPtr)
        : key(key), nodePtr(NULL), entryPtr(entryPtr)
        { }

    double key;
    Node_P nodePtr;
    Entry_P entryPtr;

    bool operator< (const MinHeapEntry& e) const;
}MinHeapEntry;

inline bool MinHeapEntry::operator< (const MinHeapEntry& e) const
{
    return key > e.key;
}

typedef vector<MinHeapEntry> MinHEntry_V;

/*******************************************************************************
 * MinHeap
 ******************************************************************************/
typedef struct MinHeap
{
    MinHEntry_V hv;

    void push(MinHeapEntry e);
    MinHeapEntry pop();
    bool isEmpty() { return hv.empty(); }
    void clear() { hv.clear(); }
}MinHeap;

inline void MinHeap::push(MinHeapEntry e)
{
    hv.push_back(e);
    push_heap(hv.begin(), hv.end());
}

inline MinHeapEntry MinHeap::pop()
{
    MinHeapEntry e = hv.front();
    pop_heap(hv.begin(), hv.end());
    hv.pop_back();
    return e;
}

/*******************************************************************************
 * MaxHeapEntry
 ******************************************************************************/
typedef struct MaxHeapEntry
{
    MaxHeapEntry(double key, Node_P nodePtr)
        : key(key), nodePtr(nodePtr), entryPtr(NULL)
        { }

    MaxHeapEntry(double key, Entry_P entryPtr)
        : key(key), nodePtr(NULL), entryPtr(entryPtr)
        { }

    double key;
    Node_P nodePtr;
    Entry_P entryPtr;

    bool operator< (const MaxHeapEntry& e) const;
}MaxHeapEntry;

inline bool MaxHeapEntry::operator< (const MaxHeapEntry& e) const
{
    return key < e.key;
}

typedef vector<MaxHeapEntry> MaxHEntry_V;

/*******************************************************************************
 * MaxHeap
 ******************************************************************************/
typedef struct MaxHeap
{
    MaxHEntry_V hv;

    void push(MaxHeapEntry e);
    MaxHeapEntry pop();
    bool isEmpty() { return hv.empty(); }
    void clear() { hv.clear(); }
}MaxHeap;

inline void MaxHeap::push(MaxHeapEntry e)
{
    hv.push_back(e);
    push_heap(hv.begin(), hv.end());
}

inline MaxHeapEntry MaxHeap::pop()
{
    MaxHeapEntry e = hv.front();
    pop_heap(hv.begin(), hv.end());
    hv.pop_back();
    return e;
}

/*******************************************************************************
 * RStarTreeUtil
 ******************************************************************************/
class RStarTreeUtil
{
    public:
        RStarTreeUtil() { }

        static bool validate(RStarTree& tree);

        static size_t rangeQuery2(RStarTree& tree, Point& center, double r, size_t k); // Exclusive outer range query
        static size_t rangeQuery2SmallTree(RStarTree& tree, Point& center, double r, size_t k); // Exclusive outer range query

        // range query in ring(minr, maxr, q)
        static void rangeQuery(RStarTree* tree, Point q, double minr, double maxr, std::vector<Data_P>& outIter);
        // incremental nearest neighbor retrieval
        static MinHeapEntry iNearestNeighbour(MinHeap& heap, Point q);
        static MinHeapEntry iNearestNeighbourKw(MinHeap& heap, Point q, const vector<string>& kw={"default"});
        static bool find(RStarTree& tree, Coord point[DIM]);
        static bool isEnclosed(Point& point, Node_P node); // test if point is inclusive enclosed by node
        static double minDis2(const Point& point, const Mbr& mbr); // min distance from point to node
        static double maxDis2(const Point& point, const Mbr& mbr); // max distance from point to node
        static double dis2(Point& point, Mbr& mbr); // max distance from point to point
};

inline bool RStarTreeUtil::isEnclosed(Point& point, Node_P node)
{
    Mbr& mbr = node->mbrn;
    for(size_t dim = 0; dim < DIM; dim++)
        if(mbr.coord[dim][0] > point.coord[dim] || mbr.coord[dim][1] < point.coord[dim])
            return false;
    return true;
}

inline double RStarTreeUtil::minDis2(const Point& point, const Mbr& mbr)
{
    double dis = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        if(mbr.coord[dim][0] > point.coord[dim] || mbr.coord[dim][1] < point.coord[dim])
        {
            double diff = min(abs(mbr.coord[dim][0] - point.coord[dim]), abs(mbr.coord[dim][1] - point.coord[dim]));
            dis += diff * diff;
        }
    return dis;
}

inline double RStarTreeUtil::maxDis2(const Point& point, const Mbr& mbr)
{
    double dis = 0;
    for(size_t dim = 0; dim < DIM; dim++)
    {
        double diff = max(abs(mbr.coord[dim][0] - point.coord[dim]), abs(mbr.coord[dim][1] - point.coord[dim]));
        dis += diff * diff;
    }
    return dis;
}

inline double RStarTreeUtil::dis2(Point& point, Mbr& mbr)
{
    double dis = 0;
    for(size_t dim = 0; dim < DIM; dim++)
    {
        double diff = mbr.coord[dim][0] - point.coord[dim];
        dis += diff * diff;
    }
    return dis;
}

}
