#include "IRStarTree.h"

namespace irstar {
/*******************************************************************************
 * Mbr
 ******************************************************************************/

Mbr::Mbr()
{
    init();
}

Mbr::Mbr(Coord c[DIM])
{
    for(size_t dim = 0; dim < DIM; dim++)
    {
        coord[dim][0] = c[dim];
        coord[dim][1] = c[dim];
    }
}

Mbr::Mbr(Coord c[DIM][2])
{
    for(size_t dim = 0; dim < DIM; dim++)
    {
        coord[dim][0] = c[dim][0];
        coord[dim][1] = c[dim][1];
    }
}

Mbr::Mbr(Coord xmin, Coord xmax, Coord ymin, Coord ymax)
{
    coord[0][0] = xmin;
    coord[0][1] = xmax;
    coord[1][0] = ymin;
    coord[1][1] = ymax;
}

void Mbr::print()
{
    for(size_t dim = 0; dim < DIM; dim++)
        cout << "<" << coord[dim][0] << "," << coord[dim][1] << ">";
    cout << endl;
}

Mbr Mbr::getMbr(Mbr &mbr1, Mbr &mbr2)
{
    Mbr mbr;
    for(size_t dim = 0; dim < DIM; dim++)
    {
        mbr.coord[dim][0] = min(mbr1.coord[dim][0], mbr2.coord[dim][0]);
        mbr.coord[dim][1] = max(mbr1.coord[dim][1], mbr2.coord[dim][1]);
    }
    return mbr;
}

double Mbr::getOverlap(Mbr &mbr1, Mbr &mbr2)
{
    double overlap = 1;
    for(size_t dim = 0; dim < DIM; dim++)
    {
        Coord maxMin = max(mbr1.coord[dim][0], mbr2.coord[dim][0]);
        Coord minMax = min(mbr1.coord[dim][1], mbr2.coord[dim][1]);
        if(maxMin >= minMax)
            return 0;
        overlap *= minMax - maxMin;
    }
    return overlap;
}

double Mbr::getEnlarge(Mbr &base, Mbr &add)
{
    Mbr mbr = Mbr::getMbr(base, add);
    return mbr.getArea() - base.getArea();
}

double Mbr::getMargin() const
{
    double margin = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        margin += coord[dim][1] - coord[dim][0];
    return margin;
}

double Mbr::getArea() const
{
    double area = 1;
    for(size_t dim = 0; dim < DIM; dim++)
        area *= coord[dim][1] - coord[dim][0];
    return area;
}

inline double Mbr::getCenter(size_t dim)
{
    return (coord[dim][0] + coord[dim][1]) * 0.5;
}

void Mbr::enlarge(Mbr &add)
{
    for(size_t dim = 0; dim < DIM; dim++)
    {
        coord[dim][0] = min(coord[dim][0], add.coord[dim][0]);
        coord[dim][1] = max(coord[dim][1], add.coord[dim][1]);
    }
}

void Mbr::init()
{
    for(size_t dim = 0; dim < DIM; dim++)
    {
        coord[dim][0] = INF_P;
        coord[dim][1] = INF_N;
    }
}

/*******************************************************************************
 * RTreeNode
 ******************************************************************************/

LeafNodeEntry::LeafNodeEntry(Mbr mbre, Data_P data, unordered_map<string, int> kw)
    : mbre(mbre), data(data), kw(kw)
{
}

LeafNodeEntry::LeafNodeEntry(Mbr mbre, Data_P data)
    : mbre(mbre), data(data)
{
  kw = unordered_map<string, int>({{"default", 1}});
}

void LeafNodeEntry::print(string ident)
{
    cout << ident << " {" << *(int*)this->data << "}@" << (void*)this;
    mbre.print();
}

RTreeNode::RTreeNode(size_t level)
    : level(level), aggregate(0), children(NULL), entries(NULL), parent(NULL)
{
    if(level)
        children = new Node_P_V();
    else
        entries = new Entry_P_V();
}

RTreeNode::~RTreeNode()
{
    if(children)
    {
        for(size_t ic = 0; ic < children->size(); ic++)
            delete children->at(ic);
        delete children;
    }
    if(entries)
        delete entries;
}

bool RTreeNode::operator< (const RTreeNode& node) const
{
    return value < node.value;
}

void RTreeNode::print(string ident)
{
    cout << ident << level << " #" << aggregate << "@" << (void*)this;
    mbrn.print();
    if(level) // If this is not a leaf node.
        for(size_t ic = 0; ic < children->size(); ic++)
        {
            Node_P nodePtr = children->at(ic);
            nodePtr->print(ident + "| ");
        }
    else // If this is a leaf node.
        for(size_t ie = 0; ie < entries->size(); ie++)
        {
            Entry_P entryPtr = entries->at(ie);
            entryPtr->print(ident + "| ");
        }
}

size_t RTreeNode::size()
{
    return level ? children->size() : entries->size();
}

void RTreeNode::insert(Node_P childPtr)
{
    children->push_back(childPtr);
    childPtr->parent = this;
    enlargeMbr(childPtr->mbrn);
    increaseAggregate(childPtr->aggregate, childPtr->kw);
}

void RTreeNode::insert(Entry_P entryPtr)
{
    entryPtr->parent = this;
    entries->push_back(entryPtr);
    // cout << "before enlarge " << endl;
     //this->print(">");
    enlargeMbr(entryPtr->mbre);
    // cout << "after enlarge " << endl;
     //this->print("<");
    increaseAggregate(1, entryPtr->kw);
}

void RTreeNode::take(RTreeNode &source, size_t from)
{
    if(source.level) // If source is not a leaf node.
    {
        Node_P_V &sourceV = *source.children;
        for(size_t ic = from; ic < sourceV.size(); ic++)
            insert(sourceV.at(ic));
    }
    else // If source is a leaf node.
    {
        Entry_P_V &sourceV = *source.entries;
        for(size_t ie = from; ie < sourceV.size(); ie++)
            insert(sourceV.at(ie));
    }
    source.resize(from);
}

void RTreeNode::resize(size_t from)
{
    kwcnt cnt;
    if(level) // If this is not a leaf node.
    {
        children->resize(from);
        size_t newAggregate = 0;
        for(size_t ic = 0; ic < children->size(); ic++) {
            newAggregate += children->at(ic)->aggregate;
            cnt.update(children->at(ic)->kw);
        }
        decreaseAggregate(aggregate - newAggregate, kw - cnt);
    }
    else // If this is a leaf node.
    {
        for (size_t i=from; i<entries->size(); i++)
          cnt.update(entries->at(i)->kw);
        entries->resize(from);
        decreaseAggregate(aggregate - from, cnt);
    }
    //cout << "===========================" << endl;
    //this->print("(");
    adjustMbr();
}

void RTreeNode::adjustMbr()
{
    mbrn.init();
    if(level) // If this is not a leaf node.
        for(size_t ic = 0; ic < children->size(); ic++)
            mbrn.enlarge(children->at(ic)->mbrn);
    else // If this is a leaf node.
        for(size_t ie = 0; ie < entries->size(); ie++)
            mbrn.enlarge(entries->at(ie)->mbre);
    if(parent)
        parent->adjustMbr();
}

void RTreeNode::enlargeMbr(Mbr &mbr)
{
    mbrn.enlarge(mbr);
    if(parent)
        parent->enlargeMbr(mbr);
}

void RTreeNode::decreaseAggregate(size_t amount, kwcnt cnt)
{
    this->kw.remove(cnt);
    aggregate -= amount;
    // assert(validateAgg() == true);
    if(parent)
        parent->decreaseAggregate(amount, cnt);
    // assert(validateAgg() == true);
}

void RTreeNode::increaseAggregate(size_t amount, kwcnt cnt)
{
    this->kw.update(cnt);
    aggregate += amount;
    // assert(validateAgg() == true);
    // cout << "cur node: " << endl;
    //this->print("@");
    if(parent)
        parent->increaseAggregate(amount, cnt);
    // assert(validateAgg() == true);
}

void RTreeNode::prepareEnlargementSort(Mbr &add)
{
    for(size_t ic = 0; ic < children->size(); ic++)
    {
        Node_P childPtr = children->at(ic);
        childPtr->value = Mbr::getEnlarge(childPtr->mbrn, add);
    }
}

void RTreeNode::prepareDisSort()
{
    Coord center[DIM];
    for(size_t dim = 0; dim < DIM; dim++)
        center[dim] = mbrn.getCenter(dim);
    if(level)
        for(size_t ic = 0; ic < children->size(); ic++)
        {
            Node_P childPtr = children->at(ic);
            childPtr->value = 0;
            for(size_t dim = 0; dim < DIM; dim++)
            {
                double sideLength = childPtr->mbrn.getCenter(dim) - center[dim];
                childPtr->value += sideLength * sideLength;
            }
        }
    else
        for(size_t ie = 0; ie < entries->size(); ie++)
        {
            Entry_P entryPtr = entries->at(ie);
            entryPtr->value = 0;
            for(size_t dim = 0; dim < DIM; dim++)
            {
                double sideLength = entryPtr->mbre.getCenter(dim) - center[dim];
                entryPtr->value += sideLength * sideLength;
            }
        }
}

/*******************************************************************************
 * Sorting
 ******************************************************************************/

AxisSort::AxisSort(size_t axis)
    : axis(axis)
{
}

bool AxisSort::operator() (Node_P child1, Node_P child2)
{
    if(child1->mbrn.coord[axis][0] == child2->mbrn.coord[axis][0])
        return (child1->mbrn.coord[axis][1] < child2->mbrn.coord[axis][1]);
    return child1->mbrn.coord[axis][0] < child2->mbrn.coord[axis][0];
}

bool AxisSort::operator() (Entry_P entry1, Entry_P entry2)
{
    if(entry1->mbre.coord[axis][0] == entry2->mbre.coord[axis][0])
        return (entry1->mbre.coord[axis][1] < entry2->mbre.coord[axis][1]);
    return entry1->mbre.coord[axis][0] < entry2->mbre.coord[axis][0];
}

/*******************************************************************************
 * RTree
 ******************************************************************************/

size_t RStarTree::pageSize = PAGESIZE;
size_t RStarTree::maxChild = MAX_CHILD;
size_t RStarTree::minChild = (size_t)maxChild * SPLIT_FACTOR >=1?maxChild*SPLIT_FACTOR:1;

RStarTree::RStarTree()
    : root(NULL)
{
    root = new RTreeNode(LEAF_LEVEL);
}

RStarTree::~RStarTree()
{
    if(root)
        delete root;
}

void RStarTree::print()
{
    root->print("");
}

/** Invoke Insert starting with the leaf level as a parameter, to insert a new data rectangle.
  */
void RStarTree::insertData(Entry_P entryPtr)
{
    /// Bits of overflowLevel indicate first time overflow of the coresponding level.
    size_t overflowLevel = -1; // initialize 32 bits of 1s.
    overflowLevel <<= (root->level); // set the lower (root->level) bits to 0s.
    insert(NULL, entryPtr, LEAF_LEVEL, overflowLevel);
}

/** Invoke ChooseSubtree, with the level as a parameter, to find an appropriate node N, in which to place the new entry E.
  *
  * If N has less than M entries, accommodate E in N.
  * If N has M entries, invoke OverflowTreatment with the level of N as a parameter [for reinsertion or split].
  *
  * If OverflowTreatment was called and a plit was performed, propagate OverflowTreatment upwards if necessary.
  * If OverflowTreatment caused a split of the root, create a new root.
  *
  * Adjust all covering rectangles in the insertion path such that they are minimum bounding boxes enclosing their children rectangles.
  *
  * --------------------------------------------------
  *
  * [OverflowTreatment]
  * If the level is not the root level and this is the first call of OverflowTreatment in the given level during the insertion of one data rectangle, then invoke ReInsert.
  * Else invoke Split.
  */
void RStarTree::insert(Node_P childPtr, Entry_P entryPtr, size_t desiredLevel, size_t &overflowLevel)
{
    /// Bits of overflowLevel indicate first time overflow of the coresponding level.
    Mbr &mbr = (desiredLevel ? childPtr->mbrn : entryPtr->mbre);
    Node_P node = chooseSubtree(mbr, desiredLevel);
    // cout << "before insert" << endl;
    // root->print("..");
    if(desiredLevel) // If desiredLevel is not leaf node level
        node->insert(childPtr);
    else // If desiredLevel is leaf node level
        node->insert(entryPtr);

    // cout << "after insert" << endl;
    // root->print("..");

    while(node) // while node is not NULL
    {
        while(node->size() > RStarTree::maxChild)
        {
            if(overflowLevel & (1<<node->level))
                split(*node);
            else
                reInsert(*node, overflowLevel);
        }
        node = node->parent;
    }
}



void RStarTree::delete_entry(Entry_P entryPtr){
    Node_P nptr = entryPtr->parent;
    assert(nptr != nullptr);
    //nptr->print("(");
    Entry_P_V& entries = *(nptr->entries);
    bool flag = false;
    for (int i=0; i<(int)entries.size(); i++) {
        if (entries[i] == entryPtr) {
            entries[i] = entries.back();
            entries[entries.size()-1] = entryPtr;
            flag = true;
            break;
        }
    }
    assert(flag);
    entryPtr->parent = nullptr;
    if(entries.size() >=1){
        nptr->resize(entries.size() - 1);
    }

}

void RStarTree::update(Entry_P entryPtr, Mbr cur, Mbr after) {
  Node_P nptr = entryPtr->parent;
  assert(nptr != nullptr);
  Entry_P_V& entries = *(nptr->entries);
  // cout << "in update before swap" << endl;
  // nptr->print("(");
  // swap to back
  bool flag = false;
  for (int i=0; i<(int)entries.size(); i++) {
    if (entries[i] == entryPtr) {
      entries[i] = entries.back();
      entries[entries.size()-1] = entryPtr;
      flag = true;
      break;
    }
  }
  assert(flag);
  entryPtr->parent = nullptr;

  // cout << "in update after swap" << endl;
  // nptr->print(")");
  // remove
  if(entries.size() >=1){
      nptr->resize(entries.size() - 1);
  }
  // auto dup = root->checkDuplicate();
  // assert(dup == nullptr);
  // reinsert
  entryPtr->mbre = after;
  insertData(entryPtr);

  // dup = root->checkDuplicate();
  // assert(dup == nullptr);
}

/** Set N to be the root.
  *
  * If N is a leaf, return N.
  * If the childpointers in N point to leaves [determin the minimum overlap cost].
  *     Choose the entry in N whose rectangle needs least overlap enlargement to include the new data rectangle.
  *     Resolve ties by choosing the entry whose rectangle needs least area enlargement.
  *     Then the entry with the rectangle of smallest area.
  * If the childpointers in N do not point to leaves [determin the minimum area cost].
  *     Choose the entry in N whose rectangle needs least area enlargement to include the new data rectangle.
  *     Resolve ties by choosing the entry with the rectangle of smallest area.
  *
  * Set N to be the childnode pointed to by the childpointer of the chosen entry and repeat.
  *
  * --------------------------------------------------
  *
  * [determine the nearly minimum overlap cost]
  * Sort the rectangles in N in increasing order of their area enlargement needed to include the new data rectangle.
  * Let A be the group of the first p entries.
  * From the entries in A, considering all entries in N, choose the entry whose rectangle needs least overlap enlargement.
  * Resolve ties as described above.
  */
Node_P RStarTree::chooseSubtree(Mbr &mbr, size_t desiredLevel)
{
    Node_P node = root;
    while(node->level > desiredLevel)
    {
        /// since (node->level > desiredLevel), hence node->level > 0, node is not a leaf node
        Node_P_V &children = *node->children;
        node->prepareEnlargementSort(mbr);
        sort(children.begin(), children.end());
        //sort(children.begin(), children.end(), EnlargementSort(mbr));

        size_t selectedIndex = 0;
        if(node->level == 1) // if the child pointers point to leaf nodes
        {
            size_t p = min(node->size(), NEAR_MINIMUM_OVERLAP_FACTOR);
            double minOverlapEnlarge = INF_P;
            for(size_t ic = 0; ic < p; ic++)
            {
                Mbr base = children.at(ic)->mbrn;
                Mbr newMbr = Mbr::getMbr(base, mbr);
                double overlapBase = 0, overlap = 0;
                for(size_t ico = 0; ico < node->size(); ico++)
                    if(ico != ic)
                    {
                        overlapBase += Mbr::getOverlap(base, children.at(ico)->mbrn);
                        overlap += Mbr::getOverlap(newMbr, children.at(ico)->mbrn);
                    }
                double overlapEnlarge = overlap - overlapBase;
                if(overlapEnlarge - EPS < minOverlapEnlarge)
                {
                    minOverlapEnlarge = overlapEnlarge;
                    selectedIndex = ic;
                }
            }
        }

        node = children.at(selectedIndex);
    }
    return node;
}

/** For all M+1 entries of a node N, compute the distance between the centers of their rectangles and the center of the bounding rectangle of N.
  * Sort the entries in decreasing order of their distances computed.
  * Remove the first p entries from N and adjust the bounding rectangle of N.
  * In the sort, starting with the maximum distance (= far reinsert) or minimum distance (= close reinsert), invoke Insert to reinsert the entries.
  */
void RStarTree::reInsert(RTreeNode &node, size_t &overflowLevel)
{
    overflowLevel |= (1<<node.level);

    node.prepareDisSort();
    if(node.level)
        sort(node.children->begin(), node.children->end());
    else
        sort(node.entries->begin(), node.entries->end());

    size_t p = RStarTree::maxChild * REINSERT_FACTOR;
    size_t reinsertFromIndex = node.size() - p;

    if(node.level) // If source is not a leaf node.
    {
        Node_P_V children;
        for(size_t ic = reinsertFromIndex; ic < node.children->size(); ic++)
            children.push_back(node.children->at(ic));
        node.resize(reinsertFromIndex);
        for(size_t ic = 0; ic < children.size(); ic++)
            insert(children.at(ic), NULL, node.level, overflowLevel);
    }
    else // If source is a leaf node.
    {
        Entry_P_V entries;
        for(size_t ie = reinsertFromIndex; ie < node.entries->size(); ie++)
            entries.push_back(node.entries->at(ie));
        node.resize(reinsertFromIndex);
        for(size_t ie = 0; ie < entries.size(); ie++)
            insert(NULL, entries.at(ie), node.level, overflowLevel);
    }
}

/** Invoke ChooseSplitAxis to determine the axis, perpendicular to which the split is performed.
  *
  * Invoke ChooseSplitIndex to determin the best distribution into two groups along that axis.
  *
  * Distribute the entries into two groups.
  */
void RStarTree::split(RTreeNode &node)
{
    size_t axis = chooseSplitAxis(node);
    size_t splitIndex = chooseSplitIndex(node, axis);

    Node_P newNode = new RTreeNode(node.level);
    newNode->take(node, splitIndex);
    if(!node.parent)
    {
        root = new RTreeNode(node.level + 1);
        root->insert(&node);
    }
    node.parent->insert(newNode);
}

/** For each axis:
  * Sort the entries by the lower then by the upper value of their rectangles;
  * And determin all distributions ( [m, M-m] elements in each partition );
  * Compute S, the sum of all margin-values of the different distributions.
  *
  * Choose the axis with the minimum S as split axis.
  *
  * Return:
  * 0 : choose x
  * 1 : choose y
  */
size_t RStarTree::chooseSplitAxis(RTreeNode &node)
{
    double minS = INF_P;
    size_t axis = 0;
    for(size_t dim = 0; dim < DIM; dim++)
    {
        if(node.level)
            sort(node.children->begin(), node.children->end(), AxisSort(dim));
        else
            sort(node.entries->begin(), node.entries->end(), AxisSort(dim));
        double S = computeS(node);
        if(S - EPS < minS)
        {
            minS = S;
            axis = dim;
        }
    }
    return axis;
}

double RStarTree::computeS(RTreeNode &node)
{
    double S = 0;
    size_t size = node.size(), last = size - 1;

    vector<Mbr> partition1(size);
    vector<Mbr> partition2(size);
    //Mbr partition1[size], partition2[size];
    if(node.level)
    {
        Node_P_V &children = *node.children;
        partition1[0] = Mbr::getMbr(children.at(0)->mbrn, children.at(0)->mbrn);
        for(size_t ic = 1; ic < size - minChild; ic++)
            partition1[ic] = Mbr::getMbr(partition1[ic - 1], children.at(ic)->mbrn);
        partition2[last] = Mbr::getMbr(children.at(last)->mbrn, children.at(last)->mbrn);
        for(size_t ic = last - 1; ic >= minChild; ic--)
            partition2[ic] = Mbr::getMbr(partition2[ic + 1], children.at(ic)->mbrn);
    }
    else
    {
        Entry_P_V &entries = *node.entries;
        partition1[0] = Mbr::getMbr(entries.at(0)->mbre, entries.at(0)->mbre);
        for(size_t ie = 1; ie < size - minChild; ie++)
            partition1[ie] = Mbr::getMbr(partition1[ie - 1], entries.at(ie)->mbre);
        partition2[last] = Mbr::getMbr(entries.at(last)->mbre, entries.at(last)->mbre);
        for(size_t ie = last - 1; ie >= minChild; ie--)
            partition2[ie] = Mbr::getMbr(partition2[ie + 1], entries.at(ie)->mbre);
    }

    // is : first element of a valid second partition
    for(size_t is = minChild; is <= size - minChild; is++)
        S += partition1[is - 1].getMargin() + partition2[is].getMargin();

    return S;
}

/** Along the chosen axis, choose the distribution with the minimum overlap-value.
  *
  * Resolve ties by choosing the distribution with minimum area-value.
  *
  * Return:
  * Index of the first element in the second partition
  */
size_t RStarTree::chooseSplitIndex(RTreeNode &node, size_t axis)
{
    size_t size = node.size(), last = size - 1;
    vector<Mbr> partition1(size);
    vector<Mbr> partition2(size);
    //Mbr partition1[size], partition2[size];
    if(node.level)
    {
        sort(node.children->begin(), node.children->end(), AxisSort(axis));
        Node_P_V &children = *node.children;
        partition1[0] = Mbr::getMbr(children.at(0)->mbrn, children.at(0)->mbrn);
        for(size_t ic = 1; ic < size - minChild; ic++)
            partition1[ic] = Mbr::getMbr(partition1[ic - 1], children.at(ic)->mbrn);
        partition2[last] = Mbr::getMbr(children.at(last)->mbrn, children.at(last)->mbrn);
        for(size_t ic = last - 1; ic >= minChild; ic--)
            partition2[ic] = Mbr::getMbr(partition2[ic + 1], children.at(ic)->mbrn);
    }
    else
    {
        sort(node.entries->begin(), node.entries->end(), AxisSort(axis));
        Entry_P_V &entries = *node.entries;
        partition1[0] = Mbr::getMbr(entries.at(0)->mbre, entries.at(0)->mbre);
        for(size_t ie = 1; ie < size - minChild; ie++)
            partition1[ie] = Mbr::getMbr(partition1[ie - 1], entries.at(ie)->mbre);
        partition2[last] = Mbr::getMbr(entries.at(last)->mbre, entries.at(last)->mbre);
        for(size_t ie = last - 1; ie >= minChild; ie--)
            partition2[ie] = Mbr::getMbr(partition2[ie + 1], entries.at(ie)->mbre);
    }

    double minOverlap = INF_P, minArea = INF_P;
    size_t splitIndex = minChild;
    // is : first element of a valid second partition
    for(size_t is = minChild; is <= size - minChild; is++)
    {
        double overlap = Mbr::getOverlap(partition1[is - 1], partition2[is]);
        if(overlap - EPS < minOverlap) // smaller or roughly equal overlap
        {
            double area = partition1[is - 1].getArea() + partition2[is].getArea();
            if(overlap + EPS > minOverlap) // roughly equal overlap
            {
                if(area < minArea)
                {
                    minArea = area;
                    splitIndex = is;
                }
            }
            else // smaller overlap
            {
                minOverlap = overlap;
                minArea = area;
                splitIndex = is;
            }
        }
    }

    return splitIndex;
}
}
