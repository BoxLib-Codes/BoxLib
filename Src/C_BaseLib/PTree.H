
#ifndef BL_TREE_H
#define BL_TREE_H

#include <cstddef>

#include <BLassert.H>
#include <PArray.H>
#include <PList.H>
#include <BoxLib.H>

enum PTreePolicy
{
    PTreeNoManage = 0,
    PTreeManage   = 1
};

template <class T> class PTreeIterator; // Forward Declaration

template <class T> class PTree;

/*
  An PTree of Pointers to Objects of Type T.

  When a PPTree<T> object is destructed, should only the pointers be
  destructed or should the objects to which they point also be deleted?
  The PPTree<T> class allows for both scenarios.  By default, it will
  only remove the pointers.  If the array is constructed with the
  PPTreeManage argument it will also delete all the objects pointed to
  in the array.
  
  Nodes are specified by ID's. An ID is an array of ints specifying the 
  position of the Node relative to its parents. The first int must always
  be 0 to represent the root node. For example, [0,1] represents the node
  at position 1 among the root node's children.
   
  We intentionally provide no methods for the deletion of single nodes 
  from the tree; either all children of a specified node must be 
  cleared, or none.

  This is a concrete class, not a polymorphic one.

  PTree<T>s do NOT provide a copy constructor or an assignment operator.
*/

template <class T>
class PTree
{
    friend class PTreeIterator<T>;
public:
    //
    // Creates a new empty tree with a root node and policy.
    //
    PTree(PTreePolicy is_managed = PTreeNoManage);
    //
    // Deconstructs the tree, deleting data if managed.
    //
    ~PTree();
    
    //
    // Returns the root data.
    //
    T* getRoot();
    
    //
    // Returns the data associated with the specified node
    //
    T* getData(const Array<int> id);
    
    //
    // 0's the pointer and returns the data associated with the specified node.
    // The user is now responsible for deleting the pointer.
    //
    T* removeData(const Array<int> id);
    
    //
    // Clears the data associated with the specified node. If managed,
    // the data pointed to will be deleted.
    //
    void clearData(const Array<int> id);
    
    //
    // Seths the data associated with the specified node
    //
    void setData(const Array<int> id, T* data);
    
    //
    // Adds a child to the specified node. Sets the child data if provided.
    //
    Array<int> addChildToNode(const Array<int> id);
    Array<int> addChildToNode(const Array<int> id, T* data);

    //
    // Gets the data from the children of the node -- only direct children, not descendents.
    //
    void getChildrenOfNode(const Array<int> id, PList<T>& children);
    
    //
    // Extracts the data from children of the node and removes the children (and all descendents).
    //
    void extractChildrenOfNode(const Array<int> id, PList<T>& extracted);
    
    //
    // Returns a PTreeIterator to root or the specified node. If restrict_level
    // is set, the iterator will skip over all but nodes at the specified level.
    //
    PTreeIterator<T> getIteratorAtRoot(int restrict_level = -1);
    PTreeIterator<T> getIteratorAtRoot(int restrict_level = -1) const;
    PTreeIterator<T> getIteratorAtNode(const Array<int> id, int restrict_level = -1);

    //
    // A hack I'm currently using to compute the finest level. Shouldn't be used in general.
    //
    int DeprecatedGetFinestLevel();
protected:
    struct Node
    {
        //
        // Creates a new node. Sets the parent, index and data if specified.
        //
        Node(PTreePolicy is_managed);
        Node(PTreePolicy is_managed, Node* parent, int index, int level);
        Node(PTreePolicy is_managed, Node* parent, int index, int level, T* data);
        
        //
        // Deconstructs this node and its children. If the nodes are managed,
        // the data will be deleted.
        //
        ~Node();
        
        //
        // Adds a new node as a child of this node. Fills in data if specified.
        // Managed policy will be the same as this node.
        //
        Array<int> addChild();
        Array<int> addChild(T* data);
        
        //
        // Returns the id for this node.
        //
        Array<int> getID();
        
        void extractChildren(PList<T>& extracted);
        
        T*              m_data;     // The data stored by this node.
        int             m_index;    // The index of this node in its parent's children.
        PArray<Node>    m_children; // The child nodes.
        Node*           m_parent;   // Pointer to the parent node.
        PTreePolicy     managed;    // The data memory management policy for the node.
        int             m_level;    // The level of this node w/in the tree (root = 0)
    };
    
    Node* root;
    Node* getNode(const Array<int> id);
    PTreePolicy managed;
};

/*
  This class provides an iterator-like interface for traversing a 
  subtree of some PTree. The traversal happens in postfix ordering--that
  is, a parent node is returned after all of its children have been 
  touched. Thus when dereferencing the iterator, the base node is the 
  very last node to be returned.
   
  If the iterator is created with restrict_level != -1, the iteration
  will automatically skip nodes not at the proper level.
   
  Standard usage is looks something like:
    for (PTreeIterator<T> it = myTree.getIteratorAtNode(id); !it.is_finished(); it++)
    {
        ...
    }

*/
template <class T>
class PTreeIterator
{
    typedef typename PTree<T>::Node iNode;
public:
    PTreeIterator(iNode* base, int r_level);
    
    void operator ++();
    T* operator *();
    
    bool is_finished();
    
protected:
    void descend();
    
    iNode* current_node;
    iNode* base_node;
    int current_index;
    bool finished;
    int restrict_level;
};

template <class T>
PTree<T>::PTree(PTreePolicy is_managed)
{
    managed=is_managed;
    root = new Node(managed);
}

template <class T>
PTree<T>::~PTree()
{
    delete root;
}

template <class T>
typename PTree<T>::Node*
PTree<T>::getNode(const Array<int> id)
{
    int N = id.size();
    BL_ASSERT(N > 0);
    BL_ASSERT(id[0] == 0);
    Node* node = root;
    for (int i = 1; i < N; i++)
    {
        BL_ASSERT(node->m_children.size() > id[i]);
        node = &(node->m_children[id[i]]);
    }
    return node;
}

template <class T>
T*
PTree<T>::getRoot()
{
    BL_ASSERT(root->m_data != 0);
    return root->m_data;
}

template <class T>
T*
PTree<T>::getData(const Array<int> id)
{
    T* data = getNode(id)->m_data;
    BL_ASSERT(data != 0);
    return data;
}

template <class T>
T*
PTree<T>::removeData(const Array<int> id)
{
    T* data = getNode(id)->m_data;
    BL_ASSERT(data != 0);
    getNode(id)->m_data = 0;
    return data;
}

template <class T>
void
PTree<T>::clearData(const Array<int> id)
{
    if (managed)
        delete (T*)getNode(id)->m_data;
    getNode(id)->m_data = 0;
}


template <class T>
void
PTree<T>::setData(const Array<int> id, T* data)
{
    BL_ASSERT(getNode(id)->m_data == 0);
    getNode(id)->m_data = data;
}

template <class T>
Array<int>
PTree<T>::addChildToNode(const Array<int> id)
{
    return getNode(id)->addChild();
}

template <class T>
Array<int>
PTree<T>::addChildToNode(const Array<int> id, T* data)
{
    return getNode(id)->addChild(data);
}

template <class T>
void
PTree<T>::getChildrenOfNode(const Array<int> id, PList<T>& children)
{
    Node* node = getNode(id);
    for (int i = 0; i < node->m_children.size(); i++)
    {
        children.push_back(node->m_children[i].m_data);
    }
}

template <class T>
void
PTree<T>::extractChildrenOfNode(const Array<int> id, PList<T>& extracted)
{
    getNode(id)->extractChildren(extracted);
}

template <class T>
PTreeIterator<T>
PTree<T>::getIteratorAtRoot(int restrict_level)
{
    PTreeIterator<T> new_it(root, restrict_level);
    return new_it;
}

template <class T>
PTreeIterator<T>
PTree<T>::getIteratorAtRoot(int restrict_level) const
{
    PTreeIterator<T> new_it(root, restrict_level);
    return new_it;
}

template <class T>
PTreeIterator<T>
PTree<T>::getIteratorAtNode(const Array<int> id, int restrict_level)
{
    return PTreeIterator<T>(getNode(id), restrict_level);
}

template <class T>
int
PTree<T>::DeprecatedGetFinestLevel()
{
    PTreeIterator<T> it = getIteratorAtRoot();
    int flev = 0;
    for ( ; !it.is_finished(); ++it)
    {
        if ((*it)->Level() > flev)
            flev = (*it)->Level() ;
    }
    return flev;
}

////////////////////////// Begin Node Methods ///////////////////////////
template <class T>
Array<int>
PTree<T>::Node::addChild()
{
    int num_children = m_children.size();
    m_children.resize(num_children+1);
    std::cout << "DEBUG: New index " << m_children.size()-1 << "\n";
    m_children.set(num_children, new Node(managed, this, num_children, m_level + 1));
    return m_children[num_children].getID();
}

template <class T>
Array<int>
PTree<T>::Node::addChild(T* data)
{
    int num_children = m_children.size();
    m_children.resize(num_children+1);
    std::cout << "DEBUG: New index " << m_children.size()-1 << "\n";
    m_children.set(num_children, new Node(managed, this, num_children, m_level + 1, data));
    return m_children[num_children].getID();
}

template <class T>
PTree<T>::Node::Node(PTreePolicy is_managed)
    :
    m_children(PArrayManage)
{
    managed = is_managed;
    m_parent = 0;
    m_index = 0;
    m_level = 0;
    m_data = 0;
}

template <class T>
PTree<T>::Node::Node(PTreePolicy is_managed, Node* parent, int index, int level)
    :
    m_children(PArrayManage)
{
    managed = is_managed;
    m_parent = parent;
    m_index = index;
    m_level = level;
    m_data = 0;
}

template <class T>
PTree<T>::Node::Node(PTreePolicy is_managed, Node* parent, int index, int level, T* data)
    :
    m_children(PArrayManage)
{
    managed = is_managed;
    m_parent = parent;
    m_index = index;
    m_level = level;
    m_data = data;
}

template <class T>
PTree<T>::Node::~Node()
{
    if (managed == PTreeManage)
    {
        delete (T*)m_data;
    }
    //m_children's destructors will be called when m_children is deleted.
}

template <class T>
Array<int>
PTree<T>::Node::getID()
{
    std::list<int> idl;
    Node* current_node = this;
    idl.push_front(current_node->m_index);
    while(current_node->m_parent != 0)
    {
        current_node = current_node->m_parent;
        idl.push_front(current_node->m_index);
    }
    // This could be cleaned up a bit.
    Array<int> id(idl.size());
    int i = 0;
    for (std::list<int>::iterator it = idl.begin(); it != idl.end(); ++it, ++i)
    {
        id[i] = *it;
    }
    return id;
}

template <class T>
void
PTree<T>::Node::extractChildren(PList<T>& extracted)
{
    for (int i = 0; i < m_children.size(); i++)
    {
        m_children[i].extractChildren(extracted);
        extracted.push_back(m_children[i].m_data);
        m_children[i].m_data = 0;
    }
    m_children.clear();
    m_children.resize(0);
}

////////////////////////Begin PTreeIterator Methods/////////////////////

template <class T>
PTreeIterator<T>::PTreeIterator(iNode* base, int r_level)
{
    base_node = base;
    restrict_level = r_level;
    BL_ASSERT(r_level == -1 || restrict_level >= base_node->m_level);
    current_index = 0;
    current_node = base_node;
    finished = false;
    descend();
    if ((restrict_level != -1 )&& (current_node->m_level != restrict_level))
        ++(*this);
}

template <class T>
T*
PTreeIterator<T>::operator *()
{
    return current_node->m_data;
}

template <class T>
inline
bool
PTreeIterator<T>::is_finished()
{
    return finished;
}

template <class T>
void
PTreeIterator<T>::operator ++()
{
    BL_ASSERT(!is_finished());
    if (current_node == base_node)
    {
        finished = true;
    }
    else
    {
        // Move up one level and record where you were.
        current_index = current_node->m_index + 1;
        current_node = current_node->m_parent;
        // If there are untraversed children at this node, go down
        // the next child. If not, stay here.
        descend();
        // If in restricted_level mode, advance again if not at the 
        // proper level
        if (restrict_level != -1 && restrict_level != current_node->m_level)
            ++(*this);
    }
}

template <class T>
void
PTreeIterator<T>::descend()
{
    if (current_node->m_children.size() > current_index)
    {
        current_node = &(current_node->m_children[current_index]);
        current_index = 0;
        descend();
    }
}


#endif /*BL_TREE_H*/
