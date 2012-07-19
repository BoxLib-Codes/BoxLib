
#ifndef BL_TREE_H
#define BL_TREE_H

#include <cstddef>

#include <BLassert.H>
#include <PArray.H>
#include <PList.H>
#include <BoxLib.H>


enum IteratorMode
{
    Prefix = 0,
    Postfix = 1
};
template <class T> class TreeIterator; // Forward Declaration

template <class T> class Tree;

/*
  An Tree of Pointers to Objects of Type T.

  When a PTree<T> object is destructed, should only the pointers be
  destructed or should the objects to which they point also be deleted?
  The PTree<T> class allows for both scenarios.  By default, it will
  only remove the pointers.  If the array is constructed with the
  PTreeManage argument it will also delete all the objects pointed to
  in the array.
  
  Nodes are specified by ID's. An ID is an array of ints specifying the 
  position of the Node relative to its parents. The first int must always
  be 0 to represent the root node. For example, [0,1] represents the node
  at position 1 among the root node's children.
   
  We intentionally provide no methods for the deletion of single nodes 
  from the tree; either all children of a specified node must be 
  cleared, or none.

  This is a concrete class, not a polymorphic one.

  Tree<T>s do NOT provide a copy constructor or an assignment operator.
*/

template <class T>
class Tree
{
    friend class TreeIterator<T>;
    friend class TreeConstIterator<T>;
public:
    //
    // Creates a new empty tree with a root node and policy.
    //
    Tree();
    //
    // Deconstructs the tree, deleting data if managed.
    //
    ~Tree();
    
    //
    // Returns the root data.
    //
    T getRoot();
    
    //
    // Returns the data associated with the specified node
    //
    T getData(const Array<int> id);
    
    //
    // Seths the data associated with the specified node
    //
    void setData(const Array<int> id, T data);
    
    //
    // Adds a child to the specified node. Sets the child data if provided.
    //
    Array<int> addChildToNode(const Array<int> id);
    Array<int> addChildToNode(const Array<int> id, T data);

    //
    // Gets the data from the children of the node -- only direct children, not descendents.
    //
    void getChildrenOfNode(const Array<int> id, PList<T>& children);
    
    //
    // Returns a TreeIterator to root or the specified node. If restrict_level
    // is set, the iterator will skip over all but nodes at the specified level.
    //
    TreeIterator<T> getIteratorAtRoot(int restrict_level = -1, IteratorMode mode = Postfix);
    TreeIterator<T> getIteratorAtNode(const Array<int> id, int restrict_level = -1, IteratorMode mode = Postfix);

    //
    // A hack I'm currently using to compute the finest level. Shouldn't be used in general.
    //
    int DeprecatedGetFinestLevel();
    
    //
    // Extracts a list representing the structure of the subtree rooted at the specified id.
    // The structure is a list giving the number of children in prefix_notation.
    //
    std::list<int> getStructure(const Array<int> id);
    
    //
    // Builds an empty subtree beneath the specified node based on the supplied structure list.
    //
    void buildFromStructure(const Array<int> id, std::list<int> structure);
    
    //
    // Returns the number of nodes in the tree
    //
    int numElements() const;
protected:
    struct Node
    {
        //
        // Creates a new node. Sets the parent, index and data if specified.
        //
        Node(TreePolicy is_managed);
        Node(TreePolicy is_managed, Node* parent, int index, int level);
        Node(TreePolicy is_managed, Node* parent, int index, int level, T* data);
        
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
        Array<int> addChild(T data);
        
        //
        // Returns the id for this node.
        //
        Array<int> getID();
        
        //
        // Recursively extracts the data from children of the node and removes the children.
        //
        void extractChildren(PList<T>& extracted);
        
        //
        // Recursively adds the structure for this node and its children to the structure list.
        //
        void getStructure(std::list<int>& structure);
        
        //
        // Recursively builds the node from the specified position in the structure list.
        //
        void buildFromStructure(std::list<int>& structure, std::list<int>::iterator& it);
        
        T               m_data;     // The data stored by this node.
        int             m_index;    // The index of this node in its parent's children.
        PArray<Node>    m_children; // The child nodes.
        Node*           m_parent;   // Pointer to the parent node.
        int             m_level;    // The level of this node w/in the tree (root = 0)
    };
    
    Node* root;
    Node* getNode(const Array<int> id);
    const Node* getNode(const Array<int> id) const;
    TreePolicy managed;
};

/*
  This class provides an iterator-like interface for traversing a 
  subtree of some Tree. The traversal happens in postfix ordering--that
  is, a parent node is returned after all of its children have been 
  touched. Thus when dereferencing the iterator, the base node is the 
  very last node to be returned.
   
  If the iterator is created with restrict_level != -1, the iteration
  will automatically skip nodes not at the proper level.
   
  Standard usage is looks something like:
    for (TreeIterator<T> it = myTree.getIteratorAtNode(id); !it.is_finished(); it++)
    {
        ...
    }

*/
template <class T>
class TreeIterator
{
    typedef typename Tree<T>::Node iNode;
public:
    TreeIterator(iNode* base, int r_level, IteratorMode mode);
    
    void operator ++();
    T operator *();
    
    bool is_finished() const ;
    
    int getLevel() const;
    
    Array<int> getID() const;
    
protected:
    void descend();
    
    iNode* current_node;
    iNode* base_node;
    int current_index;
    bool finished;
    int restrict_level;
    IteratorMode m_mode;
};

template <class T>
Tree<T>::Tree(TreePolicy is_managed)
{
    managed=is_managed;
    root = new Node(managed);
}

template <class T>
Tree<T>::~Tree()
{
    delete root;
}

template <class T>
typename Tree<T>::Node*
Tree<T>::getNode(const Array<int> id)
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
const typename Tree<T>::Node*
Tree<T>::getNode(const Array<int> id) const
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
    const Node* out_node = node;
    return out_node;
}

template <class T>
T
Tree<T>::getRoot() const
{
    return root->m_data;
}

template <class T>
T
Tree<T>::getData(const Array<int> id) const
{
    return getNode(id)->m_data;
}

template <class T>
void
Tree<T>::setData(const Array<int> id, T data)
{
    getNode(id)->m_data = data;
}

template <class T>
Array<int>
Tree<T>::addChildToNode(const Array<int> id)
{
    return getNode(id)->addChild();
}

template <class T>
Array<int>
Tree<T>::addChildToNode(const Array<int> id, T data)
{
    return getNode(id)->addChild(data);
}

template <class T>
void
Tree<T>::getChildrenOfNode(const Array<int> id, PList<T>& children)
{
    Node* node = getNode(id);
    for (int i = 0; i < node->m_children.size(); i++)
    {
        children.push_back(node->m_children[i].m_data);
    }
}

template <class T>
TreeIterator<T>
Tree<T>::getIteratorAtRoot(int restrict_level, IteratorMode mode)
{
    TreeIterator<T> new_it(root, restrict_level,mode);
    return new_it;
}

template <class T>
TreeIterator<T>
Tree<T>::getIteratorAtNode(const Array<int> id, int restrict_level, IteratorMode mode)
{
    return TreeIterator<T>(getNode(id), restrict_level,mode);
}

template <class T>
int
Tree<T>::DeprecatedGetFinestLevel()
{
    TreeIterator<T> it = getIteratorAtRoot();
    int flev = 0;
    for ( ; !it.is_finished(); ++it)
    {
        if ((*it)->Level() > flev)
            flev = (*it)->Level() ;
    }
    return flev;
}

template <class T>
std::list<int> 
Tree<T>::getStructure(const Array<int> id)
{
    // Create the list
    std::list<int> structure;
    // Recursively fill it
    getNode(id)->getStructure(structure);
    // Return it.
    return structure;
}

template <class T>
void 
Tree<T>::buildFromStructure(const Array<int> id, std::list<int> structure)
{
    Node* current_node = getNode(id);
    BL_ASSERT(current_node->m_children.size() == 0);
    std::list<int>::iterator it = structure.begin();
    current_node ->buildFromStructure(structure, it);
    BL_ASSERT(it == structure.end());
}

template <class T>
int
Tree<T>::numElements() const
{
    TreeConstIterator<T> ptree_it = getConstIteratorAtRoot();
    int cnt = 0;
    for( ; !ptree_it.is_finished(); ++ptree_it)
        cnt++;
    return cnt;
}

////////////////////////// Begin Node Methods ///////////////////////////
template <class T>
Array<int>
Tree<T>::Node::addChild()
{
    int num_children = m_children.size();
    m_children.resize(num_children+1);
    std::cout << "DEBUG: New index " << m_children.size()-1 << "\n";
    m_children.set(num_children, new Node(this, num_children, m_level + 1));
    return m_children[num_children].getID();
}

template <class T>
Array<int>
Tree<T>::Node::addChild(T data)
{
    int num_children = m_children.size();
    m_children.resize(num_children+1);
    std::cout << "DEBUG: New index " << m_children.size()-1 << "\n";
    m_children.set(num_children, new Node(this, num_children, m_level + 1, data));
    return m_children[num_children].getID();
}

template <class T>
Tree<T>::Node::Node()
    :
    m_children(PArrayManage)
{
    m_parent = 0;
    m_index = 0;
    m_level = 0;
    m_data = 0;
}

template <class T>
Tree<T>::Node::Node( Node* parent, int index, int level)
    :
    m_children(PArrayManage)
{
    m_parent = parent;
    m_index = index;
    m_level = level;
    m_data = 0;
}

template <class T>
Tree<T>::Node::Node( Node* parent, int index, int level, T data)
    :
    m_children(PArrayManage)
{
    m_parent = parent;
    m_index = index;
    m_level = level;
    m_data = data;
}

template <class T>
Tree<T>::Node::~Node()
{
}

template <class T>
Array<int>
Tree<T>::Node::getID()
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
Tree<T>::Node::extractChildren(PList<T>& extracted)
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

template <class T>
void
Tree<T>::Node::getStructure(std::list<int>& structure)
{
    int N = m_children.size();
    structure.push_back(N);
    for (int i = 0; i < N; i++)
    {
        m_children[i].getStructure(structure);
    }
}

template <class T>
void
Tree<T>::Node::buildFromStructure(std::list<int>& structure, std::list<int>::iterator& it)
{
    BL_ASSERT(it != structure.end());
    int N = *it;
    ++it;
    for (int i = 0; i < N; i++)
    {
        addChild();
        m_children[i].buildFromStructure(structure, it);
    }
}

////////////////////////Begin TreeIterator Methods/////////////////////

template <class T>
TreeIterator<T>::TreeIterator(iNode* base, int r_level, IteratorMode mode)
{
    base_node = base;
    m_mode = mode;
    restrict_level = r_level;
    BL_ASSERT(r_level == -1 || restrict_level >= base_node->m_level);
    current_index = 0;
    current_node = base_node;
    finished = false;
    if (m_mode == Postfix)
        descend();
    if ((restrict_level != -1 )&& (current_node->m_level != restrict_level))
        ++(*this);
}

template <class T>
T
TreeIterator<T>::operator *()
{
    return current_node->m_data;
}

template <class T>
inline
bool
TreeIterator<T>::is_finished() const
{
    return finished;
}

template <class T>
void
TreeIterator<T>::operator ++()
{
    BL_ASSERT(!is_finished());
    if (current_node == base_node && (m_mode == Postfix || current_index == current_node->m_children.size()))
    {
        finished = true;
    }
    else
    {
        if (m_mode == Postfix)
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
        else
        {
            ascend();
            current_index ++;
            current_node = &current_node.m_children[current_index];
        }
    }
}

template <class T>
void
TreeIterator<T>::descend()
{
    if (current_node->m_children.size() > current_index)
    {
        current_node = &(current_node->m_children[current_index]);
        current_index = 0;
        descend();
    }
}

template <class T>
void
TreeIterator<T>::ascend()
{
    if (current_node->m_children.size() == current_index && current_node != base_node)
    {
        current_index = current_node->m_index;
        current_node = current_node->m_parent;
        ascend();
    }
}

template <class T>
int
TreeIterator<T>::getLevel() const
{
    return current_node->m_level;
}

template <class T>
Array<int>
TreeIterator<T>::getID() const
{
    return current_node->getID();
}
#endif /*BL_TREE_H*/
