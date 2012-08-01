#ifndef BL_EXECUTIONTREE_H
#define BL_EXECUTIONTREE_H

#include <Tree.H>

class ExecutionTree
    :
    protected Tree<int>
{
    public:
    
    ExecutionTree(Array<int> root_id);
}

ExecutionTree::ExecutionTree(Array<int> root_id)
{
    root = new Node();
    root->m_id = root_id;
    int level = root_id.size() - 1;
    root->m_level = level;
    // I don't think we need to set index or parent here.
}

#endif /*BL_EXECUTIONTREE_H*/
