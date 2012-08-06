#ifndef _ID_H
#define _ID_H

#include <Array.H>


class ID
    :
    public Array<int>
{
    public:
        ID();
        
        ID(int length);
        
        //ID(const ID&);
        
        int level();
        
        ID parent();
        
        std::string toString() const;
};

std::ostream& operator<< (std::ostream& os, const ID& id);


#endif
