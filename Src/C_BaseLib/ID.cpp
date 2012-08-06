#include <ID.H>

ID::ID()
    :
    Array<int>()
{
    
}

ID::ID(int length)
    :
    Array<int>(length, 0)
{

}

int 
ID::level()
{
    return this->size()-1;
}

ID
ID::parent()
{
    BL_ASSERT(size() > 1);
    ID parent = *this;
    parent.resize(size() - 1);
    return parent;
}

std::string
ID::toString() const
{
    std::stringstream ss;
    for (int i = 0; i < size()-1; i++)
    {
        ss << this->operator[](i) << "_";
    }
    if (size() > 0)
        ss << this->operator[](size()-1);
    return ss.str();
}

std::ostream& 
operator<< (std::ostream& os, const ID& id) 
{
    for (int i = 0; i < id.size() - 1; i++)
    {
        os << id[i]<< "_";
    }
    if (id.size() > 0)
        os << id[id.size()-1];
    return os;
}
