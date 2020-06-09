/*
    This file is part of the Stellar library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Stellar library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Stellar library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Stellar library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONTAINER_UTILITIES_H
#define CONTAINER_UTILITIES_H

#include <algorithm>
#include <iostream>
#include <set>
#include <boost/dynamic_bitset.hpp>
using namespace std;

/**
 * @brief A procedure that counts the number of elements in a container of containers
 *
 * @param c, the container of containers
 * @return the summation of the elements
 */
template<class C> int get_num_elements_in_container_of_containers(C& c)
{
    int sum = 0;
    for(typename C::const_iterator it=c.begin(); it!=c.end(); ++it)
    {
        sum += it->size();
    }
    return sum;
}
/**
 * @brief A procedure that prints on the standard output the container content
 *
 * @param c, the container
 */
template<class C> void print_container_content(const C& c)
{
    for(typename C::const_iterator it=c.begin(); it!=c.end(); ++it)
        cout<<*it<<" ";
}
/**
 * @brief A procedure that prints on the standard output the container content followed by a user-defined caption
 *
 * @param caption, a string representing the caption
 * @param c, the container
 */
template<class C> void print_container_content(string caption, const C& c)
{
    cout<<caption;
    print_container_content(c);
    cout<<endl;
}
/**
 * @brief A procedure that prints on the standard output the container-of-containers content
 *
 * @param c, the container-of-containers
 */
template<class C> void print_container_of_containers_content(C& c)
{
    int i = 0;
    for(typename C::const_iterator it=c.begin(); it!=c.end(); ++it)
    {
        if(it->size() >0)
        {
            cout<<"  C["<<i<<"] ";
            print_container_content(*it);
            cout<<endl;
        }
        i++;
    }
}
/**
 * @brief A procedure that prints on the standard output the container-of-containers content followed by a user-defined caption
 *
 * @param caption, a string representing the caption
 * @param c, the container-of-containers
 */
template<class C> void print_container_of_containers_content(string caption, C& c)
{
    cout<<caption<<endl;
    print_container_of_containers_content(c);
}


/**
 * @brief A procedure that sorts a container
 *
 * @param c
 */
template<class C> void sort_container(C& c) { sort(c.begin(),c.end()); }
/**
 * @brief A procedure that sorts a container-of-containers
 *
 * @param c
 */
template<class C> void sort_container_of_containers(C& c)
{
    for(typename C::iterator it=c.begin(); it!=c.end(); ++it)
    {
        sort_edges_tuples(*it);
    }
}
/**
 * @brief A procedure that intersects two sets and places the resulting set in the first set
 *
 * @param s1
 * @param s2
 */
template<class C> void intersect_sets(set<C> &s1, set<C> &s2)
{
    set<C> res;
    set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(res,res.begin()));
    s1.clear();
    s1 = res;
}
/**
 * @brief A procedure that intersects two sets and places the resulting set in a third set
 *
 * @param s1
 * @param s2
 * @param res
 */
template<class C> void intersect_sets(set<C> &s1, set<C> &s2, set<C> &res)
{
    set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(res,res.begin()));
}
/**
 * @brief A procedure that intersects two containers and places the resulting set in the first container
 *        NOTA: the containers must be sorted!
 * @param v1
 * @param v2
 */
template<class C> void intersect_containers(C &v1, const C &v2)
{
    C res;
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(res));
    v1.clear();
    v1 = res;
}
/**
 * @brief A procedure that intersects two containers and ..... ......
 *        NOTA: the containers must be sorted!
 * @param v1
 * @param v2
 */
template<class C> void intersect_containers(const C &v1, const C &v2, C &res)
{
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(res));
}

/**
 * @brief A procedure that intersects two containers-of-containers and places the resulting set in the first container
 *
 * @param v1
 * @param v2
 * @param sort, true if the two containers-of-containers must be sorted
 */
///
template<class C> void intersect_container_of_containers(C &v1, C &v2, bool sort=false)
{
    if(sort)
    {
        sort_container_of_containers(v1);
        sort_container_of_containers(v2);
    }

    for(unsigned i=0; i<v1.size(); i++)
    {
        intersect_containers(v1[i],v2[i]);
    }
}
/**
 * @brief A procedure that intersects two containers-of-containers .... ..... .....
 *
 * @param v1
 * @param v2
 * @param sort, true if the two containers-of-containers must be sorted
 */
///
template<class C> void intersect_container_of_containers(C &v1, C &v2, C &res, bool sort=false)
{
    if(sort)
    {
        sort_container_of_containers(v1);
        sort_container_of_containers(v2);
    }

    for(unsigned i=0; i<v1.size(); i++)
    {
        intersect_containers(v1[i],v2[i],res[i]);
    }
}
/**
 * @brief A procedure that checks if exists an intersection of two containers
 *        NOTA: the containers must be sorted!
 * @param v1
 * @param v2
 * @return true if an intersection exist, false otherwise
 */
template<class C> bool exists_intersection_of_containers(C &v1, C &v2)
{
    C res;
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(res));
    return res.size() > 0;
}
/**
 @brief A procedure that checks if exists an intersection of two containers-of-containers
 *
 * @param v1
 * @param v2
 * @return true if an intersection exist, false otherwise
 */
template<class C> bool exists_intersection_of_container_of_containers(C &v1, C &v2)
{
    for(unsigned i=0; i<v1.size(); i++)
    {
        sort_edges_tuples(v1[i]);
        sort_edges_tuples(v2[i]);

        if(exists_intersection_of_containers(v1[i],v2[i]))
        {
            return true;
        }
    }
    return false;
}
/**
 * @brief A procedure that unifies two sets and places the resulting set in the first set
 *
 * @param s1
 * @param s2
 */
template<class C> void unify_sets(set<C> &s1, set<C> &s2)
{
    set<C> res;
    set_union(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(res,res.end()));
    s1.clear();
    s1 = res;
}
/**
 * @brief A procedure that unifies two containers and places the resulting set in the first container
 *        NOTA: the containers must be sorted!
 * @param c1
 * @param c2
 */
template<class C> void unify_containers(C &v1, C &v2)
{
    C res;
    set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(res));
    v1.clear();
    v1 = res;
}
/**
 * @brief A procedure that unifies two containers-of-containers and places the resulting set in the first container-of-containers
 *
 * @param c1
 * @param c2
 */
template<class C> void unify_container_of_containers(C &v1, C &v2)
{
    sort_container_of_containers(v1);
    sort_container_of_containers(v2);

    for(unsigned i=0; i<v1.size(); i++)
    {
        unify_containers(v1[i],v2[i]);
    }
}
/**
 * @brief A procedure that makes the difference of two containers and places the resulting set in the first container
 *        NOTA: the containers must be sorted!
 * @param c1
 * @param c2
 */
template<class C> void difference_of_containers(C &v1, C &v2)
{
    C res;
    set_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(res));
    v1.clear();
    v1 = res;
}
/**
 * @brief A procedure that makes the difference of two containers-of-containers and places the resulting set in the first container-of-containers
 *
 * @param c1
 * @param c2
 */
template<class C> void difference_of_container_of_containers(C &v1, C &v2)
{
    sort_container_of_containers(v1);
    sort_container_of_containers(v2);

    for(unsigned i=0; i<v1.size(); i++)
    {
        difference_of_containers(v1[i],v2[i]);
    }
}
/**
 * @brief A procedure that makes the difference of two sets and places the resulting set in the first set
 *
 * @param s1
 * @param s2
 */
template<class C> void difference_of_sets(set<C> &s1, set<C> &s2)
{
    set<C> res;
    set_difference(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(res,res.end()));
    s1.clear();
    s1 = res;
}
/**
 * @brief A procedure that makes the difference of two sets and places the resulting set in a third set
 *
 * @param s1
 * @param s2
 * @param res
 */
template<class C> void difference_of_sets(set<C> &v1, set<C> &v2, set<C> &res)
{
    set_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),std::inserter(res,res.end()));
}
/**
 * @brief A procedure that erase a value from a container
 *        NOTA: the procedure erase the first occurrence of the value.
 * @param c
 * @param value
 */
template<class T, class C> void erase_value_from_container(C& c, T value)
{
    for(typename C::iterator it=c.begin(); it!=c.end(); )
    {
        if(*it == value)
        {
            c.erase(it);
            return;
        }
        else
            ++it;
    }
}
/**
 * @brief A procedure that search for a value in a container
 *
 * @param c
 * @param value
 * @return true, if value is found, false otherwise
 */
template<typename T, class C> bool is_into_container(C& c, T value)
{
    return (find(c.begin(),c.end(),value) != c.end());
}


/**
  * @brief A procedure that efficiently check if a container 'a' contains container 'b'
  * This implements the 'classical' set-containment operation, that seems to be missing in C++ STL
  * NOTA: the containers must be sorted!
  *
  * @param a
  * @param b
  * @return true, if value a contains b, false otherwise
  */
template<class C> bool contains_container(C &a, C&b)
{
    boost::dynamic_bitset<> db(b.size());
    size_t a_id_pos = 0;
//    int stop;
//    print_container_content("a: ",a);
//    print_container_content("b: ",b);
    for (size_t i=0; i<b.size(); i++)
    {
        for (size_t j = a_id_pos; j<a.size(); j++)
        {
            if(b[i] == a[j])
            {
                db.set(i);
                a_id_pos = j + 1;
                break;
            }
            /// as the two containers are sorted..
            /// when I encounter the first entry in a greater than the key in b
            /// I can say that a does not contains b
            if(b[i] < a[j])
            {
//                cerr<<"[false] a does not contain b"<<endl;
//                cin>>stop;
                return false;
            }
        }
    }
//    cerr<<"a contains b: "<<(db.count() == db.size())<<endl;
//    cin>>stop;
    return db.count() == db.size();
}


#endif // CONTAINER_UTILITIES_H
