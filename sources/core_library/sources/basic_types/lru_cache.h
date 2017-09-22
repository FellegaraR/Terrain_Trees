/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Terrain Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Terrain Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LRU_CACHE_H
#define LRU_CACHE_H

#include <list>
#include <map>

using namespace std;

template<class Key, class Type> class LRU_Cache
{
public:
    typedef typename map<Key,Type>::iterator mapIt;

    LRU_Cache() {}
    LRU_Cache(unsigned cache_size) { this->max_size = cache_size; }
    inline void setMaxSize(utype cache_size) { this->max_size = cache_size; }
    inline utype getMaxSize() {  return this->max_size; }
    inline utype get_actual_size() { return lru_map.size(); }
    //do not check if the item is already present in the cache
    typename LRU_Cache<Key,Type>::mapIt insert(Key chiave, Type item);
    //
    inline typename LRU_Cache<Key,Type>::mapIt begin() { return this->lru_map.begin(); }
    inline typename LRU_Cache<Key,Type>::mapIt end() { return this->lru_map.end(); }
    typename LRU_Cache<Key,Type>::mapIt find(Key chiave);
    inline bool isFull() { return (lru_list.size() >= max_size); }
    inline void reset() {  this->lru_list.clear(); this->lru_map.clear(); }

private:
    utype max_size;
    list<Key> lru_list;
    map<Key,Type> lru_map;
};

template<class Key, class Type> typename LRU_Cache<Key,Type>::mapIt LRU_Cache<Key, Type>::find(Key chiave)
{
    typename map<Key,Type>::iterator it;
    it = lru_map.find(chiave);
    if(it != lru_map.end())//the key is into the cache
    {
        //update the list accordingly
        lru_list.remove(chiave);
        lru_list.push_front(chiave);
        return it;
    }
    else //non esiste
        return lru_map.end();
}

template<class Key, class Type> typename LRU_Cache<Key,Type>::mapIt LRU_Cache<Key,Type>::insert(Key chiave, Type item)
{
    //try to insert the element in cache only if it does not exist
    //if the cache is full remove the least recently used element
    pair<mapIt,bool> ret = lru_map.insert(make_pair(chiave,item));
    if(ret.second) // inserted
    {
        if(isFull())
        {
            Key last_key = lru_list.back();
            lru_list.pop_back();

            lru_map.erase(last_key);
        }

        lru_list.push_front(chiave);
    }
    else //the element was already into the cache, simply update the entry into the list
    {
        lru_list.remove(chiave);
        lru_list.push_front(chiave);
    }

    return ret.first;
}

#endif // LRU_CACHE_H
