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

#ifndef LRU_CACHE_H
#define LRU_CACHE_H

#include <list>
#include <unordered_map>

using namespace std;

/**
 * @brief A class defining a Least-Recent-Used cache.
 * The class accepts two templates, one for the key of the cache, and one for the items encoded by the cache.
 */
template<class Key, class Type> class LRU_Cache
{
public:
    /**
     * @brief A public alias for using the iterator associated to the map forming the cache
     *
     */
    typedef typename unordered_map<Key,Type>::iterator mapIt;

    /**
     * @brief A constructor method
     *
     */
    LRU_Cache() {}
    /**
     * @brief A constructor method
     *
     * @param cache_size, an integer referring to the maximum cache size
     */
    LRU_Cache(int cache_size) { this->max_size = cache_size; }
    /**
     * @brief A public method that set the cache maximum size
     *
     * @param cache_size, an integer referring to the maximum cache size
     */
    inline void set_max_size(int cache_size) { this->max_size = cache_size; }
    /**
     * @brief A public method that returns the cache maximum size
     *
     * @return an integer referring to the cache maximum size
     */
    inline int get_max_size() {  return this->max_size; }
    /**
     * @brief A public method that returns the current size of the cache
     *
     * @return an integer containing the current size of the cache
     */
    inline int get_actual_size() { return lru_map.size(); }
    /**
     * @brief A public method that inserts an element in the cache, and returns the iterator to the element in the cache.
     * If the element is already into the cache, then simply updates its position in the LRU list (NOTA: it does NOT update the content of the element!).
     * If the cache is full the procedure removes the least-recent-used item in cache.
     * @param key, the key of the item to insert
     * @param item, the item to insert in cache
     * @return LRU_Cache<Key, Type>::mapIt, the iterator of the inserted item
     */
    typename LRU_Cache<Key,Type>::mapIt insert(Key key, Type item);
    //
    /**
     * @brief A public method that returns the iterator of the first item in cache
     *
     * @return LRU_Cache<Key, Type>::mapIt, the iterator to the first item
     */
    inline typename LRU_Cache<Key,Type>::mapIt begin() { return this->lru_map.begin(); }
    /**
     * @brief A public method that returns the iterator of the last item in cache
     *
     * @return LRU_Cache<Key, Type>, the iterator to the last item
     */
    inline typename LRU_Cache<Key,Type>::mapIt end() { return this->lru_map.end(); }
    /**
     * @brief A public method that search for an item in cache.
     * If the item is found then the procedure returns the iterator to the item,
     * otherwise it returns the iterator to the end of the cache
     *
     * @param key, the key of the item to search
     * @return LRU_Cache<Key, Type>, the iterator to the item with key chiave
     */
    typename LRU_Cache<Key,Type>::mapIt find(Key key);
    /**
     * @brief A public method that checks if the cache is full
     *
     * @return a boolean, true if the cache is full, false otherwise
     */
    inline bool is_full() { return (lru_list.size() >= max_size); }
    /**
     * @brief A public method that reset the cache
     *
     */
    inline void reset() {  this->lru_list.clear(); this->lru_map.clear(); }

private:    
    ///A private variable representing the maximum cache size
    unsigned int max_size;
    ///A private variable representing the access order of the element in cache
    list<Key> lru_list;
    ///A private variable representing the LRU-cache
    unordered_map<Key,Type> lru_map;
};

template<class Key, class Type> typename LRU_Cache<Key,Type>::mapIt LRU_Cache<Key, Type>::find(Key key)
{
    typename LRU_Cache<Key,Type>::mapIt it;
    it = lru_map.find(key);
    if(it != lru_map.end())
    {
        lru_list.remove(key);
        lru_list.push_front(key);
        return it;
    }
    else
        return lru_map.end();
}

template<class Key, class Type> typename LRU_Cache<Key,Type>::mapIt LRU_Cache<Key,Type>::insert(Key key, Type item)
{    
    pair<mapIt,bool> ret = lru_map.insert(make_pair(key,item));
    if(ret.second)
    {
        if(is_full())
        {
            Key last_key = lru_list.back();
            lru_list.pop_back();

            lru_map.erase(last_key);
        }

        lru_list.push_front(key);
    }
    else
    {
        lru_list.remove(key);
        lru_list.push_front(key);
    }

    return ret.first;
}

#endif // LRU_CACHE_H
