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

#ifndef USAGE_H
#define USAGE_H

#include <sys/types.h>
#include <sys/resource.h>

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>

using namespace std;

/**
 * @brief The MemoryUsage class keeps track of the memory (peak) usage at run-time
 * The class is specified to work on both linux and apple systems
 */
class MemoryUsage
{

private:
    int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;

public:
    /**
     * @brief A public method that returns the current resident memory peak (RAM memory) in kilobytes
     * @param print a boolean, if true prints on the standard output a string containing the value of the memory peak
     * @return a float containing the memory peak
     */
    inline float get_Resident_Memory_in_KB(bool print=false)
    {
         ret=getrusage(who,&usage);
#ifdef __APPLE__
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0) << " KB" << endl;
         return usage.ru_maxrss/(1024.0);
#else
         if(print) cout << "Memory Usage: " << usage.ru_maxrss << " KB" << endl;
         return usage.ru_maxrss;
#endif
    }
    /**
     * @brief A public method that return the current memory peak in megabytes
     * @param print a boolean, if true prints on the standard output a string containing the value of the memory peak
     * @return a float containing the memory peak
     */
    inline float get_Resident_Memory_in_MB(bool print=false)
    {
         ret=getrusage(who,&usage);
#ifdef __APPLE__
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0*1024.0) << " MB" << endl;
         return usage.ru_maxrss/(1024.0*1024.0);
#else
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0) << " MB" << endl;
         return usage.ru_maxrss/(1024.0);
#endif
    }
    /**
     * @brief A public method that return the current memory peak in gigabytes
     * @param print a boolean, if true prints on the standard output a string containing the value of the memory peak
     * @return a float containing the memory peak
     */
    inline float get_Resident_Memory_in_GB(bool print=false)
    {
         ret=getrusage(who,&usage);
#ifdef __APPLE__
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0*1024.0*1024.0) << " GB" << endl;
         return usage.ru_maxrss/(1024.0*1024.0*1024.0);
#else
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0*1024.0) << " GB" << endl;
         return usage.ru_maxrss/(1024.0*1024.0);
#endif
    }

    /**
     * @brief A public method that returns the current virtual memory peak (RAM+SWAP memories) in kilobytes
     *
     * VmHWM has the same value of ru_maxrss
     * while VmSwap contains the value for the swap memory
     *
     * NOTA: the swap memory is not yet considered in APPLE workstations
     *
     * @return the memory peak in kilobytes
     */
    inline float get_Virtual_Memory_in_KB()
    {
#ifdef __APPLE__
        // to fix
        return this->get_Resident_Memory_in_KB(false);
#else
        //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        float result = 0;
        char line[128];
        int counter = 0;


        while (fgets(line, 128, file) != NULL)
        {
            if(strncmp(line, "VmHWM:", 6) == 0 || strncmp(line, "VmSwap:", 7) == 0)
            {
                result += get_vm_value(line);
                counter++;
            }
            if(counter == 2)
                break;
        }
        fclose(file);
        return result;
#endif
    }

    /**
     * @brief A public method that returns the current virtual memory peak (RAM+SWAP memories) in megabytes
     * NOTA: the swap memory is not yet considered in APPLE workstations
     *
     * @return the memory peak in megabytes
     */
    inline float get_Virtual_Memory_in_MB()
    {
#ifdef __APPLE__
        // to fix
        return this->get_Resident_Memory_in_MB(false);
#else
        return get_Virtual_Memory_in_KB()/(1024.0);
#endif
    }

    /**
     * @brief A public method that returns the current virtual memory peak (RAM+SWAP memories) in gigabytes
     * NOTA: the swap memory is not yet considered in APPLE workstations
     *
     * @return the memory peak in gigabytes
     */
    inline float get_Virtual_Memory_in_GB()
    {
#ifdef __APPLE__
        // to fix
        return this->get_Resident_Memory_in_GB(false);
#else
        return get_Virtual_Memory_in_KB()/(1024.0*1024.0);
#endif
    }

private:
    /**
     * @brief A private method returning the value allocated memory
     * NOTA: works only on linux workstations using the /proc/<pid>/... paradigm
     * @param line, a variable containing the current line
     * @return the memory peak in kilobytes
     */
    inline float get_vm_value(char* line)
    {
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = std::char_traits<char>::length(line);
        char *p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        return atof(p);
    }
};

#endif // USAGE_H
