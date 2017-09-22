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

#include "string_management.h"
#include "basic_wrappers.h"

#include <iostream>
#include <limits>

namespace string_management
{

string get_file_name(string path)
{
    // Strip the path, isolate file name
    string filename = strip_path(path);

    // remove the extension if any
    string::size_type point_index = filename.find_last_of('.');
    if (point_index == string::npos)
        return filename;
    else
        return filename.substr(0, point_index);
}

string get_file_extension(string path)
{
    // Strip the path, isolate file name
    string filename = strip_path(path);
    // remove the extension if any
    string::size_type point_index = filename.find_last_of('.');
    return filename.substr(point_index+1,filename.size());
}

string get_path_without_file_extension(string path)
{
    // remove the extension if any
    utype point_index = path.find_last_of('.');
    if (point_index == string::npos)
        return path;
    else
        return path.substr(0, point_index);
}

string strip_path(string path)
{
    string filename;
    string::size_type name_index_ms = path.find_last_of('\\');
    string::size_type name_index_unix = path.find_last_of('/');
    if (name_index_ms != string::npos)
        filename = path.substr(name_index_ms + 1);
    else if (name_index_unix != string::npos)
        filename = path.substr(name_index_unix + 1);
    else
        filename = path;
    return filename;
}

string get_path(string path)
{
    string p;
    string::size_type name_index_ms = path.find_last_of('\\');
    string::size_type name_index_unix = path.find_last_of('/');
    if (name_index_ms != string::npos)
        p = path.substr(0, name_index_ms);
    else if (name_index_unix != string::npos)
        p = path.substr(0, name_index_unix);
    else
        p = "";
    return p;
}

void tokenize(const string& str, vector<string>& tokens, const string& delimiter = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}

void go_to_line(ifstream& input, unsigned num)
{
    for(unsigned i=0; i < num; ++i){
        input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
}

}
