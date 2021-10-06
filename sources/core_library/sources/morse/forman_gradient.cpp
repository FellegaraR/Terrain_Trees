/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Federico Iuricich (federico.iuricich@gmail.com)

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

#include "forman_gradient.h"
#include "utilities/sorting.h"
//#include <assert.h>

Forman_Gradient::Forman_Gradient(itype num_t/*, ivect &original_v_indexes*/)
{
    forman_gradient = vector<ushort>(num_t,OUTSIDECASES);

    ushort tCaseCountOuterAllowed=0;
    CompressedToExpandedLUT compressed_to_expanded_local = CompressedToExpandedLUT(CASES);
    ExpandedToCompressedLUT expanded_to_compressed_local = ExpandedToCompressedLUT();
    compressed_to_expanded = CompressedToExpandedLUT(CASES);
    expanded_to_compressed = ExpandedToCompressedLUT();

    for(unsigned i = 0; i < 511; i++)
    {

        Triangle_Arrows ga = static_cast<Triangle_Arrows>( static_cast<unsigned>( i ));
        if( is_valid_case(ga) )
        {
            compressed_to_expanded_local[tCaseCountOuterAllowed]       = ga;
            expanded_to_compressed_local[ga]                             = tCaseCountOuterAllowed;
            tCaseCountOuterAllowed++;
        }
    }

    compressed_to_expanded = compressed_to_expanded_local;
    expanded_to_compressed = expanded_to_compressed_local;
}

bool Forman_Gradient::is_valid_case(Triangle_Arrows arrow/*, bool allowOuterFaces = false*/)
{
    // This disables the outer four bits (corresponding to the FT relations for the opposite T to each F) from being valid

    // test vertices
    if( BitTwiddle::popcount( arrow & CONTAINS_V0 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_V1 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_V2 ) > 1 )		return false;

    // test edges
    if( BitTwiddle::popcount( arrow & CONTAINS_E01 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_E02 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_E12 ) > 1 )		return false;

    // test faces
    if( BitTwiddle::popcount( arrow & CONTAINS_F012 ) > 1 )		return false;

    // valid!
    return true;
}

void Forman_Gradient::set_ET(itype t, const ivect &edge, Mesh &mesh)
{
    short e = mesh.get_triangle(t).edge_index(edge);
    ///
    TriGradient ga = convert_compressed_to_expand(t);
    ga.setEF(e);
    forman_gradient[t-1] = convert_expand_to_compressed(ga.getArrow());
}

void Forman_Gradient::set_VE(itype v, itype v2, leaf_ET &et, Mesh &mesh)
{

    ivect e = {v,v2};
  
    sort(e.begin(),e.end());
    leaf_ET::iterator it = et.find(e);

    TriGradient ga = convert_compressed_to_expand(it->second.first);

    Triangle &tri1 = mesh.get_triangle(it->second.first);
    ga.setVE(tri1.vertex_index(v), tri1.vertex_index(v2));

    forman_gradient[it->second.first-1] = convert_expand_to_compressed(ga.getArrow());

    if(it->second.second != -1)
    {
        TriGradient gb = convert_compressed_to_expand(it->second.second);
        Triangle &tri2 = mesh.get_triangle(it->second.second);
        gb.setVE(tri2.vertex_index(v), tri2.vertex_index(v2));
        forman_gradient[it->second.second-1] = convert_expand_to_compressed(gb.getArrow());
    }

}

void Forman_Gradient::set_VE(itype v, itype v2, pair<itype,itype> &et, Mesh &mesh)
{
    TriGradient ga = convert_compressed_to_expand(et.first);
    ga.setVE(mesh.get_triangle(et.first).vertex_index(v), mesh.get_triangle(et.first).vertex_index(v2));
    forman_gradient[et.first-1] = convert_expand_to_compressed(ga.getArrow());

    if(et.second != -1)
    {
        TriGradient gb = convert_compressed_to_expand(et.second);
        gb.setVE(mesh.get_triangle(et.second).vertex_index(v), mesh.get_triangle(et.second).vertex_index(v2));
        forman_gradient[et.second-1] = convert_expand_to_compressed(gb.getArrow());
    }
}

void Forman_Gradient::free_VE(itype v1, itype v2, pair<itype, itype> &et, Mesh &mesh)
{
    TriGradient ga = convert_compressed_to_expand(et.first);
    Triangle &tri = mesh.get_triangle(et.first);
    ga.clearVE(tri.vertex_index(v1), tri.vertex_index(v2));
    forman_gradient[et.first-1] = convert_expand_to_compressed(ga.getArrow());

    if(et.second != -1)
    {
        TriGradient gb = convert_compressed_to_expand(et.second);
        Triangle &tri = mesh.get_triangle(et.second);
        gb.clearVE(tri.vertex_index(v1), tri.vertex_index(v2));
        forman_gradient[et.second-1] = convert_expand_to_compressed(gb.getArrow());
    }
}
