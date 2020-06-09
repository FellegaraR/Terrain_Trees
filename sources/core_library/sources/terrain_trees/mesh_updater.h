#ifndef MESH_UPDATER_H
#define MESH_UPDATER_H

#include "basic_types/mesh.h"

/**
 * @brief The Mesh_Updater class represents an interface for updating a simplicial (or CP) mesh.
 *
 */
class Mesh_Updater
{
public:
    Mesh_Updater() {}

    /**
     * @brief A public method that compresses the vertices array of the mesh
     * The procedure removes the vertices flagged as deleted during a simplification procedure
     * and returns the new position indexes of the vertices not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     * @param surviving_vertices a vector of integers containing the non-deleted vertices
     *
     * NOTA: the size of the new_v_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     */
    template<class C, class T> void clean_vertices_array(Mesh &mesh, ivect &new_v_positions, ivect &surviving_vertices);
    /**
     * @brief A public method that compresses the triangles arrays of the mesh
     * The procedure removes the triangles flagged as deleted during a simplification procedure
     * and returns, for each top k-cell array, if it has been completely erased.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     * @param new_top_positions a vector of vectors that, at the end of the procedure, contains the new position indexes of the triangles
     * @param simplification_counters a vector referring to the number of triangles removed during a simplification procedure
     * @return a bit-vector sized as the number of triangles encoded, returning the triangles type that have been completely erased by a simplification procedure
     *
     * NOTA: the size of the new_v_positions and new_top_positions arrays is equal to the size of the original arrays of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     *
     * NOTA2: during the procedure are also updated the boundary relations of the triangles
     */
    template<class C, class T> boost::dynamic_bitset<> update_and_clean_triangles_arrays(Mesh &mesh, ivect &new_v_positions,
                                                                                         vector<ivect> &new_triangle_positions, ivect &simplification_counters);

    /**
     * @brief A public method that updates the boundary relations of the triangles, by assigning the new position indexes of the corresponding vertices,
     * and resorting the vertices array in the mesh in order to be coherent to the new vertices ordering (calls reorder_vertices_array)
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     */
    template<class C, class T> void update_triangles_boundary(Mesh &mesh, ivect &new_v_positions);

    /**
     * @brief A private method that update the vertices array of the mesh by assigning the new position indexes.
     * The procedure execute the swaps inline an thus it does not require extra storage.
     * The number of swap operation is exactly |V|, where |V| is the number of vertices of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     */
    template<class C, class T> void reorder_vertices_array(Mesh &mesh, ivect &new_v_positions);
    /**
     * @brief A public method that update the triangles arrays of the mesh by assigning the new position indexes.
     * The procedure execute the swaps inline an thus it does not require extra storage.
     * The number of swap operation is exactly |T|, where |T| is the number of triangles of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_top_positions a vector of vectors that, at the end of the procedure, contains the new position indexes of the triangles
     */
    template<class C, class T> void reorder_triangles_array(Mesh &mesh, vector<ivect> &new_triangles_positions);

private:
    /**
     * @brief A public method that compresses the top d-cells arrays of the mesh.
     * The procedure removes the top d-cells flagged as deleted during a simplification procedure
     * and returns if the corresponding array has been completely erased.
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     * @param new_top_positions a vector that, at the end of the procedure, contains the new position indexes of the top d-cells
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param counter an integer representing the number of the top d-cells deleted during the simplification procedure
     * @return a boolean, true if all the top d-cells have been deleted, false otherwise
     */
    template<class C, class T> bool update_and_clean_triangles_array(Mesh &mesh, ivect &new_v_positions, ivect &new_triangle_positions, int d, int counter);
    /**
     * @brief A public method that compresses the top d-cells array of the mesh.
     * The procedure removes the top d-cell deleted during a simplification procedure
     * and returns the new position indexes of the top d-cells not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param new_triangle_positions a vector of integers that at the end of the procedure contains the new position indexes of the top d-cells
     *
     * NOTA: the size of the new_triangle_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     */
    template<class C, class T> void clean_triangles_array(Mesh &mesh, int d, ivect &new_triangle_positions);
    /**
     * @brief A public procedure that updates the boudary relation for all the top d-cells of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     */
    template<class C, class T> void update_triangles_boundary(Mesh &mesh, int d, ivect &new_v_positions);
    /**
     * @brief A public procedure that coherently reorder all the top d-cells of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param new_top_positions a vector that, at the end of the procedure, contains the new position indexes of the top d-cells
     */
    template<class C, class T> void reorder_triangles_array(Mesh &mesh, int d, ivect new_triangles_positions);
};

template<class C, class T> void Mesh_Updater::clean_vertices_array(Mesh &mesh, ivect &new_v_positions, ivect &surviving_vertices)
{
    int v_counter = 1; // v_counter keeps the new vertex indexing while v is needed to set the new index
    new_v_positions.assign(mesh.get_vertices_num(),-1);

    vector<Vertex<C> > old_list = mesh.get_vertices_array();
    mesh.reset_vertices_vector();
    mesh.reserve_vertices_space(surviving_vertices.size());

    for(ivect_iter it=surviving_vertices.begin(); it!=surviving_vertices.end(); ++it)
    {
        mesh.add_vertex(old_list[*it-1]);
        new_v_positions[*it-1] = v_counter;
        v_counter++;
    }

    old_list.clear();
}

template<class C, class T> boost::dynamic_bitset<> Mesh_Updater::update_and_clean_triangles_arrays(Mesh &mesh, ivect &new_v_positions, vector<ivect > &new_top_positions,
                                                                                                   ivect &simplification_counters)
{
    boost::dynamic_bitset<> all_deleted = boost::dynamic_bitset<>(mesh.get_top_cells_types());
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        ivect triangle_pos;
        all_deleted[d] = update_and_clean_triangles_array(mesh,new_v_positions,triangle_pos,d,simplification_counters[d]);
        new_top_positions.push_back(triangle_pos);
    }
    return all_deleted;
}

template<class C, class T> bool Mesh_Updater::update_and_clean_triangles_array(Mesh &mesh, ivect &new_v_positions, ivect &new_top_positions, int d, int counter)
{
    // step 0: if we have removed all the top d-cells, then, simply clear the corresponding array
    if(counter==mesh.get_triangles_num())
    {
        mesh.reset_triangles();
        return true;
    }
    else
    {
        clean_triangles_array(mesh,d,new_top_positions);
        update_triangles_boundary(mesh,d,new_v_positions);
    }
    return false;
}

template<class C, class T> void Mesh_Updater::clean_triangles_array(Mesh &mesh, int d, ivect &new_triangle_positions)
{
    int t_counter = 1, t = 0; // t_counter keeps the new top indexing while t is needed to set the new index

    // (1) init the triangles position indices array
    new_triangle_positions.assign(mesh.get_top_cells_num(d),-1);
    // (2) get the new position indices for the triangles that are not deleted
    for(typename Mesh::top_d_array_iter it=mesh.top_d_array_begin(d); it!=mesh.top_d_array_end(d); ++it)
    {
        if(!mesh.is_top_cell_removed(*it))
        {
            new_triangle_positions[t] = t_counter;
            t_counter++;
        }
        t++;
    }
    // (3) reorder (inline) the triangles arrays
    reorder_triangles_array(mesh,d,new_triangle_positions);
    // (4) delete the positions from t_counter to the end of the array
    mesh.resize_top_d_cell_array(d,t_counter-1);
}

template<class C, class T> void Mesh_Updater::update_triangles_boundary(Mesh &mesh, ivect &new_v_positions)
{
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        update_triangles_boundary(mesh,d,new_v_positions);
    }
}

template<class C, class T> void Mesh_Updater::update_triangles_boundary(Mesh &mesh, int d, ivect &new_v_positions)
{
    int t_id = 1;
    for(typename Mesh::top_d_array_iter it=mesh.top_d_array_begin(d); it!=mesh.top_d_array_end(d); ++it)
    {
        if(!mesh.is_top_cell_removed(*it))
        {
            for(int v=0; v<it->get_vertices_num(); v++)
            {
                if(new_v_positions[abs(it->TV(v))-1]==-1)
                {
                    cout<<"UPDATING A TOP-SIMPLEX WITH A DELETED VERTEX: "<<abs(it->TV(v))-1<<endl;
                    cout<<t_id<<" "<<*it<<endl;
                    int a; cin>>a;
                }
                it->setTV(v, new_v_positions[abs(it->TV(v))-1]);
            }
        }
        else
            cerr<<"[WARNING] update_top_d_cells_array does not removes deleted top d-cells"<<endl;

        t_id++;
    }
}

template<class C, class T> void Mesh_Updater::reorder_vertices_array(Mesh &mesh, ivect &new_v_positions)
{
    // (2) update global vertices array
    //     -1 identifies an already set vertex
    // ** the reordering is done inline the global vertex array ** ///
    for(int i=1;i<=mesh.get_vertices_num();i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_v_positions[j] == i || new_v_positions[j] < 0)
        {
            // mark the last entry visited...
            new_v_positions[j] = -1;
            continue;
        }

        while(new_v_positions[j] != i)
        {
            mesh.vertices_swap(i,new_v_positions[j]);
            int j_prime = new_v_positions[j] -1; // -1 is needed to avoid array lookup error
            new_v_positions[j] = -1;
            j = j_prime;
        }
        // mark the last entry visited...
        new_v_positions[j] = -1;
    }
}

template<class C, class T> void Mesh_Updater::reorder_triangles_array(Mesh &mesh, vector<ivect> &new_triangles_positions)
{
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        reorder_triangles_array(mesh,d,new_triangles_positions[d]);
    }
}

template<class C, class T> void Mesh_Updater::reorder_triangles_array(Mesh &mesh, int d, ivect new_triangles_positions)
{
    for(int i=1; i<=mesh.get_top_cells_num(d); i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_triangles_positions[j] == i || new_triangles_positions[j] < 0)
        {
            // mark the last entry visited...
            new_triangles_positions[j] = -1;
            continue;
        }

        while(new_triangles_positions[j] != i && new_triangles_positions[j] > 0)
        {
            mesh.triangles_swap(d,i,new_triangles_positions[j]);
            int j_prime = new_triangles_positions[j] -1; // -1 is needed to avoid array lookup error
            new_triangles_positions[j] = -1;
            j = j_prime;
        }

        // mark the last entry visited...
        new_triangles_positions[j] = -1;
    }
}

#endif // MESH_UPDATER_H
