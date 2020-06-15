#include "terrain_trees/mesh_updater.h"


void Mesh_Updater::clean_vertices_array(Mesh &mesh, ivect &new_v_positions, ivect &surviving_vertices)
{
    int v_counter = 1; // v_counter keeps the new vertex indexing while v is needed to set the new index
    new_v_positions.assign(mesh.get_vertices_num(),-1);

    vector<Vertex> old_list = mesh.get_vertices_array();
    mesh.reset_vertices();
    mesh.reserve_vertices_space(surviving_vertices.size());

    for(ivect_iter it=surviving_vertices.begin(); it!=surviving_vertices.end(); ++it)
    {
        mesh.add_vertex(old_list[*it-1]);
        new_v_positions[*it-1] = v_counter;
        v_counter++;
    }

    old_list.clear();
}

 boost::dynamic_bitset<> Mesh_Updater::update_and_clean_triangles_arrays(Mesh &mesh, ivect &new_v_positions, ivect &new_top_positions,
                                                                                                   itype simplification_counters)
{
    boost::dynamic_bitset<> all_deleted = boost::dynamic_bitset<>(0);
 
        //ivect triangle_pos;
        all_deleted[0] = update_and_clean_triangles_array(mesh,new_v_positions,new_top_positions,simplification_counters);
      //  new_top_positions=triangle_pos;
    
    return all_deleted;
}

bool Mesh_Updater::update_and_clean_triangles_array(Mesh &mesh, ivect &new_v_positions, ivect &new_top_positions, int counter)
{
    // step 0: if we have removed all the top d-cells, then, simply clear the corresponding array
    if(counter==mesh.get_triangles_num())
    {
        mesh.reset_triangles();
        return true;
    }
    else
    {
        clean_triangles_array(mesh,new_top_positions);
        update_triangles_boundary(mesh,new_v_positions);
    }
    return false;
}

 void Mesh_Updater::clean_triangles_array(Mesh &mesh,  ivect &new_triangle_positions)
{
    int t_counter = 1, t = 0; // t_counter keeps the new top indexing while t is needed to set the new index

    // (1) init the triangles position indices array
    new_triangle_positions.assign(mesh.get_triangles_num(),-1);
    // (2) get the new position indices for the triangles that are not deleted
    for(itype i=0; i<mesh.get_triangles_num(); ++i)  // See if need to be changed to iterator
    {
        if(!mesh.is_triangle_removed(i))
        {
            new_triangle_positions[t] = t_counter;
            t_counter++;
        }
        t++;
    }
    // (3) reorder (inline) the triangles arrays
    reorder_triangles_array(mesh,new_triangle_positions);
    // (4) delete the positions from t_counter to the end of the array
    mesh.resize_triangle_array(t_counter-1);
}



void Mesh_Updater::update_triangles_boundary(Mesh &mesh, ivect &new_v_positions)
{
    int t_id = 1;
    for(itype i=0; i<mesh.get_triangles_num(); ++i)
    {
        if(!mesh.is_triangle_removed(i))
        {
            for(int v=0; v<3; v++)
            {
                if(new_v_positions[abs(mesh.get_triangle(i).TV(v))-1]==-1)
                {
                    cout<<"UPDATING A TOP-SIMPLEX WITH A DELETED VERTEX: "<<abs(mesh.get_triangle(i).TV(v))-1<<endl;
                    cout<<t_id<<" "<<i<<endl;
                    int a; cin>>a;
                }
                mesh.get_triangle(i).setTV(v, new_v_positions[abs(mesh.get_triangle(i).TV(v))-1]);
            }
        }
        else
            cerr<<"[WARNING] update_top_d_cells_array does not removes deleted top d-cells"<<endl;

        t_id++;
    }
}

 void Mesh_Updater::reorder_vertices_array(Mesh &mesh, ivect &new_v_positions)
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



 void Mesh_Updater::reorder_triangles_array(Mesh &mesh,  ivect new_triangles_positions)
{
    for(int i=1; i<=mesh.get_triangles_num(); i++)
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
            mesh.triangles_swap(i,new_triangles_positions[j]);
            int j_prime = new_triangles_positions[j] -1; // -1 is needed to avoid array lookup error
            new_triangles_positions[j] = -1;
            j = j_prime;
        }

        // mark the last entry visited...
        new_triangles_positions[j] = -1;
    }
}