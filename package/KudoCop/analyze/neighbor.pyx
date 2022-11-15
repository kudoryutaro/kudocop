# distutils: language = c++

import numpy as np
from libc.stdlib cimport malloc
from libcpp.vector cimport vector


cdef struct atom_t:
    int ind
    int pos
    double x[3]


cdef void make_cbodies(list pos, atom_t *catoms):
    cdef atom_t *catom
    for i, p in enumerate(pos):
        catom = &catoms[i]
        catom.ind = i + 1
        catom.x[0], catom.x[1], catom.x[2] = p

cdef void min_length(double x[3], double cell[3]):
    cdef:
        int dim
    for dim in range(3):
        if x[dim] >= cell[dim] * 0.5:
            x[dim] -= cell[dim]
        if x[dim] <= - cell[dim] * 0.5:
            x[dim] += cell[dim]

cdef void NeigborList(int ix[3], vector[vector[int]] &append_mesh, int mx[3], vector[int] &search_list):
    cdef:
        int dim, index

    for dim in range(3):
        if ix[dim] < 0:
            ix[dim] += mx[dim]
        if ix[dim] >= mx[dim]:
            ix[dim] -= mx[dim]

    index = ix[0] + ix[1] * mx[0] + ix[2] * mx[0] * mx[1]
    search_list.insert(search_list.end(), append_mesh[index].begin(), append_mesh[index].end())

cdef void search(atom_t *catoms, double cell[3], int mesh_id, vector[vector[int]] &append_mesh, list neighbor_list, int mx[3], double CL2):
    cdef:
        int dim, i, iid, jid
        int my_len, search_len
        int ix[3]
        int add_ix[3]
        double dx[3]
        vector[int] search_list = []

    ix[0] = mesh_id % mx[0]
    ix[1] = (mesh_id / mx[0]) % mx[1]
    ix[2] = (mesh_id / mx[0] / mx[1])

    search_list.clear()
    my_len = len(append_mesh[mesh_id])

    add_index = (
        [0,0,0],
        [1,0,0],
        [-1,1,0],
        [0,1,0],
        [1,1,0],
        [-1,0,1],
        [0,0,1],
        [1,0,1],
        [-1,-1,1],
        [0,-1,1],
        [1,-1,1],
        [-1,1,1],
        [0,1,1],
        [1,1,1],
    )

    for add in add_index:
        for dim in range(3):
            add_ix[dim] = ix[dim] + add[dim]
        NeigborList(add_ix, append_mesh, mx, search_list)

    search_len = len(search_list)
    for i in range(my_len):
        iid = search_list[i]
        for j in range(i+1,search_len):
            jid = search_list[j]
            for dim in range(3):
                dx[dim] = catoms[jid].x[dim] - catoms[iid].x[dim]
            min_length(dx, cell)
            if dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2] <= CL2:
                neighbor_list[iid].append(jid)
                neighbor_list[jid].append(iid)

def make_neighbor(data , CL_):
    cdef:
        double cell[3]
        double msx[3]
        double imx[3]
        double CL = CL_
        int mx[3], ix[3]
        int nm, total = data.total_particle
        int i, j, dim, index, summ, ind_pos
        atom_t *catoms =  <atom_t *>malloc(total * sizeof(atom_t))
        vector[vector[int]] append_mesh


    pos = data.particles['pos'].tolist()
    make_cbodies(pos, catoms)

    nm = 1
    for dim in range(3):
        cell[dim] = data.newcell[dim]
        mx[dim] = int(cell[dim]/CL)
        nm *= mx[dim]

        msx[dim] = cell[dim] / mx[dim]
        imx[dim] = 1.0 / msx[dim]

    append_mesh.resize(nm)

    for i in range(total):
        for dim in range(3):
            ix[dim] = int(catoms[i].x[dim] * imx[dim])
            if ix[dim] < 0:
                ix[dim] += mx[dim]
            if ix[dim] >= mx[dim]:
                ix[dim] -= mx[dim]

        index = ix[0] + ix[1] * mx[0] + ix[2] * mx[0] * mx[1]
        append_mesh[index].push_back(i)

    neighbor_list = [[] for _ in range(total)]

    for i in range(nm):
        search(catoms, cell, i, append_mesh, neighbor_list, mx, CL*CL)


    return neighbor_list

cdef void search_pairwise_cutoff(atom_t *catoms, double cell[3], int mesh_id, vector[vector[int]] &append_mesh, list neighbor_list, int mx[3], list atom_types, list pair_cut_off2, int max_atom_type):
    cdef:
        int dim, i, iid, jid, atom_i_type, atom_j_type
        int my_len, search_len
        int ix[3]
        int add_ix[3]
        double dx[3]
        vector[int] search_list = []

    ix[0] = mesh_id % mx[0]
    ix[1] = (mesh_id / mx[0]) % mx[1]
    ix[2] = (mesh_id / mx[0] / mx[1])

    search_list.clear()
    my_len = len(append_mesh[mesh_id])

    add_index = (
        [0,0,0],
        [1,0,0],
        [-1,1,0],
        [0,1,0],
        [1,1,0],
        [-1,0,1],
        [0,0,1],
        [1,0,1],
        [-1,-1,1],
        [0,-1,1],
        [1,-1,1],
        [-1,1,1],
        [0,1,1],
        [1,1,1],
    )

    for add in add_index:
        for dim in range(3):
            add_ix[dim] = ix[dim] + add[dim]
        NeigborList(add_ix, append_mesh, mx, search_list)

    search_len = len(search_list)
    for i in range(my_len):
        iid = search_list[i]
        atom_i_type = atom_types[iid]
        for j in range(i+1,search_len):
            jid = search_list[j]
            atom_j_type = atom_types[jid]
            for dim in range(3):
                dx[dim] = catoms[jid].x[dim] - catoms[iid].x[dim]
            min_length(dx, cell)
            if dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2] <= pair_cut_off2[atom_i_type*max_atom_type+ atom_j_type]:
                neighbor_list[iid].append(jid)
                neighbor_list[jid].append(iid)

def make_neighbor_pairwise_cutoff(data , CL_,list atom_types,  list pair_cut_off2, int max_atom_type):
    cdef:
        double cell[3]
        double msx[3]
        double imx[3]
        double CL = CL_
        int mx[3], ix[3]
        int nm, total = data.total_particle
        int i, j, dim, index, summ, ind_pos
        atom_t *catoms =  <atom_t *>malloc(total * sizeof(atom_t))
        vector[vector[int]] append_mesh


    pos = data.particles['pos'].tolist()
    make_cbodies(pos, catoms)

    nm = 1
    for dim in range(3):
        cell[dim] = data.newcell[dim]
        mx[dim] = int(cell[dim]/CL)
        nm *= mx[dim]

        msx[dim] = cell[dim] / mx[dim]
        imx[dim] = 1.0 / msx[dim]

    append_mesh.resize(nm)

    for i in range(total):
        for dim in range(3):
            ix[dim] = int(catoms[i].x[dim] * imx[dim])
            if ix[dim] < 0:
                ix[dim] += mx[dim]
            if ix[dim] >= mx[dim]:
                ix[dim] -= mx[dim]

        index = ix[0] + ix[1] * mx[0] + ix[2] * mx[0] * mx[1]
        append_mesh[index].push_back(i)

    neighbor_list = [[] for _ in range(total)]

    for i in range(nm):
        search_pairwise_cutoff(catoms, cell, i, append_mesh, neighbor_list, mx, atom_types, pair_cut_off2, max_atom_type)


    return neighbor_list

#####for debug
def make_brute_neighbor(data, CL):

    cdef:
        int iid, jid, dim, total=data.total_particle
        double dx[3]
        double cell[3]
        atom_t *catoms =  <atom_t *>malloc(total * sizeof(atom_t))

    pos = data.particles['pos'].tolist()
    make_cbodies(pos, catoms)

    for dim in range(3):
        cell[dim] = data.newcell[dim]

    neighbor_list = [[] for _ in range(total)]

    for iid in range(total):
        for jid in range(iid+1, total):
            for dim in range(3):
                dx[dim] = catoms[jid].x[dim] - catoms[iid].x[dim]
            min_length(dx, cell)
            if dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2] <= CL * CL:
                neighbor_list[iid].append(jid)
                neighbor_list[jid].append(iid)

    return neighbor_list



cdef void search_relative_coordinates(atom_t *catoms, double cell[3], int mesh_id, vector[vector[int]] &append_mesh, list neighbor_list, relative_coordinates_list, int mx[3], double CL2):
    cdef:
        int dim, i, iid, jid
        int my_len, search_len
        int ix[3]
        int add_ix[3]
        double dx[3]
        vector[int] search_list = []

    ix[0] = mesh_id % mx[0]
    ix[1] = (mesh_id / mx[0]) % mx[1]
    ix[2] = (mesh_id / mx[0] / mx[1])

    search_list.clear()
    my_len = len(append_mesh[mesh_id])

    add_index = (
        [0,0,0],
        [1,0,0],
        [-1,1,0],
        [0,1,0],
        [1,1,0],
        [-1,0,1],
        [0,0,1],
        [1,0,1],
        [-1,-1,1],
        [0,-1,1],
        [1,-1,1],
        [-1,1,1],
        [0,1,1],
        [1,1,1],
    )

    for add in add_index:
        for dim in range(3):
            add_ix[dim] = ix[dim] + add[dim]
        NeigborList(add_ix, append_mesh, mx, search_list)

    search_len = len(search_list)
    for i in range(my_len):
        iid = search_list[i]
        for j in range(i+1,search_len):
            jid = search_list[j]
            for dim in range(3):
                dx[dim] = catoms[jid].x[dim] - catoms[iid].x[dim]
            min_length(dx, cell)
            if dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2] <= CL2:
                neighbor_list[iid].append(jid)
                neighbor_list[jid].append(iid)
                relative_coordinates_list[iid].append(dx)
                for dim in range(3):
                    dx[dim] *= 1
                relative_coordinates_list[jid].append(dx)
                

def make_relative_coordinates(data , CL_):
    cdef:
        double cell[3]
        double msx[3]
        double imx[3]
        double CL = CL_
        int mx[3], ix[3]
        int nm, total = data.total_particle
        int i, j, dim, index, summ, ind_pos
        atom_t *catoms =  <atom_t *>malloc(total * sizeof(atom_t))
        vector[vector[int]] append_mesh


    pos = data.particles['pos'].tolist()
    make_cbodies(pos, catoms)

    nm = 1
    for dim in range(3):
        cell[dim] = data.newcell[dim]
        mx[dim] = int(cell[dim]/CL)
        nm *= mx[dim]

        msx[dim] = cell[dim] / mx[dim]
        imx[dim] = 1.0 / msx[dim]

    append_mesh.resize(nm)

    for i in range(total):
        for dim in range(3):
            ix[dim] = int(catoms[i].x[dim] * imx[dim])
            if ix[dim] < 0:
                ix[dim] += mx[dim]
            if ix[dim] >= mx[dim]:
                ix[dim] -= mx[dim]

        index = ix[0] + ix[1] * mx[0] + ix[2] * mx[0] * mx[1]
        append_mesh[index].push_back(i)

    neighbor_list = [[] for _ in range(total)]
    relative_coordinates_list = [[] for _ in range(total)]

    for i in range(nm):
        search_relative_coordinates(catoms, cell, i, append_mesh, neighbor_list, relative_coordinates_list, mx, CL*CL)

    return neighbor_list, relative_coordinates_list