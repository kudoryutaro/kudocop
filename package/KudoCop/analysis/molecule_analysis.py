import sys
from collections import deque


def count_mols(sdat, cut_off, lower_mol_limit, upper_mol_limit) -> dict:
    mol_counter = dict()
    sdat.get_connect_list(cut_off)

    visited = [0] * sdat.get_total_atoms()

    # number of types 種類数
    type_num_max = max(sdat.atom_type_set)

    sdat_atoms_type = sdat.atoms['type'].values

    for start_atom_idx in range(sdat.get_total_atoms()):
        if visited[start_atom_idx]:
            continue

        current_mol_counter = [0] * (type_num_max + 1)

        # 幅優先探索
        que = deque([start_atom_idx])
        while que:
            current_atom_idx = que.popleft()
            if visited[current_atom_idx]:
                continue
            visited[current_atom_idx] = 1
            current_mol_counter[sdat_atoms_type[current_atom_idx]] += 1
            for next_atom_idx in sdat.connect_list[current_atom_idx]:
                if visited[next_atom_idx]:
                    continue
                que.append(next_atom_idx)

        if lower_mol_limit <= sum(current_mol_counter) <= upper_mol_limit:
            current_mol_counter_tuple = tuple(current_mol_counter)
            if current_mol_counter_tuple not in mol_counter:
                mol_counter[current_mol_counter_tuple] = 0
            mol_counter[current_mol_counter_tuple] += 1

    return mol_counter


def get_atom_idx_from_mol(sdat, cut_off, target_mol):
    atom_idx_from_mol = []
    sdat.get_connect_list(cut_off)
    visited = [0] * sdat.get_total_atoms()

    # number of types 種類数
    type_num_max = max(sdat.atom_type_set)
    if len(target_mol) != type_num_max + 1:
        print('target_mol\'s length does not match')
        sys.exit(-1)

    sdat_atoms_type = sdat.atoms['type'].values
    for start_atom_idx in range(sdat.get_total_atoms()):
        if visited[start_atom_idx]:
            continue

        current_mol_counter = [0] * (type_num_max + 1)
        current_mol = []

        que = deque([start_atom_idx])

        while que:
            current_atom_idx = que.popleft()
            if visited[current_atom_idx]:
                continue
            visited[current_atom_idx] = 1
            current_mol.append(current_atom_idx)
            current_mol_counter[sdat_atoms_type[current_atom_idx]] += 1

            for next_atom_idx in sdat.connect_list[current_atom_idx]:
                if visited[next_atom_idx]:
                    continue
                que.append(next_atom_idx)

        if tuple(current_mol_counter) == target_mol:
            atom_idx_from_mol.append(current_mol)

    return atom_idx_from_mol
