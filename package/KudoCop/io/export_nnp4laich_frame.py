import numpy as np
import pandas as pd
from pathlib import Path
import os
from ase.neighborlist import neighbor_list as make_neighbor_list
from ase import Atoms

class ExportNNP4LaichFrame():
    def __init__(self):
        pass

    # def export_nnp4laich_frame(self, nnp4laich_frames_dir:Path, frame_dir_name:str, cut_off:float, sel:int):
    #     """nnp4laichで読み込むことの出来る形式のデータセットを作成する.
    #     Parameters
    #     ----------
    #         nnp4laich_frames_dir : Path
    #             framesのディレクトリのパス
    #         frame_dir_name : str
    #             frames内のframeの名前
    #         cut_off : float
    #             カットオフ半径
    #         sel : int
    #             最大隣接原子数
    #     Note
    #     ----
    #         frameの中にframeが入っている。frameの中にatom_types, coord, forceなどが入っている。
    #         1つのフレームは一つのframeに対応する
    #             frames_dir
    #             ├── frame.0
    #             │   ├── atom_types
    #             │   ├── coord
    #             │   └── force
    #             ├── frame.1
    #             │   ├── atom_types
    #             │   ├── coord
    #             │   └── force
    #             └── ...
    #     """
    #     assert hasattr(self, 'atoms_calc'), 'dmol3_calcで先に第一原理計算してください'
    #     nnp4laich_frames_dir = Path(nnp4laich_frames_dir)

    #     os.makedirs(nnp4laich_frames_dir, exist_ok=True)
    #     nnp4laich_frame_dir = nnp4laich_frames_dir / frame_dir_name
    #     os.makedirs(nnp4laich_frame_dir, exist_ok=True)

    #     atom_types = self.atoms['type'].values
    #     coord = self.atoms[['x', 'y', 'z']].values
    #     cell = np.array(self.cell)
    #     force = self.atoms_calc.get_forces()
    #     total_potential_energy = np.array(self.atoms_calc.get_total_energy())
    #     neighbor_list, neighbor_atom_types, shift_vectors = make_neighbor_ase(
    #         coord=coord,
    #         atom_types=atom_types,
    #         cut_off=cut_off, 
    #         cell=cell)

    #     atom_i_types = make_atom_i_types(sel=sel, atom_types=atom_types)
    #     atom_j_types = make_atom_j_types(sel=sel, atom_types=atom_types, neighbor_list=neighbor_list)
    #     atom_i_idxs = make_atom_i_idxs(sel=sel, coord=coord)
    #     atom_j_idxs = make_atom_j_idxs(sel, coord, neighbor_list)
    #     shift = make_shift(sel, coord, cell, shift_vectors, neighbor_list)

    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=total_potential_energy, 
    #                                     np_array_name='total_potential_energy')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=np.array(cut_off), 
    #                                     np_array_name='cut_off')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=np.array(sel), 
    #                                     np_array_name='sel')

    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=atom_types, 
    #                                     np_array_name='atom_types')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=coord, 
    #                                     np_array_name='coord')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=cell, 
    #                                     np_array_name='cell')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=force, 
    #                                     np_array_name='force')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=atom_i_types, 
    #                                     np_array_name='atom_i_types')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=atom_j_types, 
    #                                     np_array_name='atom_j_types')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=atom_i_idxs, 
    #                                     np_array_name='atom_i_idxs')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=atom_j_idxs, 
    #                                     np_array_name='atom_j_idxs')
        
    #     export_nnp4laich_frame_npy(nnp4laich_frame_dir=nnp4laich_frame_dir, 
    #                                     np_array=shift, 
    #                                     np_array_name='shift')
        
        


# def export_nnp4laich_frame_npy(nnp4laich_frame_dir:Path, np_array:np.array, np_array_name:str):
#     """nnp4laich_frame_dirにnp_arrayをnp_array_name.npyという名前で保存する
#     Parameter
#     ---------
#     nnp4laich_frame_dir : Path
#         atom_types.npyを作成するディレクトリのパス
#     np_array : np.array
#         保存されるarray
#     np_array_name : str
#         保存する名前
#     """
#     nnp4laich_frame_dir = Path(nnp4laich_frame_dir)
#     np.save(nnp4laich_frame_dir / np_array_name, np_array)



def make_neighbor_ase(coord:np.array, atom_types:np.array,  cut_off:float, cell:np.array):
    ase_atoms = Atoms(
        positions=coord,
        cell=cell,
        pbc=[1, 1, 1]
        )

    i_idx, j_idx, cell_offsets= make_neighbor_list(
        'ijS', ase_atoms, cutoff=cut_off, self_interaction=False
    )
    neighbor_list = [[] for _ in range(coord.shape[0])]
    neighbor_atom_types = [[] for _ in range(coord.shape[0])]
    shift_vectors = [[] for _ in range(coord.shape[0])]
    for idx in range(len(i_idx)):
        atom_i_idx = i_idx[idx]
        atom_j_idx = j_idx[idx]
        neighbor_list[atom_i_idx].append(atom_j_idx)
        neighbor_atom_types[atom_i_idx].append(atom_types[atom_j_idx])
        shift_vectors[atom_i_idx].append(cell_offsets[idx])
    
    return neighbor_list, neighbor_atom_types, shift_vectors

def make_atom_i_types(sel:int, atom_types:list):
    atom_i_types = [0] *(len(atom_types) * sel)
    for atom_i_idx in range(len(atom_types)):
        for i in range(sel):
            atom_i_types[atom_i_idx*sel + i] = atom_types[atom_i_idx]
    return atom_i_types

def make_atom_j_types(sel:int, atom_types:list, neighbor_list:list):
    atom_j_types = [0] * (len(atom_types) * sel)
    for atom_i_idx in range(len(atom_types)):
        for j in range(sel):
            if len(neighbor_list[atom_i_idx]) <= j:
                break
            atom_j_type = atom_types[neighbor_list[atom_i_idx][j]]
            atom_j_types[atom_i_idx*sel + j] = atom_j_type
    return atom_j_types

def make_atom_i_idxs(sel, coord):
    total_atoms = coord.shape[0]
    atom_i_idxs = [total_atoms] * (sel * total_atoms)
    for atom_i_idx in range(total_atoms):
        for j in range(sel):
            atom_i_idxs[atom_i_idx*sel + j] = atom_i_idx
    return atom_i_idxs

def make_atom_j_idxs(sel, coord, neighbor_list):
    total_atoms = coord.shape[0]
    atom_j_idxs = [total_atoms] * (sel * total_atoms)
    for atom_i_idx in range(total_atoms):
        for j in range(sel):
            if len(neighbor_list[atom_i_idx]) <= j:
                break
            atom_j_idxs[atom_i_idx*sel + j] = neighbor_list[atom_i_idx][j]
    return atom_j_idxs

def make_shift(sel, coord, cell, shift_vectors, neighbor_list):
    total_atoms = coord.shape[0]
    shift = [[0, 0, 0] for _ in range(sel * total_atoms)]
    for atom_i_idx in range(total_atoms):
        for j in range(sel):
            if len(neighbor_list[atom_i_idx]) <= j:
                break
            shift[atom_i_idx*sel+ j][0] = cell[0] * shift_vectors[atom_i_idx][j][0]
            shift[atom_i_idx*sel+ j][1] = cell[1] * shift_vectors[atom_i_idx][j][1]
            shift[atom_i_idx*sel+ j][2] = cell[2] * shift_vectors[atom_i_idx][j][2]
    return shift



