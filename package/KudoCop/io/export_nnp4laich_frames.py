import numpy as np
import pandas as pd
import numpy as np
from pathlib import Path
import os
from .export_nnp4laich_frame import *
from tqdm import trange
import pickle

class ExportNNP4LaichFrames():
    def __init__(self):
        pass

    def export_nnp4laich_frames(self, nnp4laich_frames_dir:Path, frame_name:str, cut_off:float, sel:int):
        """nnp4laichで読み込むことの出来る形式のデータセットを作成する.
        Parameters
        ----------
            nnp4laich_frames_dir : Path
                framesのディレクトリのパス
            frame_name : str
                保存するデータの名前は'{frame_name}.{step_idx}.pickle'で保存される
            cut_off : float
                カットオフ半径
            sel : int
                最大隣接原子数
        Note
        ----
            生成されるフォルダ
            nnp4laich_frames_dir
            ├── frame_name.0.pickle
            ├── frame_name.1.pickle
            ├── frame_name.2.pickle
            ├── frame_name.3.pickle
            └── ...
        """
        nnp4laich_frames_dir = Path(nnp4laich_frames_dir)

        if 'fx' not in self.atoms[0].columns:
            self.add_force_to_atoms()

        os.makedirs(nnp4laich_frames_dir, exist_ok=True)

        for step_idx in trange(len(self.atoms)):

            atom_types = self.atoms[step_idx]['type'].values
            coord = self.atoms[step_idx][['x', 'y', 'z']].values
            cell = np.array(self.cell)
            force = self.atoms[step_idx][['fx', 'fy', 'fz']].values
            total_potential_energy = self.potential_energy[step_idx]

            neighbor_list, neighbor_atom_types, shift_vectors = make_neighbor_ase(
                coord=coord,
                atom_types=atom_types,
                cut_off=cut_off, 
                cell=cell)

            atom_i_types = make_atom_i_types(sel=sel, atom_types=atom_types)
            atom_j_types = make_atom_j_types(sel=sel, atom_types=atom_types, neighbor_list=neighbor_list)
            atom_i_idxs = make_atom_i_idxs(sel=sel, coord=coord)
            atom_j_idxs = make_atom_j_idxs(sel, coord, neighbor_list)
            shift = make_shift(sel, coord, cell, shift_vectors, neighbor_list)

            data = {
                'atom_types':atom_types,
                'coord':coord,
                'force':force,
                'atom_i_types':atom_i_types,
                'atom_j_types':atom_j_types,
                'atom_i_idxs':atom_i_idxs,
                'atom_j_idxs':atom_j_idxs,
                'total_potential_energy':total_potential_energy,
                'shift':shift,
            }
            with open(nnp4laich_frames_dir / f'{frame_name}.{step_idx}.pickle', 'wb') as f:
                pickle.dump(data, f)

            
            