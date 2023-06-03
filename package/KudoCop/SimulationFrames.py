import sys
import pandas as pd
import numpy as np
import os
from tqdm import tqdm, trange
import pathlib
from ase.neighborlist import neighbor_list 
import ase
import pickle

from .io.import_para import ImportPara
from .SimulationDat import SimulationDat


class SimulationFrames(
    ImportPara
):
    """シミュレーションしたデータを読み込み、書き込み、分析するためのクラス
    一度に複数のdump.posとdump.bondを読み込む

    Attributes
    ----------
    sdat : List[KudoCop.SimulationDat]
        原子の座標, type, 電荷, cellなどを含むSimulationDatクラス
    atom_symbol_to_type : dict
        原子のシンボルをkey, 原子のtypeをvalueとするdict
    atom_type_to_symbol : dict
        原子のtypeをkey, 原子のシンボルをvalueとするdict
    atom_type_to_mass : dict
        原子のtypeをkey, 原子の質量(g/mol)をvalueとするdict
    """
    def __init__(self):
        # variables for para.rd
        self.atom_symbol_to_type = None
        self.atom_type_to_symbol = None
        self.atom_type_to_mass = None

    def __len__(self):
        return len(self.step_nums)

    def __getitem__(self, key):
        return self.sdat[key]
    
    def import_dumpposes(self, dir_name: str=None, step_nums: int=None, skip_num: int=None):
        assert self.atom_symbol_to_type is not None, "import para first"
        assert self.atom_type_to_mass is not None, "import para first"
        assert self.atom_type_to_symbol is not None, "import para first"
        
        if dir_name is None:
            dir_name = os.getcwd()

        file_names_in_current_dir = os.listdir(dir_name)
        if step_nums is None:
            self.step_nums = []
            for file_name in file_names_in_current_dir:
                if len(file_name) >= 9 and file_name[:9] == 'dump.pos.':
                    self.step_nums.append(int(file_name[9:]))

        else:
            self.step_nums = []
            for step_num in step_nums:
                self.step_nums.append(step_num)
        self.step_nums.sort()
        if skip_num is not None:
            self.step_nums = self.step_nums[::skip_num]

        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }
        self.sdat = [SimulationDat() for _ in range(len(self.step_nums))]

        for step_idx, step_num in enumerate(tqdm(self.step_nums)):
            self.sdat[step_idx].atom_symbol_to_type = self.atom_symbol_to_type
            self.sdat[step_idx].atom_type_to_mass = self.atom_type_to_mass
            self.sdat[step_idx].atom_type_to_symbol = self.atom_type_to_symbol
            self.sdat[step_idx].import_dumppos(f'{dir_name}/dump.pos.{step_num}')
    
    def import_vasp(self, calc_directory: str):
        self.sdat = []
        calc_directory = pathlib.Path(calc_directory)
        with open(calc_directory / "POSCAR", "r") as f:
            for _ in range(5):
                f.readline()
            atom_symbol_list = list(f.readline().split())
            atom_type_counter = list(map(int, f.readline().split()))
            atom_types = []
            for atom_type_count, atom_symbol in zip(atom_type_counter, atom_symbol_list):
                for _ in range(atom_type_count):
                    atom_types.append(self.atom_symbol_to_type[atom_symbol])

        with open(calc_directory / "OUTCAR", "r") as f:
            lines = f.readlines()
            splines = list(map(lambda l:l.split(), lines))

        
        for line_idx, spline in enumerate(splines):
            if len(spline) == 0:
                continue
            if len(spline) == 3 and spline[0] == "POSITION" and spline[1] == "TOTAL-FORCE":
                position_line_idx = line_idx
                sdat = SimulationDat()
                sdat.atom_symbol_to_type = self.atom_symbol_to_type
                sdat.atom_type_to_mass = self.atom_type_to_mass
                sdat.atom_type_to_symbol = self.atom_type_to_symbol 
                cell_line_idx = position_line_idx
                while True:
                    if len(splines[cell_line_idx]) == 6 and splines[cell_line_idx][0] == "direct" \
                        and splines[cell_line_idx][1] == "lattice":
                        break
                    cell_line_idx -= 1
                potential_energy_idx = position_line_idx
                while True:
                    if len(splines[potential_energy_idx]) >= 4 and \
                        splines[potential_energy_idx][0] == "energy" and \
                        splines[potential_energy_idx][1] == "without" and \
                        splines[potential_energy_idx][2] == "entropy":
                        break
                    potential_energy_idx -= 1
                        
                
                sdat.cell = [None, None, None]
                sdat.cell[0] = float(splines[cell_line_idx+1][0])
                sdat.cell[1] = float(splines[cell_line_idx+1+1][1])
                sdat.cell[2] = float(splines[cell_line_idx+1+2][2])
                atoms_dict = {
                    'type':atom_types,
                    'x':[],
                    'y':[],
                    'z':[],
                    'fx':[],
                    'fy':[],
                    'fz':[],
                              }
                for atom_idx in range(len(atom_types)):
                    atoms_dict['x'].append(float(splines[position_line_idx+2+atom_idx][0]))
                    atoms_dict['y'].append(float(splines[position_line_idx+2+atom_idx][1]))
                    atoms_dict['z'].append(float(splines[position_line_idx+2+atom_idx][2]))
                    atoms_dict['fx'].append(float(splines[position_line_idx+2+atom_idx][3]))
                    atoms_dict['fy'].append(float(splines[position_line_idx+2+atom_idx][4]))
                    atoms_dict['fz'].append(float(splines[position_line_idx+2+atom_idx][5]))
                sdat.atoms = pd.DataFrame(atoms_dict)
                sdat.potential_energy = float(splines[potential_energy_idx][4])
                self.sdat.append(sdat)
        self.step_nums = list(range(1, len(self.sdat) + 1))
        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }


    def export_dumpposes(self, output_folder: str, out_columns=None) -> None:
        try:
            os.makedirs(output_folder)
        except FileExistsError:
            pass

        for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[exporting dumppos]')):
            ofn = f'{output_folder}/dump.pos.{step_num}'
            self.sdat[step_idx].export_dumppos(ofn, time_step=step_num, out_columns=out_columns)


    def export_lammps_dumpposes(self, ofn: str, out_columns=None) -> None:
        """lammps形式のdumpposを出力する
        """
        if out_columns is None:
            out_columns = ['type', 'x', 'y', 'z']

        with open(ofn, 'w') as f:
            f.write('')

        for step_idx in trange(len(self.step_nums), desc='[exporting lammps dumpposes]'):
            header = []
            header.append(f'ITEM: TIMESTEP\n')
            header.append(f'{self.step_nums[step_idx]}\n')
            header.append(f'ITEM: NUMBER OF ATOMS\n')
            header.append(f'{self.sdat[step_idx].get_total_atoms()}\n')
            header.append(f'ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
            header.append(f'0.0000000000000000e+00 {self.sdat[step_idx].cell[0]:.16e} 0.0000000000000000e+00\n')
            header.append(f'0.0000000000000000e+00 {self.sdat[step_idx].cell[1]:.16e} 0.0000000000000000e+00\n')
            header.append(f'0.0000000000000000e+00 {self.sdat[step_idx].cell[2]:.16e} 0.0000000000000000e+00\n')
            header.append(f'ITEM: ATOMS id {" ".join(out_columns)}\n')

            with open(ofn, 'a') as f:
                f.writelines(header)
            
            # 1-index
            self.sdat[step_idx].atoms.index += 1
            self.sdat[step_idx].atoms.to_csv(ofn, columns=out_columns, sep=' ', header=None, mode='a')
            # 0-index
            self.sdat[step_idx].atoms.index -= 1

    def export_allegro_frames(self, output_dir:str, output_file_name:str, cut_off:float):
        """allegro用のデータセットを保存する
        Parameter
        ---------
        output_dir : str
            作成するディレクトリのパス
        output_file_name : str  
            "output_file_name.pickle"という名前で保存する
        Note
        ----
            frames : List[Dict[str, np.array]]
            frames = [
                {
                "pos":np.array, 
                "force":np.array,
                "cell":np.array,
                "atom_types":np.array,
                "edge_index":np.array,
                "edge_cell_shift":np.array,
                "potential_energy":np.array
                },...
            ]
        """
        output_dir = pathlib.Path(output_dir)
        output_file_name = pathlib.Path(output_file_name + ".pickle")
        frames_path = output_dir / output_file_name

        os.makedirs(output_dir, exist_ok=True)
        frames = []
        for step_idx in range(len(self.step_nums)):
            data = {}
            data['cell'] = np.array(self.sdat[step_idx].cell, dtype=np.float32)
            data['pos'] = np.array(self.sdat[step_idx].atoms[['x', 'y', 'z']].values, dtype=np.float32)
            data['force'] = np.array(self.sdat[step_idx].atoms[['fx', 'fy', 'fz']].values, dtype=np.float32)
            data['atom_types'] = np.array(self.sdat[step_idx].atoms['type'].values)

            ase_atoms = ase.Atoms(positions=data['pos'], cell=data['cell'], pbc=[1, 1, 1])


            i_idx, j_idx, edge_cell_shift = neighbor_list(
                'ijS', ase_atoms, cutoff=cut_off, self_interaction=False
            )
            edge_index = [[], []]
            for atom_i_idx, atom_j_idx in zip(i_idx, j_idx):
                if atom_i_idx > atom_j_idx:
                    continue
                else:
                    edge_index[0].append(atom_i_idx)
                    edge_index[1].append(atom_j_idx)
            edge_index = np.array(edge_index)
            # data['edge_cell_shift'] = edge_cell_shift.astype(np.float32) * data['cell']
            data['edge_index'] = edge_index
            data['potential_energy'] = np.array(self.sdat[step_idx].potential_energy, dtype=np.float32)
            
            frames.append(data)
        
        with open(frames_path, "wb") as f:
            pickle.dump(frames, f)