import numpy as np
import pandas as pd


class WriteInput():
    def __init__(self):
        pass

    def write_file(self, sdat, ofn: str) -> None:
        # header_line
        header_line = []
        header_line.append(f"#cellx {0.0}  {sdat.cell[0]}\n")
        header_line.append(f"#celly {0.0}  {sdat.cell[1]}\n")
        header_line.append(f"#cellz {0.0}  {sdat.cell[2]}\n")

        header_line.append("\n")

        atom_type_set = set(sdat.atoms['type'])

        header_line.append(f"#masses {len(atom_type_set)}\n")
        for atom_type in atom_type_set:
            header_line.append(
                f"{atom_type} {sdat.atom_type_to_mass[atom_type]}\n")

        header_line.append("\n")

        for fix in sdat.fix_info:
            header_line.append(fix + '\n')

        for move in sdat.move_info:
            header_line.append(move + '\n')

        for press in sdat.press_info:
            header_line.append(press + '\n')

        for sumforce in sdat.sumforce_info:
            header_line.append(sumforce + '\n')

        for thermofree in sdat.thermofree_info:
            header_line.append(thermofree + '\n')

        for wall in sdat.wall_info:
            header_line.append(wall + '\n')

        header_line.append("")

        header_line.append(f"#atoms {sdat.get_total_atoms()}\n")

        with open(ofn, 'w') as ofs:
            ofs.writelines(header_line)

        # id, type, mask, pos, velo
        out_columns = ['type', 'mask', 'x', 'y', 'z']

        if 'mask' not in sdat.atoms:
            sdat.atoms['mask'] = np.zeros(sdat.get_total_atoms(),
                                          dtype=int)
        if 'vx' in sdat.atoms and 'vy' in sdat.atoms and 'vz' in sdat.atoms:
            out_columns.append('vx')
            out_columns.append('vy')
            out_columns.append('vz')

        # 1-indexed
        sdat.atoms.index = sdat.atoms.index + 1
        sdat.atoms.to_csv(ofn, columns=out_columns, mode='a', header=False,
                          sep=' ', float_format='%.6f')
        # 0-indexed
        sdat.atoms.index = sdat.atoms.index - 1

        # connect_list
        if sdat.flagconnect:
            with open(ofn, 'a') as ofs:
                print(f"#connect {len(sdat.connect_list)}", file=ofs)
                for ind, connect in enumerate(sdat.connect_list):
                    tmp_list = [str(ind+1), str(len(connect))]
                    for target in connect:
                        tmp_list.append(str(target+1))
                    print(" ".join(tmp_list), file=ofs)
