#!/usr/bin/env python3
from KudoCop import SimulationDat
from KudoCop import SimulationDats
from KudoCop.utils import analyze_file_type
import numpy as np
import pandas as pd
import argparse
from pprint import pprint

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='X-Y-Zの数をカウントする.')
    parser.add_argument('-s', '--skip_num', default=1, type=int,
                        help='何個おきにファイルを読み込むか')
    parser.add_argument('-f', '--file_name', default=None, type=str,
                        help='読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む')
    parser.add_argument('-c', '--cut_off', default=0.5, type=float,
                        help='カットオフ')
    parser.add_argument('-b', '--bond_type', default='dumpbond', type=str, choices=['dumppos', 'dumpbond'],
                        help='結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る')
    parser.add_argument('-p', '--para_file_name', default='para.rd', type=str,
                        help='para.rdのファイル名')

    parser.add_argument('--x_min', default=None, type=float,
                        help='min x')
    parser.add_argument('--y_min', default=None, type=float,
                        help='min y')
    parser.add_argument('--z_min', default=None, type=float,
                        help='min z')
    parser.add_argument('--x_max', default=None, type=float,
                        help='max x')
    parser.add_argument('--y_max', default=None, type=float,
                        help='max y')
    parser.add_argument('--z_max', default=None, type=float,
                        help='max z')

    args = parser.parse_args()

    if args.file_name is None:
        # SimulationDats
        def condition(sdats, step_idx):
            target_atoms = np.array([True] * sdats.get_total_atoms())
            if args.x_min is not None:
                target_atoms &= (sdats.atoms[step_idx]['x'] >= args.x_min)
            if args.y_min is not None:
                target_atoms &= (sdats.atoms[step_idx]['y'] >= args.y_min)
            if args.z_min is not None:
                target_atoms &= (sdats.atoms[step_idx]['z'] >= args.z_min)
            if args.x_max is not None:
                target_atoms &= (sdats.atoms[step_idx]['x'] <= args.x_max)
            if args.y_max is not None:
                target_atoms &= (sdats.atoms[step_idx]['y'] <= args.y_max)
            if args.z_max is not None:
                target_atoms &= (sdats.atoms[step_idx]['z'] <= args.z_max)
            return target_atoms

        sdats = SimulationDats(skip_num=args.skip_num,
                               para_file_name=args.para_file_name)
        df_count_triplets = sdats.count_triplets(bond_type=args.bond_type,   cut_off=args.cut_off,
                                                 condition=condition)
        print(f'bond_type: {args.bond_type}, cut_off: {args.cut_off}')
        print(f'condition: {args.x_min} <= x <= {args.x_max}')
        print(f'condition: {args.y_min} <= y <= {args.y_max}')
        print(f'condition: {args.z_min} <= z <= {args.z_max}')

        print(df_count_triplets)
    else:
        # SimulationDat
        sdat = SimulationDat()
        sdat.import_para(args.para_file_name)
        import_file_type = analyze_file_type(args.file_name)
        sdat.import_file(args.file_name, import_file_type)
        if import_file_type == 'dumppos':
            try:
                sdat.import_dumpbond(args.file_name.replace('pos', 'bond'))
            except:
                pass

        def condition(sdats):
            target_atoms = np.array([True] * sdats.get_total_atoms())
            if args.x_min is not None:
                target_atoms &= (sdats.atoms['x'] >= args.x_min)
            if args.y_min is not None:
                target_atoms &= (sdats.atoms['y'] >= args.y_min)
            if args.z_min is not None:
                target_atoms &= (sdats.atoms['z'] >= args.z_min)
            if args.x_max is not None:
                target_atoms &= (sdats.atoms['x'] <= args.x_max)
            if args.y_max is not None:
                target_atoms &= (sdats.atoms['y'] <= args.y_max)
            if args.z_max is not None:
                target_atoms &= (sdats.atoms['z'] <= args.z_max)
            return target_atoms

        dict_count_triplets = sdat.count_triplets(
            bond_type=args.bond_type, cut_off=args.cut_off, condition=condition)
        print(f'bond_type: {args.bond_type}, cut_off: {args.cut_off}')
        print(f'condition: {args.x_min} <= x <= {args.x_max}')
        print(f'condition: {args.y_min} <= y <= {args.y_max}')
        print(f'condition: {args.z_min} <= z <= {args.z_max}')
        pprint(dict_count_triplets)
