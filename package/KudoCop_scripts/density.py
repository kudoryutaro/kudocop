#!/usr/bin/env python3
from KudoCop import SimulationDat
from KudoCop import SimulationDats
from KudoCop.utils import analyze_file_type
import numpy as np
import pandas as pd
import argparse

pd.set_option('display.max_rows', None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='密度を計算する')
    parser.add_argument('-f', '--file_name', default=None, type=str,
                        help='読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む')
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
    assert args.file_name is not None, 'set file_name'
    sdat = SimulationDat()
    sdat.import_para(args.para_file_name)
    import_file_type = analyze_file_type(args.file_name)
    sdat.import_file(args.file_name, import_file_type)

    if args.x_max is None:
        if sdat.cell[0] is not None:
            x_max = sdat.cell[0]
        else:
            x_max = sdat.atoms['x'].max()
    else:
        x_max = args.x_max
    
    if args.y_max is None:
        if sdat.cell[1] is not None:
            y_max = sdat.cell[1]
        else:
            y_max = sdat.atoms['y'].max()
    else:
        y_max = args.y_max
    
    if args.z_max is None:
        if sdat.cell[2] is not None:
            z_max = sdat.cell[2]
        else:
            z_max = sdat.atoms['z'].max()
    else:
        z_max = args.z_max
        
    if args.x_min is None:
        x_min = 0
    else:
        x_min = args.x_min
    
    if args.y_min is None:
        y_min = 0
    else:
        y_min = args.y_min
    
    if args.z_min is None:
        z_min = 0
    else:
        z_min = args.z_min

    density = sdat.density(x_min=x_min, y_min=y_min, z_min=z_min,x_max=x_max, y_max=y_max, z_max=z_max)
    print(f'condition: {x_min} <= x <= {x_max}')
    print(f'condition: {y_min} <= y <= {y_max}')
    print(f'condition: {z_min} <= z <= {z_max}')
    print(f'density (g/cm^3): {density}')
