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
        description='ファイルを変換する')
    parser.add_argument('-f', '--file_name', default=None, type=str,
                        help='変換前のファイル名')
    parser.add_argument('-o', '--output_file_name', default=None, type=str,
                        help='変換後のファイル名')
    parser.add_argument('-p', '--para_file_name', default='para.rd', type=str,
                        help='para.rdのファイル名')

    args = parser.parse_args()

    sdat = SimulationDat()
    sdat.import_para(args.para_file_name)
    import_file_type = analyze_file_type(args.file_name)
    sdat.import_file(args.file_name, import_file_type)

    output_file_type = analyze_file_type(args.output_file_name)
    sdat.export_file(args.output_file_name, output_file_type)