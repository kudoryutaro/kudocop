#!/usr/bin/env python3

import numpy as np
from . import output_abstract
from .. import simulationdat as sdat

class OutputXYZ(output_abstract.OutputFile):
    def __init__(self, ofn : str):
        super().__init__()
        self.ofn = ofn

    def output_file(self, data : sdat, out_columns, time_step):
        output_list = ['type','pos']
        output_list += out_columns
        # header_line
        header_line = []
        header_line.append(str(data.total_particle))

        ## supercell parameter
        ### for orthogonal only
        lattice_line = f"Lattice=\"{data.newcell[0]} 0.0 0.0 0.0 {data.newcell[1]} 0.0 0.0 0.0 {data.newcell[2]}\""

        ## properties
        ### S=string, I=integer, R=real, L=logical(bool)
        property_line = " Properties=species:S:1:pos:R:3"

        ## time
        time = ' Time=' + str(time_step)

        # body_line
        fmt = "{el:>3s}  {p[0]:13.9f}  {p[1]:13.9f}  {p[2]:13.9f}"
        _type = data.particles['type']
        if 'element' in output_list:
            elem = np.array([data.type_to_elem[key] for key in _type])
        else:
            elem = np.array([str(key) for key in _type])
        pos = data.particles["pos"][:,:]
        ## velocity and select?
        if 'velo' in output_list and 'sel' in output_list:
            vel = data.particles['velo'][:,:]
            property_line += ":vel:R:3:select:I:1"
            fmt += "  {v[0]:13.9f}  {v[1]:13.9f}  {v[2]:13.9f}  {sl:3d}"
            body_line = [fmt.format(el=el, p=p, v=v, sl=sl)
                        for el, p, v, sl in zip(elem, pos, vel, sl)]
        ## velocity ?
        elif 'velo' in output_list:
            vel = data.particles['velo'][:,:]
            property_line += ":vel:R:3"
            fmt += "  {v[0]:13.9f}  {v[1]:13.9f}  {v[2]:13.9f}"
            body_line = [fmt.format(el=el, p=p, v=v)
                        for el, p, v in zip(elem, pos, vel)]
        ## select ?
        elif 'sel' in output_list:
            sl = data.particles['mask']
            property_line += ":select:I:1"
            fmt += "  {sl:3d}"
            body_line = [fmt.format(el=el, p=p, sl=sl)
                        for el, p, sl in zip(elem, pos, sl)]
        ## default
        else:
            body_line = [fmt.format(el=el, p=p)
                        for el, p in zip(elem, pos)]
        
        comment_line = lattice_line + property_line + time
        header_line.append(comment_line)
        with open(self.ofn, 'w') as ofp:
            ofp.write('\n'.join(header_line) + '\n' + '\n'.join(body_line))

##############################################################################
########################            補足           ###########################
##############################################################################
"""
selectを表示する際はmaskを表示する形式にした。
selectでできることは、 ovitoの「Expression selcetion」を選択後に
Boolean expressionでselect==mask値とすることで、
mask値によるselectが可能
"""
