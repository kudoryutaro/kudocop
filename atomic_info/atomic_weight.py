import sys

atom_symbol_to_weight = {
    "H":    1.00794,
    "He":   4.002602,
    "Li":   6.941,
    "Be":	9.012182,
    "B":    10.811,
    "C":    12.0107,
    "N":    14.0067,
    "O":    15.9994,
    "F":    18.9984032,
    "Ne":	20.1797,
    "Na":	22.98976928,
    "Mg":	24.305,
    "Al":	26.9815386,
    "Si":	28.0855,
    "P":    30.973762,
    "S":    32.065,
    "Cl":	35.453,
    "Ar":	39.948,
    "K":    39.0983,
    "Ca":	40.078,
    "Sc":	44.955912,
    "Ti":	47.867,
    "V":    50.9415,
    "Cr":	51.9961,
    "Mn":	54.938045,
    "Fe":	55.845,
    "Co":	58.933195,
    "Ni":	58.6934,
    "Cu":	63.546,
    "Zn":	65.38,
    "Ga":	69.723,
    "Ge":	72.64,
    "As":	74.9216,
    "Se":	78.96,
    "Br":	79.904,
    "Kr":	83.798,
    "Rb":	85.4678,
    "Sr":	87.62,
    "Y":    88.90585,
    "Zr":	91.224,
    "Nb":	92.90638,
    "Mo":	95.96,
    "Tc":	99,
    "Ru":	101.07,
    "Rh":	102.9055,
    "Pd":	106.42,
    "Ag":	107.8682,
    "Cd":	112.411,
    "In":	114.818,
    "Sn":	118.71,
    "Sb":	121.76,
    "Te":	127.6,
    "I":    126.90447,
    "Xe":	131.293,
    "Cs":	132.9054519,
    "Ba":	137.327,
    "La":	138.90547,
    "Ce":	140.116,
    "Pr":	140.90765,
    "Nd":	144.242,
    "Pm":	145,
    "Sm":	150.36,
    "Eu":	151.964,
    "Gd":	157.25,
    "Tb":	158.92535,
    "Dy":	162.5,
    "Ho":	164.93032,
    "Er":	167.259,
    "Tm":	168.93421,
    "Yb":	173.054,
    "Lu":	174.9668,
    "Hf":	178.49,
    "Ta":	180.94788,
    "W":    183.84,
    "Re":	186.207,
    "Os":	190.23,
    "Ir":	192.217,
    "Pt":	195.084,
    "Au":	196.966569,
    "Hg":	200.59,
    "Tl":	204.3833,
    "Pb":	207.2,
    "Bi":	208.9804,
    "Po":	210,
    "At":	210,
    "Rn":	222,
    "Fr":	223,
    "Ra":	226,
    "Ac":	227,
    "Th":	232.03806,
    "Pa":	231.03588,
    "U":    238.02891
}

def get_weight(atom_symbol : str) -> float:
    if atom_symbol not in atom_symbol_to_weight:
        print('element symbol assignment error') 
        sys.exit()
    
    return atom_symbol_to_weight[atom_symbol]

def get_atom_symbol_to_type(ifn : str) -> dict:
    atom_symbol_to_type = {}

    with open(ifn,'r') as ifp:
        lines = ifp.readlines()
    
    type_num = 1

    for line in lines:
        sline = line.split()
        if sline[0] in atom_symbol_to_weight:
            atom_symbol_to_type[sline[0]]   = type_num
            type_num += 1


if __name__ == '__main__':
    print(get_weight('Au'))