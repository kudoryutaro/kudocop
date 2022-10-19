from .analyze_atom import AnalyzeAtom, AnalyzeAtomForSDats
from .analyze_bond import AnalyzeBond, AnalyzeBondForSDats
from .analyze_molecule import AnalyzeMolecule, AnalyzeMoleculeForSDats
from .dmol3 import DMol3KudoCop
from .packmol import Packmol
from .laich import Laich

class Analyze(
    AnalyzeAtom,
    AnalyzeBond,
    AnalyzeMolecule,
    DMol3KudoCop,
    Packmol,
    Laich,
):
    def __init__():
        pass


class AnalyzeForSDats(
    AnalyzeAtomForSDats,
    AnalyzeBondForSDats,
    AnalyzeMoleculeForSDats
):
    def __init__():
        pass
