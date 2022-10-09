from .analyze_atom import AnalyzeAtom, AnalyzeAtomForSDats
from .analyze_bond import AnalyzeBond, AnalyzeBondForSDats
from .analyze_molecule import AnalyzeMolecule, AnalyzeMoleculeForSDats
from .dmol3 import DMol3KudoCop

class Analyze(
    AnalyzeAtom,
    AnalyzeBond,
    AnalyzeMolecule,
    DMol3KudoCop
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
