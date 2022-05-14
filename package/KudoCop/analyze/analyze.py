from .analyze_atom import AnalyzeAtom, AnalyzeAtomForSDats
from .analyze_bond import AnalyzeBond, AnalyzeBondForSDats
from .analyze_molecule import AnalyzeMolecule, AnalyzeMoleculeForSDats


class Analyze(
    AnalyzeAtom,
    AnalyzeBond,
    AnalyzeMolecule
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
