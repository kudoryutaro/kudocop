from . import simulationdat as sdat

from .io import io_functions
from .atomic_info import atomic_weight

class Stream():
    def __init__(self):
        self.sdat = sdat.SimulationDat()
    
    def import_file(self,ifn : str,ftype : str):
        importer = io_functions.import_selector(ifn , ftype)
        importer.import_file(self.sdat)
