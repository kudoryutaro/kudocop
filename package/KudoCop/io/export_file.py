from .export_input import ExportInput
from .export_dumppos import ExportDumppos
from .export_xyz import ExportXyz
from .export_car import ExportCar
from .export_dp_system import ExportDPSystem
from .export_lammps_input import ExportLammpsInput
from .export_vasp import ExportVasp

class ExportFile(
    ExportDumppos,
    ExportInput,
    ExportXyz,
    ExportCar,
    ExportDPSystem,
    ExportLammpsInput,
    ExportVasp,
):
    def __init__():
        super().__init__()

    def export_file(self, export_file_name: str, export_file_type: str):
        if export_file_type == 'input':
            self.export_input(export_file_name)

        elif export_file_type == 'dumppos':
            self.export_dumppos(export_file_name)

        elif export_file_type == 'xyz':
            self.export_xyz(export_file_name)

        elif export_file_type == 'car':
            self.export_car(export_file_name)
