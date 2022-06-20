from .import_para import ImportPara
from .import_config import ImportConfig
from .import_input import ImportInput
from .import_dumppos import ImportDumppos
from .import_dumpbond import ImportDumpbond
from .import_outmol import ImportOutmol
from .import_xyz import ImportXyz
from .import_energy import ImportEnergy
from .import_dumpbond_cg import ImportDumpbondCG
from .import_out_cg import ImportOutCg


class ImportFile(
    ImportPara,
    ImportConfig,
    ImportInput,
    ImportDumppos,
    ImportDumpbond,
    ImportXyz,
    ImportOutmol,
    ImportEnergy,
    ImportDumpbondCG,
    ImportOutCg
):
    def __init__():
        super().__init__()

    def import_file(self, import_file_name: str, import_file_type: str):
        if import_file_type == 'input':
            self.import_input(import_file_name)

        elif import_file_type == 'dumpbond':
            self.import_dumpbond(import_file_name)

        elif import_file_type == 'dumppos':
            self.import_dumppos(import_file_name)

        elif import_file_type == 'para':
            self.import_para(import_file_name)

        elif import_file_type == 'config':
            self.import_config(import_file_name)

        elif import_file_type == 'xyz':
            self.import_xyz(import_file_name)

        elif import_file_type == 'outmol':
            self.import_outmol(import_file_name)

        elif import_file_type == 'energy':
            self.import_energy(import_file_name)

        elif import_file_type == 'dumpbond_cg':
            self.import_dumpbond_cg(import_file_name)

        elif import_file_type == 'out_cg':
            self.import_out_cg(import_file_name)
