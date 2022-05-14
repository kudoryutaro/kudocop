from .import_para import ImportPara
from .import_config import ImportConfig
from .import_input import ImportInput
from .import_dumppos import ImportDumppos
from .import_dumpbond import ImportDumpbond
from .import_xyz import ImportXyz


class ImportFile(
    ImportPara,
    ImportConfig,
    ImportInput,
    ImportDumppos,
    ImportDumpbond,
    ImportXyz
):
    def __init__():
        super().__init__()
