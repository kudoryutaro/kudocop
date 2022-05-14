from .export_input import ExportInput
from .export_dumppos import ExportDumppos
from .export_xyz import ExportXyz


class ExportFile(
    ExportDumppos,
    ExportInput,
    ExportXyz
):
    def __init__():
        super().__init__()
