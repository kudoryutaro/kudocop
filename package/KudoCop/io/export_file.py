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

    def export_file(self, export_file_name: str, export_file_type: str):
        if export_file_type == 'input':
            self.export_input(export_file_name)

        elif export_file_type == 'dumppos':
            self.export_dumppos(export_file_name)

        elif export_file_type == 'xyz':
            self.export_xyz(export_file_name)
