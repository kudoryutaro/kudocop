from . import import_laich

def import_selector(ifn : str,ftype : str):
    if ftype == 'input':
        importer = import_laich.ImportLaich(ifn)
    
    return importer


