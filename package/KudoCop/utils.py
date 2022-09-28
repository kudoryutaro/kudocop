

def analyze_file_type(file_name:str) -> str:
    """ファイル名からそのファイルの種類を解析する
        Parameters
        ----------
        file_name : str
            解析対象のファイル名
        Returns
        -------
            file_type: str
            解析対象のファイルの種類
    """
    if '.car' in file_name:
        return 'car'

    elif 'dump.pos' in file_name:
        return 'dumppos'

    elif 'dump.bond' in file_name:
        return 'dumpbond'

    elif 'input' in file_name:
        return 'input'

    elif 'xyz' in file_name:
        return 'xyz'

    elif 'config' in file_name:
        return 'config'

    elif 'energy' in file_name:
        return 'energy'

    elif 'sumforce' in file_name:
        return 'sumforce'
    
