#!/usr/bin/env python3

from abc import ABCMeta
from abc import abstractmethod

class ImportFile(object, metaclass=ABCMeta):
    @abstractmethod
    def import_file(self):
        pass
