#!/usr/bin/env python3

from abc import ABCMeta
from abc import abstractmethod

class OutputFile(object, metaclass=ABCMeta):
    @abstractmethod
    def output_file(self):
        pass
