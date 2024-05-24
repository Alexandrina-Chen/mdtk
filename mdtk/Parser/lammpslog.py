# -*- encoding: utf-8 -*-
# python 3.10.12
'''
Filename         : lammpsthermoparser.py
Description      : parser for lammps thermo output file
Time             : 2024/02/20
Author           : Sijia Chen
Version          : 1.0
Email            : sijiachen@uchicago.edu
'''


# class definition
class LammpsThermoParser:
    def __init__(self, filename):
        """
        This class is used to parse the lammps log file.
        NOTICE: the parser assumes the thermo output starts with "Step" and ends with "Loop time of"!
        NOTICE: the parser stores all data's field names in lower case.

        All fields can be accessed by their field names. 
        e.g. if the instance of this class is called `log_data`, you can access the field `step` by `log_data.step`.
        Or you can access the data dictionary directly by the data attribute.
        
        Parameters
        ----------
        filename : str
            the name of the lammps log file.

        Attributes
        ----------
        filename : str
            the name of the lammps log file.
        data : dict
            the parsed data stored in a dictionary.


        """
        self.filename = filename
        self._data = self._read_log()
        # provide direct access to all fields
        for key in self._data.keys():
            setattr(self, key, self._data[key])

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value


    def _read_log(self):
        log_data = {}
        trigger = 0
        num_fields = 0
        field_names = []
        with open(self.filename, "r") as f:
            for line in f:
                if "Time step     :" in line:
                    log_data["timestep"] = float(line.split()[-1])
                if line.startswith("Step"):
                    tokens = [x.lower() for x in line.strip().split()]
                    field_names = tokens
                    for token in tokens:
                        log_data[token] = []
                    trigger += 1
                    num_fields = len(tokens)
                    continue
                if line.startswith("Loop time"):
                    break
                if trigger == 1 and len(line.split()) == num_fields: # avoid shake information lines
                    tokens = line.split()
                    for i in range(num_fields):
                        log_data[field_names[i]].append(float(tokens[i]))
        return log_data