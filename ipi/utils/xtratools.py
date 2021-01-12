class listDict(list):
    def __getitem__(self, item):
        if isinstance(item, str):
            lista = list()
            for ii in self:
                if item in ii.keys():
                    lista.append(ii[item])
            return lista
        else:
            return super(listDict, self).__getitem__(item)

    @classmethod
    def fromDict(cls, dict):
        keylist = list(dict.keys())
        if keylist:
            nb = len(dict[keylist[0]])
            for key in keylist:
                val = dict[key]
                if len(val) != nb:
                    raise ValueError(
                        "The number of beads is different for the different properties."
                    )

        lista = cls()
        if keylist:
            for nn in range(nb):
                diction = {}
                for key in keylist:
                    diction[key] = dict[key][nn]
                lista.append(diction)
        return lista

    def listofKeys(self):
        return list(self[0])
