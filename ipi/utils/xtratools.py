class listDict(list):
    def __getitem__(self, item):
        if isinstance(item,str):
            return list(ii[item] for ii in self)
        else:
            return super(listDict, self).__getitem__(item)

    @classmethod
    def fromDict(cls, dict):
        keylist = list(dict.keys())
        nb = len(dict[keylist[0]])
        for key in keylist:
            val = dict[key]
            if len(val) != nb:
                raise ValueError("The number of beads is different for the different properties.")

        lista = cls()
        for nn in range(nb):
            diction = {}
            for key in keylist:
                diction[key] = dict[key][nn]
            lista.append(diction)
        return lista