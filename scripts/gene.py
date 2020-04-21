class Gene:
    def __init__(self, name, exons=None, chro=None):
        self.name = name
        if exons is None:
            self.exons = {}
        else:
            self.exons = exons
        self.chro = chro

    def add_exon(self, number, position):
        self.exons[str(number)] = position

    def get_exons(self):
        return self.exons

    def get_name(self):
        return self.name

    def get_chro(self):
        return self.chro

    def get_gene_boundaries(self):
        minpos = 9999999999999999999999
        maxpos = 0
        strand = 0
        for _, pos in self.exons.items():
            if minpos > int(pos[0]):
                minpos = int(pos[0])
            if maxpos < int(pos[0]):
                maxpos = int(pos[0])
            if minpos > int(pos[1]):
                minpos = int(pos[1])
            if maxpos < int(pos[1]):
                maxpos = int(pos[1])
            strand = pos[2]
        return (minpos, maxpos, strand)

    def __str__(self):
        return 'Gene object. ' + self.get_name() + '. Exons: ' + str(len(self.get_exons().keys())) + '; ' + ';'.join(
            self.get_exons().keys())

    def __repr__(self):
        return 'Gene object. ' + self.get_name() + '. Exons: ' + str(len(self.get_exons().keys())) + '; ' + ';'.join(
            self.get_exons().keys())