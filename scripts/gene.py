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
        exones = sorted(self.exons.keys())
        startexon = exones[0]
        endexon = exones[len(exones)-1]
        if self.exons[startexon][0] < self.exons[startexon][1]:
            startpos = self.exons[startexon][0]
        else:
            startpos = self.exons[startexon][1]
        if self.exons[endexon][0] < self.exons[endexon][1]:
            endpos = self.exons[endexon][1]
        else:
            endpos = self.exons[endexon][0]
        return (startpos, endpos)
