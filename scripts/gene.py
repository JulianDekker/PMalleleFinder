class Gene:
    def __init__(self, name, exons=None):
        self.name = name
        if exons is None:
            self.exons = {}
        else:
            self.exons = exons

    def add_exon(self, number, position):
        self.exons[str(number)] = position

    def get_exons(self):
        return self.exons

    def get_name(self):
        return self.name

