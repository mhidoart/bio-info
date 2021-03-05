"""nucleotide class model 
"""


class Nucleotid:
    def __init__(self, chain_id, index_chain, paired, paired_type):
        self.chain_id = chain_id
        self.index_chain = index_chain
        self.paired = []
        self.paired_type = []
        self.tup = []
        # paired and paired_type are  list
        if type(paired) == type(list()):
            self.paired = paired
        else:
            self.paired.append(paired)
        if type(paired_type) == type(list()):
            self.paired_type = paired_type
        else:
            self.paired_type.append(paired_type)

    def add_tup_realtion(self, tup, relation):
        self.tup.append(Tuples_paired(tup, relation))

    def __repr__(self):
        return " ===> self.chain_id: " + str(self.chain_id) + " self.index_chain: " + str(self.index_chain) + " paired : " + str(self.paired) + \
            " self.paired_type: " + str(self.paired_type) + "\n"


class Tuples_paired:

    def __init__(self, tup, relation):
        self.tup = tup
        self.relation = relation
