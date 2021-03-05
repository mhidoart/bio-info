import networkx as nx

import pandas as pd
from models.nucleotide import Nucleotid
from pyvis.network import Network
from random import randrange

# the controller who manage nucleotides
"""
model class (the controller of nucleotides)
"""


class Nucleotide_Manager:
    def __init__(self, datasource):
        self.datasource = datasource
        self.list_nucleotide = []
        self.list_motif = []
        self.subgraph = []
        self.tuples = []
        self.grps = dict()
        self.load_model()

    def get_pair_type_LW(self, types):
        # this methode takes a text as input and return a list of relation types [tHs,cWW ,...etc]
        ll = []
        for item in types:
            sublist = str(item).split(',')
            if (len(sublist) > 1):
                for item2 in sublist:
                    if (str(item2)) not in ll:
                        ll.append(str(item2))
            else:
                if (str(item)) not in ll:
                    ll.append(str(item))
        return ll

    def get_items(self, types):
        """
        this methode takes as parameter a text which contains info separated by comas ',' then the methode returns a list of info
        """
        ll = list()
        tab = str(types).split(',')

        if (len(tab) > 1):
            return tab
        return str(types)

    """assosiat a relation type with a rank or a level  which will be translated later to a color when showing the ARN
    details:
    each relation type will be represented later by a color for example
    cWW (is the key in our dictionary self.grps) and the value is 1 wich will be interpreted as blue

    Args:
        item: a relation type

    Returns:
        nothig, but it insert values in a dictionary which is an atribute of the class (self.grps)

    """

    def add_paire_group(self, item):
        print(item)
        if(str(item) not in self.grps.keys()) and ("nan" not in str(item)):
            grp = randrange(100)
            if(grp in self.grps.values()):
                self.add_paire_group(item)
            else:
                self.grps[item] = item
    """load data from the database RNANet.db

    Args:
       nothing, the methode uses the class attribute

    Returns:
        nothing,it stors data in the attributes : self.tuples and  self.list_nucleotide

    """

    def load_model(self):
        df = pd.read_csv(self.datasource)
        for index, row in df.iterrows():
            nucleotid = Nucleotid(row['chain_id'], row['index_chain'], self.get_items(
                row['paired']), self.get_items(row['pair_type_LW']))
            if "," in str(row['pair_type_LW']):
                for t in self.get_items(row['pair_type_LW']):
                    self.add_paire_group(t)
                    print(t)
            else:
                self.add_paire_group(row['pair_type_LW'])
            self.list_nucleotide.append(nucleotid)
            # it's better to create tuples the first time when we are creating our nucleotides = earn of N iterations
            if nucleotid.index_chain != 1:
                self.tuples.append(
                    tuple((str(nucleotid.index_chain - 1), str(nucleotid.index_chain))))
                # adding the new tup to the nucleotide
                nucleotid.add_tup_realtion(
                    tuple((str(nucleotid.index_chain - 1), str(nucleotid.index_chain))), nucleotid.paired_type[0])
            pairs = self.get_pairs(nucleotid)
            cp = 0
            for pair in pairs:
                if (None != pair and not self.check_if_exist(pair.index_chain)):
                    self.tuples.append(
                        tuple((str(nucleotid.index_chain), str(pair.index_chain))))

                    nucleotid.add_tup_realtion(tuple((str(nucleotid.index_chain), str(
                        pair.index_chain))), nucleotid.paired_type[cp])
                    cp += 1

    """return a nucleotide from the class attribute (self.list_nucleotides) 

    Args:
        index: index of nucleotide in the chain ARN.

    Returns:
        Nucleotide (class model)

    """

    def get_nucleotid_by_chain_index(self, index):
        for item in self.list_nucleotide:
            if (str(item.index_chain) == str(index)):
                return item
        return None
    """
    check if a tuple (a relation of two nucleotides) already exist

    Args:
        tup: a tuple who contains index of two nucleotides.
    Returns:
        True if exist and False if it doesn't exist

    """

    def check_if_exist(self, tup):
        return tup in self.tuples
    """
    get  pairs of a nucleotide as a list of Nucleotides (class model) 

    Args:
        nucleotide: a nucleotide (class model).
    Returns:
        The return value. True for success, False otherwise.

    """

    def get_pairs(self, nucleotide):
        pairs = nucleotide.paired
        ll = []
        for item in pairs:
            ll.append(self.get_nucleotid_by_chain_index(item))
        return ll

    """
    draw the ARN using the class attributes

    Args:
        nothing

    Returns:
        a local generated web page

    """

    def new_draw(self):
        got_net = Network(height="750px", width="100%",
                          bgcolor="#222222", font_color="white")
        got_net.heading = "Nucleotids"
        # set the physics layout of the network
        got_net.barnes_hut()
        cp = 0
        for nuk in self.list_nucleotide:
            for tuples in nuk.tup:
                if tuples.relation in self.grps.keys():
                    relation_type = self.grps.get(tuples.relation)
                else:
                    relation_type = 0
                relation_type = str(nuk.paired_type)
                got_net.add_node(tuples.tup[0],
                                 label="("+str(tuples.tup[0])+")" + " -> " + str(tuples.tup[1]) + " : " + str(
                    relation_type), group=relation_type)
                got_net.add_node(tuples.tup[1],
                                 label="("+str(tuples.tup[1])+")" + " -> " + str(tuples.tup[0]) + " : " + str(
                    relation_type), group=relation_type)
                got_net.add_edge(tuples.tup[0], tuples.tup[1])
        cp += 1

        # got_net.add_node(src, src, title=src)
        # got_net.add_node(dst, dst, title=dst)
        # got_net.add_edge(src, dst, value=w)

        # neighbor_map = got_net.get_adj_list()

        # add neighbor data to node hover data
        '''for node in got_net.nodes:
            node["title"] += " Neighbors:<br>" + \
                "<br>".join(neighbor_map[node["id"]])
            node["value"] = len(neighbor_map[node["id"]])'''

        got_net.show("networkRNA.html")

    """
    search for a nucleotide in the list in param by the id 
    of the nucleotide in the chain ARN

    Args:
        id: id of nucleotide in the chain.
        nuks: a list of nucleotides.

    Returns:
        The return value. True for success, False otherwise.

    """

    def get_nucleotid_by_id(self, id, nuks):
        for item in nuks:
            if item.chain_id == id:
                return item
        return None
