import networkx as nx

import pandas as pd
from models.nucleotide import Nucleotid
from pyvis.network import Network
from random import randrange

"""
the manager responsible to search in the ARN for subgraphs
"""


class Search_Manager:
    def __init__(self, datasource):
        self.datasource = datasource
        self.list_nucleotide = []
        self.list_motif = []
        self.subgraph = []
        self.tuples = []
        self.grps = dict()
        self.load_model()

    def get_pair_type_LW(self, types):
        """
        # this methode takes a text as input and return a list of relation types [tHs,cWW ,...etc]
        """
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
        types = str(types)
        if types == "" or types == "nan":
            return ll
        elif "," in types:
            tab = str(types).split(',')
            return tab
        else:
            ll.append(types)
            return ll

    def add_paire_group(self, item):
        """assosiat a relation type with a rank or a level  which will be translated later to a color when showing the ARN
        details:
        each relation type will be represented later by a color for example
        cWW (is the key in our dictionary self.grps) and the value is 1 wich will be interpreted as blue

        Args:
            item: a relation type

        Returns:
            nothig, but it insert values in a dictionary which is an atribute of the class (self.grps)

        """
        if(str(item) not in self.grps.keys()) and ("nan" not in str(item)):
            grp = randrange(100)
            if(grp in self.grps.values()):
                self.add_paire_group(item)
            else:
                self.grps[item] = item

    def load_model(self):
        """load data from the database RNANet.db

        Args:
        nothing, the methode uses the class attribute

        Returns:
            nothing,it stors data in the attributes : self.tuples and  self.list_nucleotide

        """
        df = pd.read_csv(self.datasource)
        for index, row in df.iterrows():
            nucleotid = Nucleotid(row['chain_id'], row['index_chain'], self.get_items(
                row['paired']), self.get_items(row['pair_type_LW']))
            if "," in str(row['pair_type_LW']):
                for t in self.get_items(row['pair_type_LW']):
                    self.add_paire_group(t)
                    # print(t)
            else:
                self.add_paire_group(row['pair_type_LW'])
            self.list_nucleotide.append(nucleotid)
            # it's better to create tuples the first time when we are creating our nucleotides = earn of N iterations
            if nucleotid.index_chain != 1:
                self.tuples.append(
                    tuple((str(nucleotid.index_chain - 1), str(nucleotid.index_chain))))
                # adding the new tup to the nucleotide
                if(len(nucleotid.paired_type) != 0):
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

    def get_neighbour(self, index):
        """
        get index of dirct attached pairs

        Args:
            index: index of nucleotide in the list attribute.

        Returns:
            Nucleotide (class model)

        """
        for item in self.list_nucleotide:
            if(index in item.paired):
                return item

    def get_nucleotid_by_chain_index(self, index):
        """return a nucleotide from the class attribute (self.list_nucleotides)
        Args:
            index: index of nucleotide in the chain ARN.

        Returns:
            Nucleotide (class model)

        """
        for item in self.list_nucleotide:
            if (str(item.index_chain) == str(index)):
                return item
        return None

    def check_if_exist(self, tup):
        """
        check if a tuple (a relation of two nucleotides) already exist

        Args:
            tup: a tuple who contains index of two nucleotides.
        Returns:
            True if exist and False if it doesn't exist

        """
        return tup in self.tuples

    def get_pairs(self, nucleotide):
        """
        get  pairs of a nucleotide as a list of Nucleotides (class model)

        Args:
            nucleotide: a nucleotide (class model).
        Returns:
            The return value. True for success, False otherwise.

        """
        pairs = nucleotide.paired
        ll = []
        for item in pairs:
            ll.append(self.get_nucleotid_by_chain_index(item))
        return ll

    def new_draw(self):
        """
        draw the ARN using the class attributes

        Args:
            nothing

        Returns:
            a local generated web page

        """
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

    def get_index_of_nucleotid(self, id, nuks):
        """
        get the index of a given nucleotide from a given list
        Args:
            id: id of nucleotide in the chain.
            nuks: a list of nucleotides.

        Returns:
            the index in the given list

        """
        # print("looking for id= " + str(id))
        for i in range(len(nuks)):
            # print("chain index x : "+str(item.index_chain))
            if str(nuks[i].index_chain) == str(id):
                return i
        return None

    def get_nucleotid_by_id(self, id, nuks):
        """
        search for a nucleotide in the list in param by the id 
        of the nucleotide in the chain ARN

        Args:
            id: id of nucleotide in the chain.
            nuks: a list of nucleotides.

        Returns:
            The return value. True for success, False otherwise.

        """
        # print("looking for id= " + str(id))
        for item in nuks:
            # print("chain index x : "+str(item.index_chain))
            if str(item.index_chain) == str(id):
                return item
        return None

    def compare_nucleotid(self, nuk_motif, nuk_chain):
        """
        compare nucleotides in terms of relationtypes (used searching for a subgraph)

        Args:
            param1: The first parameter.
            param2: The second parameter.

        Returns:
            The return value. True for success, False otherwise.

        """
        isSame = True
        for item in nuk_motif.paired_type:
            if item not in nuk_chain.paired_type:
                isSame = False
        return isSame


# algorithm subgraph area


    def is_condidat(self, node):
        """detect all possible entries to start researching for a subgraph
        Args:
            node: a nucleotide

        Returns:
            The return value. True for success, False otherwise.

        """
        list_condidat = []
        for nuk in self.list_nucleotide:
            isSame = True
            for item in node.paired_type:
                if item not in nuk.paired_type:
                    isSame = False
            if(isSame == True):
                list_condidat.append(nuk)
        return list_condidat

    def neighbors(self, node, list_nuks):
        """get all neighbors of a nucleotide
        Args:
            node: a nucleotide.
            list_nuks: list of nucleotide.

        Returns:
            a list of nucleotides

        """
        ll = list()
        for item in node.paired:
            nuk = self.get_nucleotid_by_id(item, list_nuks)
            if nuk != None:
                ll.append(nuk)
        nuk = self.get_nucleotid_by_id(node.index_chain + 1, list_nuks)
        if nuk != None:
            ll.append(nuk)
        return ll

    def nuk_exist_in_list(self, nuk, list_nuk):
        """check if a given nucleotide already exist in a given list
        Args:
            nuk: nucleotide.
            param2: list of nucleotides.

        Returns:
            True or False

        """
        for item in list_nuk:
            if nuk.chain_id == item.chain_id and nuk.index_chain == item.index_chain:
                return True
        return False

    def sub_graph(self, node_id, motif_id, stack):
        """search for a subgraph
        Args:
            motif_id: id of a nucleotide from motif.
            node_id: id of a nucleotide from ARN chain.
        Returns:
            stack: a graph of nucleotide (representation contigue).
        """
        # print("node_id : " + str(node_id) + " motif_id: " +              str(motif_id) + " = " + str(stack))
        # if(len(stack) == len(self.list_motif)):
        #    return stack
        nuk = self.get_nucleotid_by_id(node_id, self.list_nucleotide)
        nuk_motif = self.get_nucleotid_by_id(motif_id, self.list_motif)
        if self.compare_nucleotid(nuk_motif, nuk):
            if not self.nuk_exist_in_list(nuk, stack):
                stack.append(nuk)
            nuk_neighbors = self.neighbors(nuk, self.list_nucleotide)
            motif_neighbors = self.neighbors(nuk_motif, self.list_motif)
            # print("nuk_nei : " + str(nuk_neighbors))
            # print("motif_nei : " + str(motif_neighbors))
            if len(motif_neighbors) != 0:
                for item in motif_neighbors:
                    for ne in nuk_neighbors:
                        if self.compare_nucleotid(item, ne):
                            self.sub_graph(
                                ne.index_chain, item.index_chain, stack)
        return stack

    def draw_list(self):
        """
        draw the ARN using the class attributes

        Args:
            nothing

        Returns:
            a local generated web page

        """
        got_net = Network(height="750px", width="100%",
                          bgcolor="#222222", font_color="white")
        got_net.heading = "motif detected"
        # set the physics layout of the network
        got_net.barnes_hut()
        cp = 0
        for nuk in self.list_motif:
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

    def test_subgraph(self):
        """
        triger searching subgraph algorithm
        """
        self.list_motif = []
        n1 = Nucleotid(0, 1, [4], ["tHS"])
        n2 = Nucleotid(1, 2, [3], ["cWW"])
        n3 = Nucleotid(1, 3, [4], [])
        self.list_motif.append(n1)
        self.list_motif.append(n2)
        self.list_motif.append(n3)
        # self.search_for_sub(self.list_motif, 0, [], 0)
        condidat = self.is_condidat(n1)
        print(" size of motif: " + str(len(self.list_motif)))
        print("condidat: " + str(condidat))
        for c in condidat:
            res = self.sub_graph(
                c.index_chain, n1.index_chain, [])
            # if(len(res) == len(self.list_motif)):
            print(
                "========================================================================================")
            print(res)
            print(
                "========================================================================================")
        self.draw_list()

        # res = self.is_condidat(n2)
        # self.detect_all_motif()
        # print(self.subgraph)
        # print(self.list_nucleotide)
