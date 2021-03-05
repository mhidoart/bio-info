import sqlite3
import pandas as pd

from dao.dao import Dao
from managers.nucleotide_manager import Nucleotide_Manager
from managers.chain_manager import Chain_Manager
from managers.search_Manager import Search_Manager


dao = Dao('RNANet.db')

dao.load_chains()
chain_manager = Chain_Manager("chains.csv")

dao.load_chain_by_id(105)
nucleotide_manager = Nucleotide_Manager("nucleotides.csv")
nucleotide_manager.new_draw()
search_manager = Search_Manager("nucleotides.csv")
search_manager.test_subgraph()
