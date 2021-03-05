import sqlite3
import pandas as pd

"""
data access layer
responsible for getting data from RNANet.db
"""


class Dao:

    def __init__(self, datasource):
        self.datasource = datasource

    def load_chain_by_id(self, chain_id):
        """
        load an ARN chain based on a given chain ID

        Args:
            chain_id: number.

        Returns:
            nothing it writes a csv file nucleotides.csv

        """
        connection = sqlite3.connect(self.datasource)
        df = pd.read_sql_query('SELECT chain_id,index_chain, paired,pair_type_LW FROM nucleotide where chain_id =?;',
                               connection, params=[str(chain_id)])
        df.to_csv("nucleotides.csv")
        connection.close()

    def get_pair_type_LW(self, chain_id):
        """get all relation types (cWW,tHs...etc) existing in a given ARN chain

        Args:
            chain_id: number

        Returns:
            nothing, it wirites a csv file pair_types.csv

        """
        connection = sqlite3.connect(self.datasource)
        df = pd.read_sql_query('SELECT Distinct pair_type_LW FROM nucleotide ',
                               connection, params=[str(chain_id)])
        df.to_csv("pair_types.csv")
        connection.close()

    def load_chains(self):
        """get all ARN (id and name) existing in thee database
        Returns:
            nothing, it wirites a csv file chains.csv

        """
        connection = sqlite3.connect(self.datasource)
        df = pd.read_sql("""SELECT chain_id, chain_name, eq_class
                                FROM chain ;""", con=connection)
        df.to_csv("chains.csv")
        connection.close()
