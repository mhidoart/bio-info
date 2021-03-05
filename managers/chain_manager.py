import pandas as pd
from models.chain import Chain

"""
manager of the existing ARN chains in the database
"""
class Chain_Manager:

    def __init__(self, datasource):
        self.datasource = datasource
        self.list_chain = []
        self.load_model()
    """
    load data from the database 

    Args:
        nothing this methode use class attribute datasource

    Returns:
        nothing it stores data in the class attribute list_chain

    """
    def load_model(self):
        df = pd.read_csv(self.datasource)
        for index, row in df.iterrows():
            chain = Chain(row['chain_id'], row['chain_name'], row['eq_class'])
            self.list_chain.append(chain)

