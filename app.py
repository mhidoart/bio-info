from flask import Flask, render_template, redirect
from managers.chain_manager import Chain_Manager
from managers.nucleotide_manager import Nucleotide_Manager
from dao.dao import Dao

app = Flask(__name__)


@app.route('/')
def hello():
    dao = Dao('RNANet.db')
    dao.load_chains()
    chain_manager = Chain_Manager("chains.csv")
    result = chain_manager.list_chain
    return render_template('index.html', result=result)


@app.route('/chain/id/<chain_id>')
def chain_by_id(chain_id):
    dao = Dao('RNANet.db')
    dao.load_chain_by_id(chain_id)
    nucleotide_manager = Nucleotide_Manager("nucleotides.csv")
    nucleotide_manager.new_draw()
    print(nucleotide_manager.grps)
    return redirect('/', code=302)


if __name__ == '__main__':
    app.secret_key = 'A0Zr98j/3yX R~XHH!jmN]LWX/,?RT'
    app.run(debug=True)
