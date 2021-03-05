# Encodage des structures sous forme de graphe d’interactions et recherche de motifs

### Réalisé par :

  - AMMY DRISS SOUFIANE 
  - ASSABBANE MEHDI 
  - CHOUBBY IBTISSAM 

# Contexte du projet
Dans le cadre du Module de Bioinformatique, nous avons comme projet la réalisation d’une application qui permet :
      o L’extraction des données depuis une base de données SQLITE ou bien depuis des fichiers d’extension csv.
      o Exploiter ces données afin de les représenter sous forme de graph.
      o Recherche et  détection des motifs dans une chaine. 

# Langage et bibliothèques utilisés 

    -Pour la réalisation du projet on a choisi d’utiliser le langage de programmation python vu la flexibilité qu’il offre et le grand nombre des ressources existantes sur le web.
    -Pour l’extraction des données depuis la base de données SQLITE on a utilisé le package sqlite3 de python, on a aussi utilisé le package Pandas pour manipuler les fichiers d’extension csv.
    -Pour la représentation des graphes on a utilisé la bibliothèque Networkx.
    -On a aussi utilisé la bibliothèque Pyvis afin de présenter le graph dansune page web.
    -On a utilisé le Framework Flask qui est framework web qui nous a permis de rendre toute l’application web afin qu’elle soit plus simple à utiliser. 

# Instructions d’installation 
    -Copier la base de données SQL dans le répertoire du projet.
    -Accéder au projet via la ligne de commande.
    -Installer les packages utilisées à l’aide du package manager de python
## « pip »:
    ```
    pip install pandas sql3 flask networkx pyvis
    
    ```
    -Pour lancer l’application éxecuter la commande suivante :
    ```
    python –m flask run 

    ```
    ou just :
    ```
    python app.py
    
    ```
    et pour lancer l'algo de recherche 
    ```
    python main.py
    
    ```
