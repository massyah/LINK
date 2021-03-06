LINK: Latent Inference of NetworKs
========================================

LINK is a tool to infer and analyze networks reprensenting signaling pathways. LINK is based on 

1. a _background interactome_ (a network of all known PPIs)
2. _a corpus_: the textual content of published documents annotating the PPIs of the interactome

To infer networks, LINK receives as input _a seed_, conisting of a starting seed of proteins and/or description of interactions between them. Interactions can be provided either as 

1. Edges from the backgroun interactome
2. PMIDs corresponding to documents annotating PPI
3. Free text 


An overview of LINK is given in the next figure :

![LINK Overview](figures/1b_algorithmic_overview_v3.png)



Two example reconstructions are provided in the [experiments](experiments) folder:

1. [Cisplatin network](experiments/131202_Daniel_Cisplatin_network/build_cisplatin_network.py)
2. [IL-2 reconstruction](experiments/131205_examples_for_BioInformatics/rebuild_IL2.py) (compared to the NetPath reference pathway)

Results for the example Il-2 reconstruction provided with the code are :

Tag      | Edges | TP Edges | REF Edges |   Edges Prec   |   Edges Recall  | Nodes | TP Nodes | REF Nodes |   Nodes Prec   |  Nodes Recall
---- |---- |---- |---- |---- |---- |---- |---- |---- |---- |----
    MST      |   9   |    7     |     96    | 0.777777777778 | 0.0729166666667 |   10  |    10    |     43    |      1.0       | 0.232558139535
Expanded 20  |   20  |    16    |     96    |      0.8       |  0.166666666667 |   15  |    15    |     43    |      1.0       | 0.348837209302
Expanded 20  |   20  |    16    |     96    |      0.8       |  0.166666666667 |   15  |    15    |     43    |      1.0       | 0.348837209302
Expanded 50  |   50  |    38    |     96    |      0.76      |  0.395833333333 |   31  |    27    |     43    | 0.870967741935 | 0.627906976744
Expanded 70  |   70  |    47    |     96    | 0.671428571429 |  0.489583333333 |   37  |    28    |     43    | 0.756756756757 | 0.651162790698
Expanded 100 |  100  |    54    |     96    |      0.54      |      0.5625     |   51  |    33    |     43    | 0.647058823529 | 0.767441860465
