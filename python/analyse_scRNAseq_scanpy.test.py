import numpy as np
import pandas as pd
import scanpy as sc
import vaeda

sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

sc.read_10x_h5()












# import scanpy as sc
# from matplotlib.pyplot import rc_context
#
# sc.set_figure_params(dpi=100, color_map="viridis_r")
# sc.settings.verbosity = 0
# sc.logging.print_header()
#
# pbmc = sc.datasets.pbmc68k_reduced()
#
# # compute clusters using the leiden method and store the results with the name `clusters`
# sc.tl.leiden(
#     pbmc,
#     key_added="clusters",
#     resolution=0.5,
#     n_iterations=2,
#     flavor="igraph",
#     directed=False,
# )
#
# # create a dictionary to map cluster to annotation label
# cluster2annotation = {
#     "0": "Monocytes",
#     "1": "NK",
#     "2": "T-cell",
#     "3": "Dendritic",
#     "4": "Dendritic",
#     "5": "Plasma",
#     "6": "B-cell",
#     "7": "Dendritic",
#     "8": "Other",
# }
#
# # add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function
# pbmc.obs["cell type"] = pbmc.obs["clusters"].map(cluster2annotation).astype("category")
# # scale and store results in layer
# pbmc.layers["scaled"] = sc.pp.scale(pbmc, copy=True).X
#
# import matplotlib.pyplot as plt
#
# with rc_context({"figure.figsize": (9, 2.5)}):
#     sc.pl.rank_genes_groups_violin(pbmc, n_genes=20, jitter=False)
#
#


