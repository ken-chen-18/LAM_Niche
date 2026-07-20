import scanpy as sc
import pandas as pd
import sys

data_ids = ["LAM3_Multiome", "LAM50_Multiome"]
for data_id in data_ids:
    sc_obj = sc.read_10x_h5("../input/data/" + data_id + "_filtered_feature_bc_matrix.h5")
    with open("../txt/" + data_id + '_scrub_info.txt', 'w') as sys.stdout:
         sc.external.pp.scrublet(sc_obj)
    sc_obj.obs[['doublet_score', 'predicted_doublet']].to_csv(data_id + '_doublet_scores_predictions.csv')
    sc.external.pl.scrublet_score_distribution(sc_obj, save="../txt/" + '_doublets_' + data_id + '.png')
