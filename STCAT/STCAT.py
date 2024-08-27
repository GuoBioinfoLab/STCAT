import pandas as pd
import os
import scanpy as sc
from . import annotate
import json
import numpy as np
import warnings
from multiprocessing import Process, Queue
from . import logger
warnings.filterwarnings("ignore")

def STCAT(adata):

    data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    # Distinguish between T cell/non T cell
    adata.var['CD3'] = adata.var_names.isin(['CD3D','CD3E','CD3G'])
    sc.pp.calculate_qc_metrics(adata, qc_vars=['CD3'], percent_top=None, log1p=False, inplace=True)
    adata_T = adata[adata.obs.total_counts_CD3 > 0]
    adata_None_T = adata[adata.obs.total_counts_CD3 == 0]
    adata_None_T.obs['Prediction'] = 'None T'

    logger.info(f'🧮 T cell number : {adata_T.X.shape[0]}')
    logger.info(f'🧮 None T cell number : {adata_None_T.X.shape[0]}')
    
    # pre - processing
    # Calculate the amount of CD4 expression in each cell
    adata_T.var['CD4'] = adata_T.var_names.isin(['CD4'])
    sc.pp.calculate_qc_metrics(adata_T, qc_vars=['CD4'], percent_top=None, log1p=False, inplace=True)
    
    # Calculate the expression of CD8 in each cell
    adata_T.var['CD8'] = adata_T.var_names.isin(['CD8A','CD8B'])
    sc.pp.calculate_qc_metrics(adata_T, qc_vars=['CD8'], percent_top=None, log1p=False, inplace=True)

    if float(adata_T.X[:1000].max()).is_integer():
        logger.info(f"💬 The input file seems a raw count matrix.")
        logger.info(f"⚙ Do normalized!")
        sc.pp.normalize_total(adata_T, target_sum=1e4)
        sc.pp.log1p(adata_T)

    else:
        logger.info(f"💬 The input file seems not a raw count matrix.")
        
        if np.abs(np.expm1(adata_T.X[0]).sum()-10000) > 1:
            logger.info(f"❗ Invalid expression matrix, expect all genes and log1p normalized expression to 10000 counts per cell. The prediction result may not be accurate")
        else:
            logger.info(f"✔ All genes and log1p normalized expression to 10000 counts per cell.")

    adata_T.X = adata_T.X.astype(np.float64)
    np.random.seed(0)
    model = os.path.join(data_path, 'T_cell.pkl')
    predictions_voting = annotate.annotate(adata_T, model = model, majority_voting = True, min_prop =0.3)
    df = predictions_voting.predicted_labels.rename(columns={'predicted_labels': 'predicted_labels_'+ str(0), 'over_clustering': 'over_clustering_' + str(0), 'majority_voting': 'majority_voting_'+ str(0)})
    adata_T.obs = pd.merge(adata_T.obs, df, how='outer', left_index=True, right_index=True)

    type_dict = {}
    for over_clustering in predictions_voting.predicted_labels['over_clustering'].unique():
        type_dict[over_clustering] = predictions_voting.predicted_labels[predictions_voting.predicted_labels['over_clustering'] == over_clustering]['majority_voting'][0]

    adata_T_UN_D = adata_T[adata_T.obs['majority_voting_0'].isin(['Double Negative','Double Positive','MAIT','NKT','Tgd']), :]
    adata_T_UN_D.obs['Prediction'] = adata_T_UN_D.obs['majority_voting_0']

    # core
    next_step = []
    for key in type_dict.keys():
        if type_dict[key] in ['Double Negative','Double Positive','MAIT','NKT','Tgd']:
            continue
        else:
            next_step.append(key)

    # Divide CD4/8
    adata_T_sub = adata_T[adata_T.obs['over_clustering_0'].isin(next_step), :]

    CD8_COUNTS = adata_T_sub.obs.total_counts_CD8
    CD4_COUNTS = adata_T_sub.obs.total_counts_CD4
    adata_T_sub.obs['CD4_OR_CD8'] = np.where(CD8_COUNTS > CD4_COUNTS, 'CD8', 'CD4')
    adata_T_sub.obs['CD4_OR_CD8'] = adata_T_sub.obs['CD4_OR_CD8'].astype('category')
    
    global adata_T_sub_CD4
    global adata_T_sub_CD8

    # extreme cases
    cd4_or_cd8_values = adata_T_sub.obs['CD4_OR_CD8']
    cd4_count = (cd4_or_cd8_values == 'CD4').sum()
    cd8_count = (cd4_or_cd8_values == 'CD8').sum()
    total_cells = len(cd4_or_cd8_values)
    cd4_percentage = cd4_count / total_cells
    cd8_percentage = cd8_count / total_cells
    if cd4_percentage > 0.95:
        adata_T_sub.obs['CD4_OR_CD8'] = 'CD4'
    elif cd8_percentage > 0.95:
        adata_T_sub.obs['CD4_OR_CD8'] = 'CD8'

    adata_T_sub_CD4 = adata_T_sub[adata_T_sub.obs['CD4_OR_CD8'] == 'CD4']
    adata_T_sub_CD8 = adata_T_sub[adata_T_sub.obs['CD4_OR_CD8'] == 'CD8']

#############################################################################################################################################################################################################
    # process of CD4
    logger.info(f'⚡ Process of CD4')
    def process_CD4(q1):
        global adata_T_sub_CD4
        if adata_T_sub_CD4.X.shape[0] > 50:
            data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
            model = os.path.join(data_path, 'cd4.pkl')
            predictions_voting = annotate.annotate(adata_T_sub_CD4, model = model, majority_voting = True, min_prop =0.3)
            df = predictions_voting.predicted_labels.rename(columns={'predicted_labels': 'predicted_labels_'+ 'CD4', 'over_clustering': 'over_clustering_' + 'CD4', 'majority_voting': 'majority_voting_'+ 'CD4'})
            adata_T_sub_CD4.obs = pd.merge(adata_T_sub_CD4.obs, df, how='outer', left_index=True, right_index=True)

            # rank_genes_groups
            cell_counts = adata_T_sub_CD4.obs['over_clustering_CD4'].value_counts()
            single_cell_groups = cell_counts[cell_counts == 1].index
            adata_T_sub_CD4 = adata_T_sub_CD4[~adata_T_sub_CD4.obs['over_clustering_CD4'].isin(single_cell_groups)]
            adata_T_sub_CD4.uns['log1p'] = {'base': None}
            sc.tl.rank_genes_groups(adata_T_sub_CD4, 'over_clustering_CD4', method='wilcoxon' , key_added = 'model_CD4',n_genes=1000)
            
            marker = os.path.join(data_path, 'marker_CD4.txt')
            with open(marker, 'r') as txtfile_high:
                marker_CD4 = json.load(txtfile_high)

            df_0 = sc.get.rank_genes_groups_df(adata_T_sub_CD4, group=None, key='model_CD4',log2fc_min=0.05,pval_cutoff=0.5)

            markers = {}
            for group in adata_T_sub_CD4.obs['over_clustering_CD4'].unique():
                sub_df = df_0[df_0["group"] == group]
                top_60 = sub_df["names"].head(60).tolist()
                markers[group] = top_60

            type_dict_CD4 = {}
            for over_clustering in predictions_voting.predicted_labels['over_clustering'].unique():
                type_dict_CD4[over_clustering] = predictions_voting.predicted_labels[predictions_voting.predicted_labels['over_clustering'] == over_clustering]['majority_voting'][0]

            for i in list(single_cell_groups):
                del type_dict_CD4[i]

            # core
            for key in type_dict_CD4.keys():
                gene_list = markers[key]

                if type_dict_CD4[key] == 'Heterogeneous':
                    max_ratio = 0
                    max_key = None
                    for combine_key in marker_CD4.keys():
                        temp_list = marker_CD4[combine_key]
                        temp_intersection = list(set(gene_list) & set(temp_list))
                        temp_ratio = len(temp_intersection) / len(temp_list)
                        if temp_ratio > max_ratio:
                            max_ratio = temp_ratio
                            max_key = combine_key
                    if max_key is not None:
                        type_dict_CD4[key] = max_key
                else:
                    combine_list = marker_CD4[type_dict_CD4[key]]
                    intersection = list(set(gene_list) & set(combine_list))
                    if len(intersection) >= len(combine_list) / 2:
                        continue
                    else:
                        max_ratio = 0
                        max_key = None
                        for combine_key in marker_CD4.keys():
                            temp_list = marker_CD4[combine_key]
                            temp_intersection = list(set(gene_list) & set(temp_list))
                            temp_ratio = len(temp_intersection) / len(temp_list)
                            if temp_ratio > max_ratio:
                                max_ratio = temp_ratio
                                max_key = combine_key
                        if max_key is not None:
                            type_dict_CD4[key] = max_key           
            for key in type_dict_CD4.keys():
                if type_dict_CD4[key] in ['CD4 Tn', 'CD4 Tn quiescence', 'CD4 Tn adhesion', 'CD4 Tn IFN-response', 'CD4 Tn regulating'] and sum(marker in markers[key] for marker in ['CCR7', 'SELL', 'TCF7', 'LEF1']) >2:
                    continue
                elif type_dict_CD4[key] in ['CD4 Tn', 'CD4 Tn quiescence', 'CD4 Tn adhesion', 'CD4 Tn IFN-response', 'CD4 Tn regulating'] and any(marker in markers[key] for marker in marker_CD4['CD4 Tcm']):
                    type_dict_CD4[key] = 'CD4 Tcm'
            adata_T_sub_CD4.obs['Prediction'] = adata_T_sub_CD4.obs['over_clustering_CD4'].map(type_dict_CD4)
        elif 0 < adata_T_sub_CD4.X.shape[0] <= 50:
            data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
            model = os.path.join(data_path, 'cd4.pkl')
            predictions_voting = annotate.annotate(adata_T_sub_CD4, model = model, majority_voting = True, min_prop =0.3)
            df = predictions_voting.predicted_labels.rename(columns={'predicted_labels': 'predicted_labels_'+ 'CD4'})
            df['predicted_labels_CD4'] = df['predicted_labels_CD4'].value_counts().idxmax()
            adata_T_sub_CD4.obs = pd.merge(adata_T_sub_CD4.obs, df, how='outer', left_index=True, right_index=True)
            adata_T_sub_CD4.obs['Prediction'] = adata_T_sub_CD4.obs['predicted_labels_CD4']
        else:
            adata_T_sub_CD4 = adata_T_sub_CD4
        q1.put(adata_T_sub_CD4)
        logger.info(f'✅ Process of CD4 done!')
##############################################################################################################################################################################################################
    # process of CD8
    logger.info(f'⚡ Process of CD8')
    def process_CD8(q2):
        global adata_T_sub_CD8    
        if adata_T_sub_CD8.X.shape[0] > 50:
            data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
            model = os.path.join(data_path, 'cd8.pkl')
            predictions_voting = annotate.annotate(adata_T_sub_CD8, model = model, majority_voting = True, min_prop =0.3)
            df = predictions_voting.predicted_labels.rename(columns={'predicted_labels': 'predicted_labels_'+ 'CD8', 'over_clustering': 'over_clustering_' + 'CD8', 'majority_voting': 'majority_voting_'+ 'CD8'})
            adata_T_sub_CD8.obs = pd.merge(adata_T_sub_CD8.obs, df, how='outer', left_index=True, right_index=True)

            # rank_genes_groups
            cell_counts = adata_T_sub_CD8.obs['over_clustering_CD8'].value_counts()
            single_cell_groups = cell_counts[cell_counts == 1].index
            adata_T_sub_CD8 = adata_T_sub_CD8[~adata_T_sub_CD8.obs['over_clustering_CD8'].isin(single_cell_groups)]
            adata_T_sub_CD8.uns['log1p'] = {'base': None}
            sc.tl.rank_genes_groups(adata_T_sub_CD8, 'over_clustering_CD8', method='wilcoxon' , key_added = 'model_CD8',n_genes=1000)
            
            marker = os.path.join(data_path, 'marker_CD8.txt')
            with open(marker, 'r') as txtfile_high:
                marker_CD8 = json.load(txtfile_high)

            df_0 = sc.get.rank_genes_groups_df(adata_T_sub_CD8, group=None, key='model_CD8',log2fc_min=0.05,pval_cutoff=0.5)

            markers = {}
            for group in adata_T_sub_CD8.obs['over_clustering_CD8'].unique():
                sub_df = df_0[df_0["group"] == group]
                top_60 = sub_df["names"].head(60).tolist()
                markers[group] = top_60

            type_dict_CD8 = {}
            for over_clustering in predictions_voting.predicted_labels['over_clustering'].unique():
                type_dict_CD8[over_clustering] = predictions_voting.predicted_labels[predictions_voting.predicted_labels['over_clustering'] == over_clustering]['majority_voting'][0]

            for i in list(single_cell_groups):
                del type_dict_CD8[i]

            # core
            for key in type_dict_CD8.keys():
                gene_list = markers[key]
                if type_dict_CD8[key] == 'Heterogeneous':
                    max_ratio = 0
                    max_key = None
                    for combine_key in marker_CD8.keys():
                        temp_list = marker_CD8[combine_key]
                        temp_intersection = list(set(gene_list) & set(temp_list))
                        temp_ratio = len(temp_intersection) / len(temp_list)
                        if temp_ratio > max_ratio:
                            max_ratio = temp_ratio
                            max_key = combine_key
                    if max_key is not None:
                        type_dict_CD8[key] = max_key
                else:
                    combine_list = marker_CD8[type_dict_CD8[key]]
                    intersection = list(set(gene_list) & set(combine_list))
                    if len(intersection) >= len(combine_list) / 2:
                        continue
                    else:
                        max_ratio = 0
                        max_key = None
                        for combine_key in marker_CD8.keys():
                            temp_list = marker_CD8[combine_key]
                            temp_intersection = list(set(gene_list) & set(temp_list))
                            temp_ratio = len(temp_intersection) / len(temp_list)
                            if temp_ratio > max_ratio:
                                max_ratio = temp_ratio
                                max_key = combine_key
                        if max_key is not None:
                            type_dict_CD8[key] = max_key
            for key in type_dict_CD8.keys():
                if type_dict_CD8[key] in ['CD8 Tn'] and sum(marker in markers[key] for marker in ['GZMA', 'GZMB', 'GZMK', 'GZMH','CREM', 'FAM177A1', 'LDHA', 'OAZ1', 'CMC1', 'CLDND1', 'SARAF', 'FTH1', 'TRAT1', 'MLLT3', 'GGA2', 'CYSTM1', 'DUSP2', 'CST7', 'DKK3', 'ITM2C', 'GPR183', 'TRMO', 'ATP6V0C'])>1:
                    type_dict_CD8[key] = 'CD8 Tcm'
                elif type_dict_CD8[key] in ['CD8 Tcm'] and 'IL7R' not in markers[key]:
                    type_dict_CD8[key] = 'CD8 Tem'
            adata_T_sub_CD8.obs['Prediction'] = adata_T_sub_CD8.obs['over_clustering_CD8'].map(type_dict_CD8)
        elif 0 < adata_T_sub_CD8.X.shape[0] <= 50:
            data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
            model = os.path.join(data_path, 'cd8.pkl')
            predictions_voting = annotate.annotate(adata_T_sub_CD8, model = model, majority_voting = True, min_prop =0.3)
            df = predictions_voting.predicted_labels.rename(columns={'predicted_labels': 'predicted_labels_'+ 'CD8'})
            df['predicted_labels_CD8'] = df['predicted_labels_CD8'].value_counts().idxmax()
            adata_T_sub_CD8.obs = pd.merge(adata_T_sub_CD8.obs, df, how='outer', left_index=True, right_index=True)
            adata_T_sub_CD8.obs['Prediction'] = adata_T_sub_CD8.obs['predicted_labels_CD8']
        else:
            adata_T_sub_CD8 = adata_T_sub_CD8
        q2.put(adata_T_sub_CD8)
        logger.info(f'✅ Process of CD8 done!')
################################################################################################################################################################## 
    def start_processes():
        q1 = Queue()
        q2 = Queue()
        p1 = Process(target=process_CD4,args=(q1,))
        p2 = Process(target=process_CD8,args=(q2,))
        p1.start()
        p2.start()
        adata_T_sub_CD4 = q1.get()
        adata_T_sub_CD8 = q2.get()
        p1.join()
        p2.join()
        return adata_T_sub_CD4, adata_T_sub_CD8

    adata_T_sub_CD4, adata_T_sub_CD8 = start_processes()        
    adata = sc.AnnData.concatenate(adata_None_T, adata_T_UN_D, adata_T_sub_CD4, adata_T_sub_CD8, index_unique = None)
    remove_obs = ['n_genes_by_counts', 'total_counts', 'total_counts_CD3', 'pct_counts_CD3', 'total_counts_CD4', 'pct_counts_CD4', 'total_counts_CD8', 'pct_counts_CD8', 'predicted_labels_0', 
                         'over_clustering_0', 'majority_voting_0', 'CD4_OR_CD8', 'predicted_labels_CD4', 'over_clustering_CD4', 'majority_voting_CD4', 'predicted_labels_CD8', 'over_clustering_CD8', 
                         'majority_voting_CD8', 'batch']
    remove_var = ['CD3', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'CD4-1', 'CD8-1', 'CD4-2', 'CD8-2', 'CD4-3', 'CD8-3']
    adata.obs = adata.obs.drop(columns=remove_obs)
    adata.var = adata.var.drop(columns=remove_var)
    logger.info(f'✅ STCAT done!')
    return adata