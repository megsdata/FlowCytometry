import sys
import os
import yaml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import joypy
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import flowkit as fk
import umap.umap_ as umap
import anndata as ad
import seaborn as sns
import readfcs
import pynndescent
import scanpy as sc
from importlib.metadata import version

#config = "config.yaml"
config = yaml.safe_load(open("config/config.yaml"))

def get_docs():
    print("Numpy version: ", version('numpy'))

def load_data():
    base_dir = config.paths.base_dir
    sample_path = os.path.join(base_dir, config.paths.sample_path)
    wsp_path = os.path.join(base_dir, config.paths.wsp_name)
    sample_group = config.paths.sample_group
    # using flowkit's workspace to load samples
    wsp = fk.Workspace(wsp_path, fcs_samples=sample_path)
    sample_groups = wsp.get_sample_groups()
    sample_ids = sorted(wsp.get_sample_ids(sample_group))
    wsp.analyze_samples(sample_group)
    # create a list of events
    df_processed_events = []
    for sample_id in sample_ids:
        # The get_gate_events method returns all events if not given a gate
        df = wsp.get_gate_events(sample_id)
        df_processed_events.append(df)
    df_processed_events = pd.concat(df_processed_events)
    df_processed_events.head()

def get_gate_membership():
    # DataFrame to hold gate membership data, columns are gates, rows are events
    df_gate_membership = pd.DataFrame()

    # choose a sample to get the gates from, we assume all samples have the same gate tree
    for gate_name, gate_path in wsp.get_gate_ids(sample_ids[0]):
        results = []
        for sample_id in sample_ids:
            result = wsp.get_gate_membership(
                sample_id,
                gate_name=gate_name,
                gate_path=gate_path
            )
            results.append(result)

        results = np.concatenate(results)
        df_gate_membership[':'.join(list(gate_path) + [gate_name])] = results
    return df_gate_membership

def create_anndata():
    data = df_processed_events.iloc[:, 1:].values
    adata = ad.AnnData(data, dtype=data.dtype)
    adata.obs_names = np.arange(adata.shape[0]).astype('str')
    adata.var_names = df_processed_events.columns[1:]
    adata.obs['sample_id'] = pd.Categorical(df_processed_events.sample_id)
    adata.obs['sample_id'].index
    df_gate_membership.index = adata.obs['sample_id'].index
    adata.obsm['gate_membership'] = df_gate_membership
    gate_hierarchy = {sample_id: str(wsp.get_gate_hierarchy(sample_id)) for sample_id in sample_ids}
    adata.uns['gate_hierarchy'] = gate_hierarchy
    transforms = {sample_id: str(wsp.get_transforms(sample_id)) for sample_id in sample_ids}
    adata.uns['transforms'] = transforms
    adata

def print_info(adata):
    pass


def main():
    get_docs()
    df_gate_membership = load_data()
    
    
    
if __name__ == '__main__':
    main()
    sys.exit(0)