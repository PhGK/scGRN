#from dask_jobqueue import SLURMCluster
#cluster= SLURMCluster(
#    project="myproj",
#    queue='regular',
#    cores=32)


#cluster.scale(jobs=2)
#client = Client(cluster)

from dask.distributed import Client

if __name__ == '__main__':
    client = Client(n_workers=20, threads_per_worker=1, memory_limit = '35GB')
    
    import dask.dataframe as dd
    import pandas as pd
    import time
    import numpy as np    
    t0 = time.time()
    print(t0)    
    data = dd.read_parquet('./results/LRP_values_au/epi2000_bih/LRP_*parquet')
    #print(data)
    #print(data.iloc[:,0].compute())    
    descriptors = pd.read_csv('./data/epi_top2000.csv').loc[:,['cell_id', 'cell_type_epi']]
    #print(data.compute())


    descripted_data = data.merge(descriptors, left_on= 'sample_name',right_on='cell_id')
    #mean_interactions = descripted_data.groupby(['source_gene', 'target_gene','cell_type_epi'])['LRP'].agg(['mean']).reset_index()
    #mean_interactions = descripted_data_filtered.groupby(['source_gene', 'target_gene', 'cell_type_epi']).agg({'counts':'size', 'LRP':'mean'}).reset_index()

    #mean_interactions.compute().to_csv('/fast/users/keylpg_c/work/plots_statistics/results/mean_interactions.csv')
    t1 = time.time()
    total = t1-t0
    print(total)
    descripted_data_filtered = descripted_data #[(descripted_data['inpv']>0) & (descripted_data['tinpv']>0)]
    #mean_interactions_all_types = descripted_data_filtered.groupby(['source_gene', 'target_gene'])['LRP'].agg(['mean']).reset_index()
    
    descripted_data_filtered['counts'] = 1
    mean_interactions_all_types = descripted_data_filtered.groupby(['source_gene', 'target_gene']).agg({'counts':'size', 'LRP':'mean'}).reset_index()
    
    mean_interactions_all_types.compute().to_csv('./results/mean_interactions_all_types.csv')

