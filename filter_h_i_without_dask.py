#from dask_jobqueue import SLURMCluster
#cluster= SLURMCluster(
#    project="myproj",
#    queue='regular',
#    cores=32)


#cluster.scale(jobs=2)
#client = Client(cluster)

#from dask.distributed import Client

if __name__ == '__main__':
    import sys, os
    print(sys.argv)
    filename = sys.argv[1] #int(sys.argv[1])
    
    import dask.array as da
    import dask.dataframe as dd
    import pandas as pd
    import time
    import numpy as np
    t0 = time.time()
    #filenames = os.listdir('/fast/users/keylpg_c/work/plots_statistics/results/epi2000_bih/use_data/')
    #filenames.sort()
    #filename = filenames[i]   
    data = pd.read_parquet('./results/LRP_values_au/epi2000_bih/' + filename)
    
    descriptors = pd.read_csv('./data/epi_top2000.csv').loc[:,['cell_id', 'cell_type_epi']]
    #sample_names = descriptors['cell_id'].unique()
    #print(descriptors)
    #descriptors.compute()
    described_data = data.merge(descriptors, left_on= 'sample_name',right_on='cell_id')
    #print(sample_names[i])
    #described_data_subset = described_data[described_data['sample_name'] == sample_names[i]]

    described_data_filtered = described_data.nlargest(500, columns=['LRP'])#.compute() #.reset_index()
    #print(described_data_filtered.head())
    
    try:
        os.mkdir('./results/highest_interactions/')
    except:
        pass
    
    #print(described_data_filtered['cell_id'])
    #sample_name = described_data_filtered['cell_id'].iloc[0]
    described_data_filtered.to_csv('./results/highest_interactions/high_int' +'_'+ filename +'.csv')
    print(time.time()-t0)
