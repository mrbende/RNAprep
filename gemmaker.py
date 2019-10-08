# This code will take a directory containing GTEX-TCGA samples and combine them into a single GEM

print('Loading modules...')
import pandas as pd
import numpy as np
import os

# This is the absolute path to the directory containing your files
dir_path = '/scratch2/mrbende/KIDNEY/'

# This is the name of the output labels file for tsne plotting
labels_file = 'kidney-labels.txt'

# This is the name of the output GEM
GEM_file = 'kidney_FPKM.txt'

lst = os.listdir(dir_path)
labels = []
kidney_GEM = pd.DataFrame()
print('Beginning to combine GEMs...')
for i in lst:
    print('Appending ', i)
    GEM = pd.read_csv(dir_path + i, sep='\t', index_col=0)
    # The downloaded GTEX-TCGA data contains both Entrez and Hugo Gene IDs. We drop the Entrez ID here
    GEM.drop('Entrez_Gene_Id', axis=1, inplace=True)
    GEM = GEM.rename_axis(None)
    # The label for your sample will be the name of the file, without '.txt'
    label = i.replace('.txt', '')
    samples = np.loadtxt(list(GEM), dtype=str)
    labels_lst = list(np.repeat(label, len(samples)))
    labels.extend(labels_lst)
    
    # Here, the next GEM in the directory is added to the complete GEM. Genes are realigned to match
    kidney_GEM = pd.merge(kidney_GEM, GEM, left_index=True, right_index=True, how='outer')
    

# Saving to files
print('Saving to files...')
labels_tofile = np.asarray(labels)
np.savetxt(labels_file, labels, delimiter='\n', fmt='%s')
df = pd.DataFrame(kidney_GEM)
df.to_csv(GEM_file, sep='\t')
print('Done!')
