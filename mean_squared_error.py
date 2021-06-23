import os
import glob
import pathlib
import pandas as pd
import numpy as np 
from sklearn.metrics import mean_squared_error
# for standard error 
from scipy import stats
from scipy.stats import sem
import researchpy as rp

path = str(pathlib.Path(__file__).parent.absolute())

# import data 
conc_average = pd.read_csv(path + '/read_out_precision_average_conc_p.csv', sep=',')
conc_random = pd.read_csv(path + '/read_out_precision_random_conc_p.csv', sep=',')
conc_cilium = pd.read_csv(path + '/read_out_precision_cilium_conc_p.csv', sep=',')

print(conc_average['conc_average'])

# caluclating MSE
rmse_average_cilium = mean_squared_error(conc_average['conc_CV_average'], conc_cilium['conc_CV_cilium'], squared=False)
rmse_average_random = mean_squared_error(conc_average['conc_CV_average'], conc_random['conc_CV_random'], squared=False)
rmse_random_cilium = mean_squared_error(conc_cilium['conc_CV_cilium'], conc_random['conc_CV_random'], squared=False)

print(rmse_average_cilium)
print(rmse_average_random)
print(rmse_random_cilium)

# calclating max distance 
max_difference_cilium_average = max(abs(conc_average['conc_CV_average'] - conc_cilium['conc_CV_cilium']))
max_difference_random_average = max(abs(conc_average['conc_CV_average'] - conc_random['conc_CV_random']))
max_difference_random_cilium =  max(abs(conc_cilium['conc_CV_cilium'] - conc_random['conc_CV_random']))

print(max_difference_cilium_average)
print(max_difference_random_average)
print(max_difference_random_cilium)


# count how often the elements in one method are higher 
print(len(conc_cilium['conc_CV_cilium']))

print(len([i for i, j in zip(conc_average['conc_CV_average'], conc_cilium['conc_CV_cilium']) if i > j]))
print(len([i for i, j in zip(conc_random['conc_CV_random'], conc_cilium['conc_CV_cilium']) if i > j]))
print(len([i for i, j in zip(conc_random['conc_CV_random'], conc_average['conc_CV_average']) if i > j]))

# Calculate average difference between the two 

diff_av_cil = conc_cilium['conc_CV_cilium'] - conc_average['conc_CV_average'] 
print('mean difference average/cilium')
print(np.mean(diff_av_cil))
print('SE av - cil')
print(sem(diff_av_cil))
print(np.std(diff_av_cil)/np.sqrt(len(conc_cilium)))

diff_av_rand = conc_random['conc_CV_random'] - conc_average['conc_CV_average'] 

print(diff_av_rand)

print('mean difference average/random')
print(np.mean(diff_av_rand))
print('SE av - rand')
print(sem(diff_av_rand))
print(np.std(diff_av_rand)/np.sqrt(len(conc_average)))

diff_cilium_rand = conc_random['conc_CV_random'] - conc_cilium['conc_CV_cilium'] 
print('mean difference cilium/random')
print(np.mean(diff_cilium_rand))
print('SE cil - rand')
print(sem(diff_cilium_rand))
print(np.std(diff_cilium_rand)/np.sqrt(len(conc_random)))

print('t test average - random')
print(stats.ttest_1samp(diff_av_rand, 0))

print('t test cilium - random')
print(stats.ttest_1samp(diff_cilium_rand, 0))

print('t test cilium - average')
print(stats.ttest_1samp(diff_av_cil , 0))