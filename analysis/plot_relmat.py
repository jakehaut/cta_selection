#from ggplot import *
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import statsmodels.stats.api as sms
#usage: python plot_relmat.py [snps relation file] [replicant number] [output filename] [bounds to use - sigma, or 95%CI otherwise]

#Read in the observed SNP relatedness matrix
snp_relmat_file = open('{0}_interest.rel'.format(sys.argv[1]),'r')
snp_relmat = list(snp_relmat_file)
snp_relmat_file.close()

#Read in the observed SNP list - the one provided by the user
snp_relmat_ids_file = open('{0}_interest.rel.id'.format(sys.argv[1]),'r')
snp_relmat_ids = list(snp_relmat_ids_file)
snp_relmat_ids_file.close()

#Add individual ids to list
snp_rel_list_ids = []
for a in snp_relmat_ids:
	a = a.split('\n')[0].split('\t')
	snp_rel_list_ids.append(a[0])

#Add Observed SNP relatedness to list
snp_rel_list = []
for a in snp_relmat:
    a = a.split('\n')[0].split('\t')
    snp_rel_list.append(a)

#Make a dictionary of all observed SNPs Individual ID:Relatedness value
snp_rel_dict = {}
for n in range(len(snp_rel_list)):
    snp_rel_dict[snp_rel_list_ids[n]] = snp_rel_list[n]

#Number of replicants that we need to iterate over
reps = int(sys.argv[2])

#Set up the various data frames - one for the observed snps, and 2 for the replicant snps (one for the standard deviation, one for the mean)
snp_rel_df = pd.DataFrame.from_dict(snp_rel_dict)
repsnpAll_rel_mean_df = pd.DataFrame(columns=snp_rel_df.columns,index=snp_rel_df.index)
repsnpAll_rel_std_df = pd.DataFrame(columns=snp_rel_df.columns,index=snp_rel_df.index)

#For each of the individuals in the population we're looking at, add the ID followed by a number of lists equal to the number of indidivuals
#This will store the relatedness values from each replicant SNP set for each pair of individuals as a list
curr_snp_list = []
for i in range(0,len(snp_rel_dict.keys())):
	curr_snp_list.append([snp_rel_list_ids[i]])
	for j in range(0,len(snp_rel_dict.keys())):
		curr_snp_list[i].append([])
#Loop through each replicant SNP relatedness matrix, create a dictionary of the values, and then add each value to the corrisponding list in the curr_snp_list list
for r in range(0,reps):
	repsnp_relmat_filename = '{0}_{1}.rel'.format(sys.argv[1],r)

	repsnp_relmat_file = open(repsnp_relmat_filename,'r')
	repsnp_relmat = list(repsnp_relmat_file)
	repsnp_relmat_file.close()

	repsnp_relmat_ids_file = open('{0}.id'.format(repsnp_relmat_filename),'r')
	repsnp_relmat_ids = list(repsnp_relmat_ids_file)
	repsnp_relmat_ids_file.close()

	repsnp_rel_list_ids = []
	for a in repsnp_relmat_ids:
		a = a.split('\n')[0].split('\t')
		repsnp_rel_list_ids.append(a[0])

	repsnp_rel_list = []
	for a in repsnp_relmat:
	    a = a.split('\n')[0].split('\t')
	    repsnp_rel_list.append(a)

	repsnp_rel_dict = {}
	for n in range(len(repsnp_rel_list)):
		repsnp_rel_dict[repsnp_rel_list_ids[n]] = [float(i) for i in repsnp_rel_list[n]]

	for i in range(0,len(snp_rel_dict.keys())):
		for j in range(0,len(snp_rel_dict.keys())):
			curr_snp_list[i][j+1].append(repsnp_rel_dict[curr_snp_list[i][0]][j])
#After we've added all the replicant values, for each individual we convert that list of lists to an array then find the mean and std for each entry and add the value to the corrisponding data frame
for sl in curr_snp_list:
	repsnpAll_rel_mean_df[sl[0]] = np.mean(np.array(sl[1:]),axis=1)
	repsnpAll_rel_std_df[sl[0]] = np.std(np.array(sl[1:]),axis=1)

#Add a column with the list of IDs, corrisponding to the individual being compared in that row
#Then melt the data frame along that compared individual - gives us a dataframe with 3 columns, the individual, the value at each row, and the compared individual for that row.
id_names = pd.Series(snp_rel_list_ids)
snp_rel_df['ind.comp'] = id_names.values
obs_df = pd.melt(snp_rel_df,id_vars=['ind.comp'])
obs_df.columns = ['ind.comp','ind','observed.value']
#Repeat for the replicant mean dataframe
id_names = pd.Series(repsnp_rel_list_ids)
repsnpAll_rel_mean_df['ind.comp'] = id_names.values
rep_df = pd.melt(repsnpAll_rel_mean_df,id_vars=['ind.comp'])
rep_df.columns = ['ind.comp','ind','replicate.value']
#Merge the replicant and observed dataframe, matching them by the individual and compared individual
df_or = pd.merge(obs_df,rep_df,on=['ind.comp','ind'])
#Melting the standard deviation dataframe
id_names = pd.Series(repsnp_rel_list_ids)
repsnpAll_rel_std_df['ind.comp'] = id_names.values
rep_std_df = pd.melt(repsnpAll_rel_std_df,id_vars=['ind.comp'])
rep_std_df.columns = ['ind.comp','ind','replicate.std.value']
#Merging the dataframes, ending up with a 5 column data frame with individuals, replicant value, observed value, and standard deviation of replicant value
df_or_std = pd.merge(df_or,rep_std_df,on=['ind.comp','ind'])

#This is for testing purposes, so I can play around with the data and don't have to rerun it to replot the data generated for 1000 replicants (for example)
#This block can be commented out without any issues
out_name = sys.argv[3]
f = open("{0}_df_fortest.txt".format(out_name),'w')
for index, row in df_or_std.iterrows():
	f.write('{0}\t'.format(index))
	f.write('{0}\t'.format(row['ind.comp']))
	f.write('{0}\t'.format(row['ind']))
	f.write('{0}\t'.format(row['observed.value']))
	f.write('{0}\t'.format(row['replicate.value']))
	f.write('{0}\n'.format(row['replicate.std.value']))
f.close()


#Plotting the data
#This is sorting based on replicate value, so that matplotlib can handle the data. It gives the order of the array, which is used to order the other data that is plotted
order = np.argsort(df_or_std['replicate.value'])
t = np.array(df_or_std['replicate.value'])[order] #replicated values
sigma = np.array(df_or_std['replicate.std.value'])[order] #standard deviation of replicate values
X = np.array(df_or_std['observed.value'])[order].astype(np.float) #observed values

#If the user passes 'sigma' as the 4th argument, plot the 3-sigma range, otherwise plot the 95% CI
if(sys.argv[4] == 'sigma'):
	lower_bound = t-(3*sigma)
	upper_bound = t+(3*sigma)
	bounds_label = '3 sigma range'
else:
	ci = st.norm.interval(.95, loc=t, scale=sigma)
	lower_bound = ci[0].astype(np.float)
	upper_bound = ci[1].astype(np.float)
	bounds_label = '95% CI range'

#Checking how many and what percent of the observed values fall outside of the bounds (sigma/CI)
#Makes those into text which will be printed below the graphs
vals_oor = np.where(np.logical_or(X>upper_bound, X<lower_bound))[0]
percent_oor = (len(vals_oor)/len(X))*100
oor_text = 'number o-r outside of {2} = {0}, {1}% of total o-r'.format(len(vals_oor),percent_oor,bounds_label)
vals_above = len(np.where(X>upper_bound)[0])
vals_below = len(np.where(X<lower_bound)[0])
oor_above_text = 'number of values above bounds = {0}, {1}% of total values'.format(vals_above,(100*vals_above/len(X)))
oor_below_text = 'number of values below bounds = {0}, {1}% of total values'.format(vals_below,(100*vals_below/len(X)))

#Plot the observed vs replicate graph, and saves it
fig, ax = plt.subplots(1)
ax.scatter(t, X, label='rep-obs', color='blue',s=3, alpha=0.75)
ax.scatter(t, t, label='rep-rep', color='black',s=3, alpha=0.75)
ax.fill_between(t, lower_bound, upper_bound, facecolor='#D3D3D3', alpha=0.6,
				label=bounds_label)
ax.legend(loc='best')
ax.set_xlabel('replicate value')
ax.set_ylabel('observed value')
plt.annotate(oor_text, (0,0), (0, -40), xycoords='axes fraction', textcoords='offset points', va='top')
plt.annotate(oor_above_text, (0,0), (0, -55), xycoords='axes fraction', textcoords='offset points', va='top')
plt.annotate(oor_below_text, (0,0), (0, -70), xycoords='axes fraction', textcoords='offset points', va='top')
ax.set_aspect('equal')
ax.grid()
fig.savefig("{0}.png".format(out_name), bbox_inches='tight')

#Plot the (observed-replicate) vs replicate graph, and saves it
fig, ax = plt.subplots(1)
ax.scatter(t, (X-t), label='obs-rep', color='blue',s=5) 
ax.fill_between(t,(lower_bound-t), (upper_bound-t), facecolor='#D3D3D3', alpha=0.5,label=bounds_label)
ax.legend(loc='best')
ax.set_xlabel('replicate value')
ax.set_ylabel('observed - replicate value')
plt.annotate(oor_text, (0,0), (0, -40), xycoords='axes fraction', textcoords='offset points', va='top')
plt.annotate(oor_above_text, (0,0), (0, -55), xycoords='axes fraction', textcoords='offset points', va='top')
plt.annotate(oor_below_text, (0,0), (0, -70), xycoords='axes fraction', textcoords='offset points', va='top')
ax.set_aspect('equal')
ax.grid()
fig.savefig("{0}_minusrep.png".format(out_name), bbox_inches='tight')