import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import statsmodels.stats.api as sms
#usage: python plot_relmat.py [snps relation file] [replicant number] [output filename] [bounds to use - sigma, or 95%CI otherwise]


snp_ibd_file = open('{0}_interest.genome'.format(sys.argv[1]),'r')
snp_ibd = list(snp_ibd_file)
snp_ibd_file.close()


obs_ibd_full = []
snp_ibd_list_ids = []
ibd_list_id1 = []
ibd_list_id2 = []
obs_ibd_pihat = []
for i in snp_ibd:
    i = i.split('\n')[0].split(' ')
    i = [x for x in i if x != '']
    obs_ibd_full.append(i)
    snp_ibd_list_ids.append([i[1],i[3]])
    ibd_list_id1.append(i[1])
    ibd_list_id2.append(i[3])
    obs_ibd_pihat.append(i[9])
ibd_col_headers = obs_ibd_full[0]
obs_ibd_full = obs_ibd_full[1:]
snp_ibd_list_ids = snp_ibd_list_ids[1:]
obs_ibd_pihat = obs_ibd_pihat[1:]
ibd_list_id1 = ibd_list_id1[1:]
ibd_list_id2 = ibd_list_id2[1:]

ibd_dict = {'ID1':ibd_list_id1, 'ID2':ibd_list_id2,'PiHat':obs_ibd_pihat}
obs_pihat_df = pd.DataFrame.from_dict(ibd_dict)


repsnpAll_pihat_mean_df = pd.DataFrame(columns=obs_pihat_df.columns,index=obs_pihat_df.index)
repsnpAll_pihat_std_df = pd.DataFrame(columns=obs_pihat_df.columns,index=obs_pihat_df.index)
repsnpAll_pihat_std_df['ID1'] = obs_pihat_df['ID1']
repsnpAll_pihat_std_df['ID2'] = obs_pihat_df['ID2']
repsnpAll_pihat_std_df['PiHat'] = 0
repsnpAll_pihat_mean_df['ID1'] = obs_pihat_df['ID1']
repsnpAll_pihat_mean_df['ID2'] = obs_pihat_df['ID2']
repsnpAll_pihat_mean_df['PiHat'] = 0

allrefsnp_pihat_dict = {}
for x in range(len(ibd_list_id2)):
    curr_key = '{0} {1}'.format(ibd_list_id1[x],ibd_list_id2[x])
    allrefsnp_pihat_dict[curr_key] = []

#Number of replicants that we need to iterate over
reps = int(sys.argv[2])

for r in range(0,reps):
    repsnp_ibd_filename = '{0}_{1}.genome'.format(sys.argv[1],r)

    repsnp_ibd_file = open(repsnp_ibd_filename,'r')
    repsnp_ibd = list(repsnp_ibd_file)
    repsnp_ibd_file.close()

    rep_ibd_list_id1 = []
    rep_ibd_list_id2 = []
    rep_obs_ibd_pihat = []
    rep_snp_ibd_list_ids = []
    for i in repsnp_ibd:
        i = i.split('\n')[0].split(' ')
        i = [x for x in i if x != '']
        rep_snp_ibd_list_ids.append([i[1],i[3]])
        rep_ibd_list_id1.append(i[1])
        rep_ibd_list_id2.append(i[3])
        rep_obs_ibd_pihat.append(i[9])
        
    rep_snp_ibd_list_ids = rep_snp_ibd_list_ids[1:]
    rep_obs_ibd_pihat = rep_obs_ibd_pihat[1:]
    rep_ibd_list_id1 = rep_ibd_list_id1[1:]
    rep_ibd_list_id2 = rep_ibd_list_id2[1:]

    
    for x in range(len(rep_ibd_list_id2)):
        curr_key = '{0} {1}'.format(ibd_list_id1[x],ibd_list_id2[x])
        allrefsnp_pihat_dict[curr_key].append(float(rep_obs_ibd_pihat[x]))
    

for r in allrefsnp_pihat_dict.items():
    id1 = r[0].split(' ')[0]
    id2 = r[0].split(' ')[1]
    rmean = np.mean(np.array(r[1]))
    rstd = np.std(np.array(r[1]))
    repsnpAll_pihat_std_df.loc[((repsnpAll_pihat_std_df.ID1==id1) & (repsnpAll_pihat_std_df.ID2==id2)),['PiHat']] = rstd
    repsnpAll_pihat_mean_df.loc[((repsnpAll_pihat_mean_df.ID1==id1) & (repsnpAll_pihat_mean_df.ID2==id2)),['PiHat']] = rmean
repsnpAll_pihat_mean_df.columns = ['ID1','ID2','rep_PiHat']
repsnpAll_pihat_std_df.columns = ['ID1','ID2','rep_std_PiHat']
obs_pihat_df.columns = ['ID1','ID2','obs_PiHat']

#obs_pihat_df.columns = ['ID1','ID2','obs_PiHat']
obs_rep_df = pd.merge(obs_pihat_df,repsnpAll_pihat_mean_df,on=['ID1','ID2'])
obs_rep_std_df = pd.merge(obs_rep_df,repsnpAll_pihat_std_df,on=['ID1','ID2'])

out_name = sys.argv[3]

#This is for testing purposes, so I can play around with the data and don't have to rerun it to replot the data generated for 1000 replicants (for example)
#This block can be commented out without any issues
out_name = sys.argv[3]
f = open("{0}_df_fortest.txt".format(out_name),'w')
for index, row in obs_rep_std_df.iterrows():
	f.write('{0}\t'.format(index))
	f.write('{0}\t'.format(row['ID1']))
	f.write('{0}\t'.format(row['ID2']))
	f.write('{0}\t'.format(row['obs_PiHat']))
	f.write('{0}\t'.format(row['rep_PiHat']))
	f.write('{0}\n'.format(row['rep_std_PiHat']))
f.close()

order = np.argsort(obs_rep_std_df['rep_PiHat'])
t = np.array(obs_rep_std_df['rep_PiHat'])[order].astype(np.float) #replicated values
sigma = np.array(obs_rep_std_df['rep_std_PiHat'])[order].astype(np.float) #standard deviation of replicate values
X = np.array(obs_rep_std_df['obs_PiHat'])[order].astype(np.float) #observed values

#If the user passes 'sigma' as the 4th argument, plot the 3-sigma range, otherwise plot the 95% CI
#if(sys.argv[4] == 'sigma'):
#    lower_bound = t-(3*sigma)
#    upper_bound = t+(3*sigma)
#    bounds_label = '3 sigma range'
#else:
ci = st.norm.interval(.95, loc=t, scale=sigma)
lower_bound = ci[0].astype(np.float)
upper_bound = ci[1].astype(np.float)
bounds_label = '95% CI range'
lower_bound = [0 for x in lower_bound if(x<0)]

#Checking how many and what percent of the observed values fall outside of the bounds (sigma/CI)
#Makes those into text which will be printed below the graphs
vals_oor = np.where(np.logical_or(X>upper_bound, X<lower_bound))[0]
percent_oor = round(((len(vals_oor)/len(X))*100),3)
oor_text = 'number o-r outside of {2} = {0}, {1}% of total o-r'.format(len(vals_oor),percent_oor,bounds_label)
vals_above = len(np.where(X>upper_bound)[0])
vals_below = len(np.where(X<lower_bound)[0])
vabove_percent = round((100*vals_above/len(X)),3)
vbelow_percent = round((100*vals_below/len(X)),3)
oor_above_text = 'number of values above bounds = {0}, {1}% of total values'.format(vals_above,vabove_percent)
oor_below_text = 'number of values below bounds = {0}, {1}% of total values'.format(vals_below,vbelow_percent)

#Plot the observed vs replicate graph, and saves it
fig, ax = plt.subplots(1)
ax.scatter(t, X, label='rep-obs', color='blue',s=3, alpha=0.75)
ax.scatter(t, t, label='rep-rep', color='black',s=3, alpha=0.75)
ax.fill_between(t, lower_bound, upper_bound, facecolor='#D3D3D3', alpha=0.6, label=bounds_label)
ax.legend(loc='best')
ax.set_xlabel('replicate value')
ax.set_ylabel('observed value')
plt.annotate(oor_text, (0,0), (0, -40), xycoords='axes fraction', textcoords='offset points', va='top')
plt.annotate(oor_above_text, (0,0), (0, -55), xycoords='axes fraction', textcoords='offset points', va='top')
plt.annotate(oor_below_text, (0,0), (0, -70), xycoords='axes fraction', textcoords='offset points', va='top')
plt.legend(loc=2, bbox_to_anchor=(1.05, 1.0))
plt.xlim([(min(t)-0.01),(max(t)+0.01)])
ax.set_aspect('equal')
ax.grid()
fig.savefig("{0}.png".format(out_name), bbox_inches='tight')
