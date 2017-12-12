import random
import sys
from multiprocessing import Manager, Pool
import threading
import multiprocessing as mp
import numpy as np
import os
#usage:
# python generate_cta_snps.py [input snp file] [snp frequency file] [number replications] [epsilon]

#Read in the observed SNP list
input_snp_file = open(sys.argv[1],'r')
snps_to_check = input_snp_file.readline().split(' ')[:-1]
input_snp_file.close()
#Frequency of the observed SNPs
obssnp_freq_file = open('{0}_interest.frq'.format(sys.argv[2]),'r')
obsnp_freqs = list(obssnp_freq_file)
obssnp_freq_file.close()
#Will be a dictionary, where the key is the frequency (of an observed SNP), the value will be a list of two items:
#the first item is the number of observed SNPs with that frequency, and the second is a shared list that will contain SNPs (none of which are in the observed SNPs)
mafs_to_check_dict = {}
#IDs of SNPs from the frequency file - we record this in case one of them turns out not to have a frequency for whatever reason
obsnp_ids = []
#This object handles the shared dict/lists that are shared among the workers in the pool
#manager = Manager()

num_maf_found_dict = {}

#Loop through the list of observed SNP frequencies, and make a dictionary out of the frequencies
for line in obsnp_freqs:
	h = line.split('\n')[0].split(' ')
	h = [x for x in h if x != '']
	try:
		curr_freq = float(h[4])
		obsnp_ids.append(h[1])
		if(curr_freq not in mafs_to_check_dict.keys()):
			mafs_to_check_dict[curr_freq] = [h[1]]
			num_maf_found_dict[curr_freq] = 0
			#mafs_to_check_dict[curr_freq] = [1,manager.list()]
		else:
			mafs_to_check_dict[curr_freq].append(h[1])
			#mafs_to_check_dict[curr_freq] = [(mafs_to_check_dict[curr_freq][0]+1),manager.list()]
	except ValueError:
		print('{0} does not have a valid frequency, so it will not be included'.format(h))


rep = int(sys.argv[3]) #number of replicant snp lists to generate
eps = float(sys.argv[4]) #epsilon, the tolerance for frequencies (adding a SNP to the list if it falls within eps +/- an observed frequency)



#Adding snp	frequency to dictionary; this is called by the pool
def freq_file_reader(snp_list):
	#Make a dictionairy that will store the number of times that this worker added a snp to a particular frequency 
	curr_dict = {x:0 for x in shared_dict.keys()}
	for line in snp_list:
		a = line.split('\n')[0].split(' ')
		a = [x for x in a if x != '']
		try:
			#If the frequency isn't valid, we skip it
			curr_snp = float(a[4])
			#For each SNP, check it against all of the observed frequencies to check if it's within range
			for f in curr_dict.keys():
				try:
					curr_key = float(f)
					#Don't add the snp to this dictionary entry if this worker has already added 1000 snps to it, or if the dictionary already has 30k entries for the frequency
					if(curr_dict[curr_key] == 1000 or len(shared_dict[curr_key][1]) == 30000):
						break
					if((curr_snp >= (curr_key-eps)) and (curr_snp <= (curr_key+eps)) and (a[1] not in obsnp_ids)):
						shared_dict[curr_key][1].append(a[1])
						curr_dict[curr_key]+=1
				except ValueError:
					print('{0} is not a valid frequency (this should not occur)'.format(f))
		except ValueError:
			print('{0} does not have a valid frequency, so it will not be included'.format(a))

#This splits up the list of frequencies into a list of lists, where each list has a certain number of frequencies in it
def split_freqlist(snp_freq_list):
	num_threads = 1
	for i in range(2,int(math.sqrt(len(snp_freq_list)))):
		if(len(snp_freq_list)%i == 0):
			num_threads = i
	#If we can't split up the list at all, add this entry to the list (which will not be added to any frequency lists)
	#Then it runs this function again. This way we can at least split the list in 2, if not more.
	if(num_threads == 1):
		snp_freq_list.append('CHR SNP A1 A2 MAF NCHROBS\n')
		print('the list was not able to be split so a value was added to it and it will be rechecked for splitting')
		return split_freqlist(snp_freq_list)
	snpf_split = []
	t_range = int(len(snp_freq_list)/num_threads)
	for n in range(num_threads):
		thread_start = t_range*n
		thread_end = t_range*(n+1)
		snpf_split.append(snp_freq_list[thread_start:thread_end])
	return snpf_split


def read_random_line(f, chunk_size=16):
	with open(f, 'r') as f_handle:
		f_handle.seek(0, os.SEEK_END)
		size = f_handle.tell()
		i = random.randint(0, size)
		while True:
			i -= chunk_size
			if i < 0:
				chunk_size += i
				i = 0
			f_handle.seek(i, os.SEEK_SET)
			chunk = f_handle.read(chunk_size)
			i_newline = chunk.rfind('\n')
			if i_newline != -1:
				i += i_newline + 1
				break
			if i == 0:
				break
		f_handle.seek(i, os.SEEK_SET)
		return f_handle.readline()

#This is how we run pools
# if __name__ == '__main__':
# 	freq_file = open('{0}.frq'.format(sys.argv[2]),'r')
	#snp_freqs = list(freq_file) #This is a list of all of the snp frequencies for the entire genome - a big list
	#freq_file.close()
	#np.random.shuffle(snp_freqs) #Shuffle it so we don't preferentially pick from a particular part  of the genome
	#shared_dict = manager.dict(mafs_to_check_dict) #The shared dict will be correctly handled by the manager, so all of the workers in the pool can safely access it
	# snp_maf_split = split_freqlist(snp_freqs) #split the list of frequencies
	# if(len(snp_maf_split[0]) > 2000000): #If the list is too big when it is sent to the workers, it'll crash the program, so this check just splits it again if it is too large
	# 	print("too many snp frequencies, so splitting the list again")
	# 	for f in snp_maf_split: #For each split freq list, split it again and send just that split subset to a pool of workers
	# 		curr_snps = split_freqlist(f)
	# 		with Pool(4) as p:
	# 			p.map(freq_file_reader, curr_snps)
	# else:
	# 	with Pool(4) as p:
	# 		p.map(freq_file_reader, snp_maf_split)
	# print("made frequency dictionary - generating replicate snp files")
	# freq_dict = dict(shared_dict.items()) #Make a non-shared dictionary from the shared one, and make regular lists from the shared lists
	# for x in freq_dict.values():
	# 	x[1] = list(x[1])




#Writing the SNPs from the observed frequency file to a new file, this way if there are any without frequencies it won't crash the program when it goes to generate the relatedness matrices
obsnps_file = open('refsnps_observed.txt','w')
for snp in obsnp_ids:
	obsnps_file.write(snp+' ')
obsnps_file.close()

#making a new file for each of the observed snps. These will be populated by a number of random snps that 
for snp in obsnp_ids:
	curr_obssnp_filename = "{0}_observed.txt".format(snp)
	curr_obssnp_file = open(curr_obssnp_filename,'w')
	curr_obssnp_file.close()

#freq_file = open('{0}.frq'.format(sys.argv[2]),'r')
freq_file = '{0}.frq'.format(sys.argv[2])
frq_still_to_match = list(mafs_to_check_dict.keys())
#print(mafs_to_check_dict)
while len(frq_still_to_match)>0:
	line = read_random_line(freq_file)
	a = line.split('\n')[0].split(' ')
	a = [x for x in a if x != '']
	try:
		curr_freq = float(a[4])
		for k in frq_still_to_match:
			if((curr_freq >= (k-eps)) and (curr_freq <= (k+eps)) and (a[1] not in obsnp_ids)):
				num_maf_found_dict[k] += 1
				for osnp in mafs_to_check_dict[k]:
					curr_file = open("{0}_observed.txt".format(osnp),'a')
					curr_file.write(a[1]+' ')
					curr_file.close()
			if(num_maf_found_dict[k] == 10000):
				frq_still_to_match.remove(k)
				print('10,000 matched snps found for frequency {0}'.format(k))
	except ValueError:
		print('{0} does not have a valid frequency, so it will not be included'.format(a))


sr = random.SystemRandom()
for n in range(0,rep):
	ref_snps_file = open('refsnps{0}.txt'.format(n),'w')
	snps_in_curr_file = []
	for osnp in obsnp_ids:
		curr_osnp_matched_file = open('{0}_observed.txt'.format(osnp),'r')
		curr_osnp_matched_list = list(curr_osnp_matched_file)[0].split(' ')
		curr_osnp_matched_file.close()
		while True:
			snp_to_write = sr.choice(curr_osnp_matched_list)
			if(snp_to_write not in snps_in_curr_file):
				snps_in_curr_file.append(snp_to_write)
				break
		ref_snps_file.write(snp_to_write + ' ')
	ref_snps_file.close()
#For however many replicants we want to generate, pull out a number of random snps for a given observed frequency equal to the number of observed snps that share that frequency
# sr = random.SystemRandom()
# for n in range(0,rep):
#     ref_snps_file = open('{0}{1}{2}'.format('refsnps',n,'.txt'),'w')
#     snps_in_curr_file = []
#     for m in freq_dict.items():
#         for x in range(0,int(m[1][0])):
#             while True:
#                 snp_to_write = sr.choice(m[1][1])
#                 if(snp_to_write not in snps_in_curr_file):
#                     snps_in_curr_file.append(snp_to_write)
#                     break
#             ref_snps_file.write(snp_to_write + ' ')
#     ref_snps_file.close()