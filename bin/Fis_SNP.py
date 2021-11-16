
# returns the global Fis for each species within the file
import sys
import random
from numpy import quantile
from numpy import mean 
from numpy import std 

def calc_FIS(position):
	# position = vector of alleles found among n sequenced copies (haploid individuals) at a given position
	# Hs
	n = len(position) # number of sequenced copies
	nAlleles = [ position.count(i) for i in list(set(position)) ] # number of times each segregating alleles are found in the aligned position
	Hs = 1.0
	for i in range(len(nAlleles)):
		Hs -= (nAlleles[i]/(1.0*sum(nAlleles)))**2
	Hs *= n/(n-1.0)
	
	# Ho
	Ho = 0
	for ind in range(0, n, 2):
		if position[ind]!=position[ind+1]:
			Ho += 1.0
	Ho /= (n/2.0)
	
	# Fis
	Fis = (Hs-Ho)/Hs
	return([Hs, Ho])

## get arguments
fastaFile = sys.argv[1] # name of the fasta file

## get sequences
fasta = open(fastaFile).readlines()
seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')

## arrange sequences in order to treat them
nSequences = len(seq)
sequences = dict()

for i in range(nSequences):
	locus = seqName[i].split('|')[0]
	if locus not in sequences:
		sequences[locus] = []
	sequences[locus].append(seq[i])

# compute Fis over SNPs
## total number of SNPs
nLocus = len(sequences)
nSNPs = 0
for i in sequences:
	nSNPs += len(sequences[i][0])

## Fis
cnt = 0
Fis_SNP = [0]*nSNPs
for locus in sequences:
	L = len(sequences[locus][0]) # number of SNPs at this locus
	n = len(sequences[locus]) # number of individuals
	
	for pos_tmp in range(L):
		num = 0 # Hs-Ho
		denom = 0 # Hs
		position = []
		for ind_tmp in range(n):
			position.append(sequences[locus][ind_tmp][pos_tmp])
		res = calc_FIS(position)
		num += (res[0] - res[1])
		denom += res[0]
		Fis_SNP[cnt] = num/denom
		cnt += 1

# summary statistics = quentile of genomic distribution of SNPs
quantiles = [0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 1]
res = quantile(Fis_SNP, q=quantiles)

# output
header = "Fis_avg\tFis_std\t"
header += "\t".join([ "q{i}".format(i=i) for i in quantiles ])

output = header + "\n"
output += "{avg:.5f}\t{std:.5f}\t".format(avg=mean(Fis_SNP), std=std(Fis_SNP))
output += "\t".join([ str(round(i, 5)) for i in res ])

print(output)

