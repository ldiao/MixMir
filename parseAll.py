# parses the UTR file
import re
import numpy as np
import os
import itertools

# This script loads the UTR file (e.g. test-utrs.fa) and the gene expression file (e.g. test-exprs.txt). It computes a correlation
# matrix of all UTRs using the kmer counts in each UTR and outputs the correlation matrix (e.g. test-kin.tsv). It also outputs the
# UTR and expression data (e.g. test.ped) in a format that can be read by the program PLINK, keeping only the set of genes 
# with both UTR and expression data.  MixMir allows the use of different mixed linear model solvers and if we use the solver fastLMM 
# (set option "useFast = True") then we need to write out a tab delimited expressions file with the following data:
# gene name, gene name, expression
# The script also outputs a map file (e.g. test.map) that together with test-kin.tsv and test.ped are required by the PLINK software.

# The following is the list of options to this script:
# frac = True means fractional counts of motifs are returned, i.e. accounts for 
# 	the length of each sequence--this is the recommended setting
# frac = False means only absolute counts are returned in the matrix.
# useFast == True means we parse the data for input to fastLMM instead of GEMMA
# fastLMM requires a different format for kinship matrices, namely that header and 
# column must contain family ID and individual ID, separated by a space
# Also, fastLMM requires that phenotype data be family ID, individual ID, value
# (tab-separated)
# if doKin == False, do not create a kinship matrix
# **Note:  While the motif length used to create the kinship matrix is given by k (option kkin),
# we may want to analyze motifs of a different length, e.g. 6 for microRNAs (option kmotif)
def doAll(doKin=True,seqf='testdat/test-utrs.fa',exprf='testdat/test-exprs.txt',
		outfnkin='testdat/test-kin.tsv',
		outPedFile='testdat/test.ped',
		outMapFile='testdat/test.map',
		kkin=6,kmotif=6,frac=False,useFast=False):

	# load sequences, expressions
	utrs = loadfa(fname=seqf)
	exprs = loadmic(fname=exprf)

	# make dictionaries of transcripts and sequences, expressions
	dutrs = makeFaDict(utrs)
	dexprs = makeExprDict(exprs)

	# get genes with overlapping expression and sequence data
	genes = overlapGenes(dutrs,dexprs)

	# load motifs for kinship matrix
	kin_motifs = loadMotifs(k=kkin)
	# get motif counts for each transcript
	kin_dcounts = getCounts(genes,dutrs,kin_motifs,frac=frac)

	# clean dcounts and motifs:  remove motifs with no count information
	kin_dcounts,kin_motifs = clean(kin_dcounts,kin_motifs)

	if doKin == True:
		makeKin(dcounts=kin_dcounts,genes=genes,outfn=outfnkin,useFast=useFast)

	del kin_motifs, kin_dcounts

	# load motifs for analysis
	motifs = loadMotifs(k=kmotif)
	# get motif counts for each transcript
	dcounts = getCounts(genes,dutrs,motifs,frac=frac)

	# clean dcounts and motifs:  remove motifs with no count information
	dcounts,motifs = clean(dcounts,motifs)

	makePed(genes=genes,dcounts=dcounts,dexprs=dexprs,outfn=outPedFile)
	makeMap(motifs,outfn=outMapFile)

	# fastlmmc requires a special phenotype file
	if useFast == True:
		f = open(exprf + '.fastlmmc','w')
		for g in genes:
			f.write(g+'\t'+g+'\t'+str(dexprs[g])+'\n')
		f.close()

def load(fname,sep='\t'):
	f = open(fname,'r')
	dat = [row.strip().split(sep) for row in f.readlines()]
	f.close()
	return(dat)

def writeDat(dat,outfn,sep='\t'):
	f = open(outfn,'w')
	for row in dat:
		f.write(sep.join([str(x) for x in row])+'\n')
	f.close()

# load fasta file
def loadfa(fname='testdat/test-utrs.fa'):
	f = open(fname)
	utrs = [row.strip() for row in f.readlines()]
	f.close()

	return(utrs)

# load microarray expression file and normalize
def loadmic(fname='testdat/test-exprs.txt'):
	f = open(fname)
	exprs = [row.strip().split('\t') for row in f.readlines()]
	f.close()
	
	exprs = [[row[0],float(row[1])] for row in exprs]

	return(exprs)

# makes the output from loadfa into a dictionary
# with refseq IDs as keys and sense, sequence as entries
def makeFaDict(utrs):
	d = {}
	name,strand = getInfo(utrs[0])
	d[name] = [strand]
	seq = ''

	for row in utrs[1:]:
		if row[0] == '>':
			d[name].append(seq)

			name,strand = getInfo(row)
			d[name] = [strand]
			seq = ''
		else:
			seq = seq + row

	d[name].append(seq)
	return(d)

def makeExprDict(exprs):
	d = {}
	for row in exprs:
		if row[0] not in d:
			d[row[0]] = row[1]
		else:
			print(row[0],'already in dictionary!  Duplicate entry warning')
			return(None)

	return(d)


# string here is the '>' line in a fasta
# file, which includes the gene name, sense, etc.
def getInfo(string):
	m = re.search('(?<=refGene_)\w+_\d+',string)
	if m:
		name = m.group(0)
	else:
		print('No refseq ID detected',string)
		return(None)
	
	m = re.search('(?<=strand=)[+-]',string)
	if m:
		strand = m.group(0)
	else:
		print('No strand info detected',string)
		return(None)

	return(name,strand)

# gets utrs and exprs data only for those genes with both
# utrs and exprs both as dictionary
def overlapGenes(dutrs,dexprs):
	both = set(dutrs.keys()) & set(dexprs.keys())
	both = list(both)

	# remove short UTRs (short defined as < 10nts)
	both = [b for b in both if len(dutrs[b][1]) >= 10]
	both.sort()
	return(both)

###############################################################
############### Getting counts information ####################
############ and making kinship matrices ######################
###############################################################
def loadMotifs(k=6):
	nts = ['A','C','T','G']
	l = itertools.product(nts,repeat=k)
	motifs = [''.join(x) for x in l]
	return(motifs)

# for each gene, gets counts of each motif in the utr
def getCounts(genes,dutrs,motifs,frac=False):
	dcounts = {}
	for g in genes:
		c = [dutrs[g][1].count(m) for m in motifs]
		dcounts[g] = c

	if frac == True:
		for x in dcounts:
			dcounts[x] = fractionate(dcounts[x])

	return(dcounts)

def writeCounts(geneNames,motifs,dcounts,outfn):
	f = open(outfn,'w')
	f.write('\t'.join(motifs)+'\n')
	for gene in geneNames:
		newrow = [gene] + [str(x) for x in dcounts[gene]]
		f.write('\t'.join(newrow) + '\n')
	f.close()

# for a list of counts, c, returns
# fractional c
def fractionate(c):
	s = float(sum(c))
	if s != 0:
		c = [round(x/s,4) for x in c]
	return(c)

#####################################################################
########### Clean dcounts dictionary #############
#####################################################################
def clean(dcounts,motifs):
	remove = []
	for i in range(len(motifs)):
		if sum([dcounts[x][i] for x in dcounts]) == 0:
			remove.append(i)

	print('Removing '+str(len(remove))+' motifs with no information')
	for x in dcounts:
		dcounts[x] = [dcounts[x][i] for i in range(len(motifs)) if i not in remove]

	motifs = [motifs[i] for i in range(len(motifs)) if i not in remove]

	return(dcounts,motifs)


#####################################################################
########### Similarity matrix creation from counts data #############
#####################################################################
def makeKin(dcounts,genes,outfn,useFast=False):
	mat = [dcounts[g] for g in genes]
	K = np.corrcoef(mat)
	np.savetxt(outfn,K,delimiter='\t',fmt='%.4f')

	# if using fastLMM, have to insert header column and row
	# we will do this using bash commands
	if useFast == True:
		rowheader = [x+' '+x for x in genes]
		colheader = ['var'] + rowheader
	
		f = open('temp-colheader.txt','w')
		f.write('\t'.join(colheader)+'\n')
		f.close()

		f = open('temp-rowheader.txt','w')
		f.write('\n'.join(rowheader))
		f.close()

		os.system('paste temp-rowheader.txt '+outfn+' > temp1.txt')
		os.system('cat temp-colheader.txt temp1.txt > '+outfn+'.fastlmmc')
		os.system('rm temp-rowheader.txt temp-colheader.txt temp1.txt '+outfn)

# make .ped files
# the .ped files we make have the following columns:
# Individual ID, Phenotype, Genotype
# Genotype is denoted by either AA or TT (0 is reserved for missing genotype)
# tags running plink will then be:
# --no-fid (no family ID
# --no-parents
# --no-sex
# need a different ped file for each of dimers, trimers, ... sixmers, etc.
def makePed(genes,dcounts,dexprs,outfn=None):
	# first row of counts is the header, with motifs
	ped = []
	for g in genes:
		newrow = [g,dexprs[g]] + map(isZero,dcounts[g])
		ped.append(newrow)

	if outfn == None:
		return(ped)
	else:
		writeDat(ped,outfn=outfn,sep='\t')

# make .map files
# each row is:  [0, geneName, 0]
# must use --map3 in plink
def makeMap(motifs,outfn=None):
	mapdat = [[0,m,0] for m in motifs]
	
	if outfn == None:
		return(mapdat)
	else:
		writeDat(mapdat,outfn=outfn,sep='\t')

# fills in genotype of ped files
def isZero(x):
	if float(x) != 0:
		return('A A')
	else:
		return('T T')



if __name__ == "__main__":
	doAll()


