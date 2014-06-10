# analyze results of the mixed linear model solver, compare with miReduce and also the microRNAs
import re
import numpy as np
import scipy.stats as stat

import parseAll as pa


# This script takes as input either a GEMMA or a FaST-LMM output file, a miRNA fasta file and P-value cutoff. 
# It outputs the results of motifs matched to miRNAs to produce final output of MixMir.
# The columns in the output file are:
#	- Rank of motif
# 	- Motif
#	- Nominal P-value (i.e. not corrected for multiple testing, rounded to 8 decimal places) 
#	- Fixed effect coefficient (rounded to 8 decimal places)
#	- Number of UTRs containing motif
# 	- miRNAs matched, if any
def doAll(resFile='output/test.assoc.txt',mirFile='testdat/testmirs.fa',
		seqFile='testdat/test-utrs.fa',N = 20,useFast=False,outfn='MixMir-results.txt'):
    '''Parses output file from GEMMA or FaST-LMM and returns the top N results matched to miRNAs given in mirFile.
    
    Keyword arguments:
    resFile -- output results file, either from GEMMA or FaST-LMM
    mirFile -- miRNA fasta file
    seqFile -- UTRs fasta file originally analyzed
    N -- Top N motifs and matched miRNAs to be returned
    useFast -- if analysis was performed using GEMMA, must be set to False; if analysis
                was performed using FaST-LMM, must be set to True
    outfn -- filename for MixMir results
    '''
    loadMirs(mirFile=mirFile,userev2=True)

    if useFast == False:	
        res = loadGemmaRes(f=resFile,reverse=True)
    elif useFast == True:
        res = loadFastRes(f=resFile,reverse=True)
    
    matches = matchMotifList([row[0] for row in res[1:N+1]])
    
    # compute the number of UTRs that contain each motif
    Nutrs = utrCounts([row[0] for row in res[1:N+1]],seqf=seqFile,reverse=True)
    out = expand(res[1:N+1], Nutrs, matches=matches)

    print('Writing MixMir results to '+str(outfn))
    writeDat(out,outfn)

def writeDat(dat,outfn,sep='\t'):
    ofile = open(outfn,'w')
    for row in dat:
        ofile.write(sep.join([str(x) for x in row])+'\n')
    ofile.close()

def read(f,sep='\t'):
    ifile = open(f,'r')
    dat = [row.strip().split(sep) for row in ifile.readlines()]
    ifile.close()
    return(dat)

# match motifs from the mixed linear model solver to the microRNA sequences
# columns in output:
#	- Rank
# 	- Motif
#	- P-value (rounded to 8 decimal places)
#	- Bonferroni-corrected P-value (rounded to 8 decimal places)
#	- Fixed effect coefficient (rounded to 8 decimal places)
#	- Number of UTRs containing motif
# 	- miRNAs matched
def expand(res,Nutrs,matches):
    '''Parses the top N results to a nice readable format, including motif rank, miRNA matched, and P-value'''
    out = [['Rank','Motif','P-value','P-value (Bonf)','Coef','NUTRs','miRNAs Matched']]

    motifs = [row[0] for row in res]
    matched = {row[1]:[row[0],row[2]] for row in matches[1:]} # matched motifs
    pvals = {row[0]:str(round(row[1],8)) for row in res} # p-values
    pvals_bonf = {row[0]:str(round(len(res)*row[1],8)) for row in res} # Bonferroni-corrected p-values
    coefs = {row[0]:str(round(row[2],8)) for row in res} # fixed effect coefficients
    Nutrs = {row[0]:str(row[1]) for row in Nutrs} # number of UTRs

    for m in motifs:
        if m not in matched:
            out.append([str(motifs.index(m)+1),m,pvals[m],pvals_bonf[m],coefs[m],Nutrs[m],''])
        else:
            out.append([matched[m][0],m,pvals[m],pvals_bonf[m],coefs[m],Nutrs[m],matched[m][1]])
                    
    return(out)


# Loads miRNA data from .fa file
# Creates 4 global dictionaries, accounting for offset
# 	6mers matches to miRNAs, including A1 matching
# userev == True:  
def loadMirs(mirFile,userev2 = False):
    '''Load miRNA seed sequences for exact and offset matching from mirFile'''
    global mirs1, mirs2, mirs3, mirsA1
    
    f = open(mirFile)
    dat = [row.strip() for row in f.readlines()]
    f.close()

    mirs1 = {}
    mirs2 = {}
    mirs3 = {}
    mirsA1 = {}
    
    # info rows
    k = [i for i in range(len(dat)) if '>' in dat[i]]

    for i in k:
        m = re.search('>[\w\-]+',dat[i])
        mir = m.group(0)[1:]
        seq = dat[i+1]
        seq1 = seq[:6]
        seq2 = seq[1:7]
        seq3 = seq[2:8]
        seqA1 = 'T' + seq[:5]

        if seq1 in mirs1:
            mirs1[seq1].append(mir)
        else:
            mirs1[seq1] = [mir]
        
        if seq2 in mirs2:
            mirs2[seq2].append(mir)
        else:
            mirs2[seq2] = [mir]

        if seq3 in mirs3:
            mirs3[seq3].append(mir)
        else:
            mirs3[seq3] = [mir]
        
        if seqA1 in mirsA1:
            mirsA1[seqA1].append(mir)
        else:
            mirsA1[seqA1] = [mir]

# load results from Gemma .assoc.txt file
# if reverse = True:  get reverse complement of motifs
def loadGemmaRes(f='output/test.assoc.txt',reverse=True):
    '''Reads GEMMA output'''
    gem = read(f)
    gem = gem[1:]

    if reverse == False:
        gem = [[row[1],float(row[-1]),float(row[4])] for row in gem]
    else:
        gem = [[revcomp(row[1]),float(row[-1]),float(row[4])] for row in gem]

    gem.sort(key = lambda x:x[1])
    gem.insert(0,['Motif','P-value','Coefficient'])
    return(gem)


# load results from FastLMM results file
# if reverse = True:  get reverse complement of motifs
def loadFastRes(f='test.out.txt',reverse=True):
    '''Reads FaST-LMM output'''
    fast = read(f)
    fast = fast[1:]

    if reverse == False:
        fast = [[row[0],float(row[4]),float(row[9])] for row in fast]
    else:
        fast = [[revcomp(row[0]),float(row[4]),float(row[9])] for row in fast]

    fast.sort(key = lambda x:x[1])
    fast.insert(0,['Motif','P-value','Coefficient'])
    return(fast)


####################################################################
####################################################################
####################################################################
# L is a list of motifs to match
def matchMotifList(L,matchPos=['1','2','3','A1']):
    global mirs1, mirs2, mirs3, mirsA1
    c = 0
    allMatches = []
    for m in L:
        matches = matchMotif(m,matchPos=matchPos)
        if matches != []:
            c += 1
            rank = str(L.index(m) + 1)
            allMatches.append([rank] + flattenMatches(matches))

    print(c,'of',len(L),'entries matched')
    allMatches.insert(0,['Rank','Motif','miRNAs Matched'])

    return(allMatches)

# flatten a list of lists of matches for a given motif to output 
def flattenMatches(matches):
    # motif name
    newmatches = [matches[0][1]]
    mirs = ''
    for row in matches:
        newmirs = ', '.join(['['+row[0]+']'+mir for mir in row[2:]])
        if mirs != '':
            mirs = mirs + ', ' + newmirs
        else:
            mirs = newmirs
    newmatches.append(mirs)
    return(newmatches)

# given a motif, matches to microRNA database
def matchMotif(motif,matchPos=['1','2','3','A1']):
    matches = []
    if '1' in matchPos:
        if motif in mirs1:
            newmatch = ['1',motif] + mirs1[motif]
            matches.append(newmatch)
    if '2' in matchPos:
        if motif in mirs2:
            newmatch = ['2',motif] + mirs2[motif]
            matches.append(newmatch)
    if '3' in matchPos:
        if motif in mirs3:
            newmatch = ['3',motif] + mirs3[motif]
            matches.append(newmatch)
    if 'A1' in matchPos:
        if motif in mirsA1:
            newmatch = ['A1',motif] + mirsA1[motif]
            matches.append(newmatch)
    return(matches)

####################################################################
####################################################################
####################################################################
# given a list of top N motifs, returns the number of UTRs in which
# each given motif appears
# Needs to load fasta file of UTR sequences
def utrCounts(topmotifs,seqf='testdat/test-utrs.fa',reverse=True):
    utrs = pa.loadfa(fname=seqf)
    dutrs = pa.makeFaDict(utrs)
    dutrs = {x:dutrs[x][1].upper() for x in dutrs}

    if reverse == True:
        topmotifs1 = [revcomp(x).upper() for x in topmotifs]
    else:
        topmotifs1 = [x.upper() for x in topmotifs]

    d = {}
    for i in range(len(topmotifs1)):
        c = 0
        for u in dutrs:
            if topmotifs1[i] in dutrs[u]:
                c += 1
        d[topmotifs[i]] = c

    out = [[x,d[x]] for x in topmotifs]

    return(out)

####################################################################
####################################################################
####################################################################

# gets reverse complement
def revcomp(seq):
    seq = seq[::-1].upper()
    seq = map(comp,seq)
    seq = ''.join(seq)
    return(seq)

# gets complementary nucleotide
def comp(x):
    if x == 'A':
            return('T')
    elif x == 'T':
            return('A')
    elif x == 'G':
            return('C')
    elif x == 'C':
            return('G')

# determines all possible miRNA matches
def isMatch(m):
    L = [m]
    
    m1 = is53p(m)
    if m1 not in L:
        L.append(m1)
    
    m2 = isDup(m1)
    if m2 not in L:
        L.append(m2)
    
    m3 = isClose(m2)
    if m3 not in L:
        L.append(m3)
    return(L)

# determines if m is essentially a duplicate, i.e.
# mir-xx-1 mir-xx-2 mir-xx-3 etc.
def isDup(m):
    s = re.search('\w+-.+-.+-[123]',m)
    if s:
        return(m[:-2])
    else:
        return(m)

# see if mirna ends in -5p or -3p
def is53p(m):
    if m[-3:] == '-5p':
        return(m[:-3])
    elif m[-3:] == '-3p':
        return(m[:-3])
    else:
        return(m)

# see if miRNA name is "close" to another miRNA name, i.e.
# let-7a, let-7b, let-7c, etc.
def isClose(m): 
    s = re.search('\d+[a-z]$',m)
    if s:
        return(m[:-1])
    else:
        return(m)







