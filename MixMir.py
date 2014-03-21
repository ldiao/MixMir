# python wrapper for mixMIR
# Things to make sure you have:  gemma version 0.94 (allows
# 	eigen decomposition of the kinship matrix to be saved).
# Things which must be in the same folder as this script:
# 	- gemma version 0.94 executable
#	- sequence/utrs fasta file:  refseq IDs
#	- expression file: tab-separated file of gene names and expression values
#	- mirna fasta file

# import helper files
import parseAll
import analyze
import os,re,sys,argparse

def doAll(doKin=True):
	print('Parsing files for PLINK')	
	runParse(doKin=doKin)

	print('Finished parsing files!\nRunning PLINK...........................')
	runPlink()

	print('PLINK files generated!\nSolving MLM.................')
	if useFast == False:
		runGEMMA(firstTime=True)
	elif useFast == True:
		phenoFile=exprf+'.fastlmmc'
		runFAST(phenoFile)

	print('MLM solve complete!\nAnalyzing results.......................')
	getResults()

	print('Analysis complete!')
	cleanFiles()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Does all the parsing, including:
# - Making the PED file and the MAP file for PLINK 
def runParse(frac=False,doKin=True):
	outPedFile = plinkf + '.ped'
	outMapFile = plinkf + '.map'
	parseAll.doAll(doKin=doKin,seqf=seqf,exprf=exprf,outfnkin=kinf,outPedFile=outPedFile,outMapFile=outMapFile,k=k,frac=frac,useFast=useFast)

# Parse .ped and .map files through PLINK
# to create .bed
def runPlink():
	com = './plink --noweb --make-bed --file '+plinkf+' --no-fid --no-parents --no-sex --map3 --out '+plinkf
	os.system(com)

# Run GEMMA
# firstTime == True implies it's the first time you run GEMMA, 
# 	so you have to perform the eigen decomposition.  Otherwise we
# 	can utilize a pre-decomposed matrix.
def runGEMMA(firstTime=True):
	oFile = os.path.basename(plinkf)
	
	if firstTime == True:
		com = './gemma-0.94 -bfile '+plinkf+' -k '+kinf+' -eigen -o '+oFile+';'
		os.system(com)
	
	# first run the eigen decomposition so that it doesn't need to be run in the future
	dFile = 'output/'+oFile+'.eigenD.txt'
	uFile = 'output/'+oFile+'.eigenU.txt'
	com = './gemma-0.94 -bfile '+plinkf+' -d '+dFile+' -u '+uFile+' -lmm 4 -o '+oFile+';'
	os.system(com)

def runFAST(phenoFile):
	oFile = os.path.basename(plinkf)
	com = './fastlmmc -bfile '+plinkf+' -sim '+kinf+'.fastlmmc'+' -pheno '+phenoFile+';'
	os.system(com)

# parse results
def getResults():
	global resf

	oFile = os.path.basename(plinkf)
	if useFast == False:
		resFile = 'output/'+oFile+'.assoc.txt'
		resf = resf + '.gemma'
	elif useFast == True:
		resFile = oFile+'.out.txt'
		resf = resf + '.fastlmmc'

	analyze.doAll(resFile=resFile,mirFile=mirf,seqFile=seqf,N=Nkeep,useFast=useFast,outfn=resf)


# move fastLMM iles to proper names
def cleanFiles():
	# move fastLMMC dumped output file
	if useFast == True:
		fastlmmc_outf = os.path.basename(plinkf)+'.out.txt'
		new_fastlmmc_outf = os.path.basename(plinkf)+'.fastlmmc-out.txt'
		
		com = 'mv ' + fastlmmc_outf + ' ' + new_fastlmmc_outf 
		os.system(com)

if __name__ == "__main__":

	global seqf,exprf,k,N,mirf,useFast,out

	parser = argparse.ArgumentParser(description='Tell me which files and format to use')
	parser.add_argument('--seqf',type=str,required=True,help='UTR sequence file in fasta format')
        parser.add_argument('--exprf',type=str,required=True,help='Phenotype file:  column 1 is ID; column 2 is expression')
        parser.add_argument('--kinf',type=str,required=False,help='If a kinship file has already been computed, can designate here')
        parser.add_argument('--k',type=str,default=6,help='Motif length for kinship matrix (default is 6)')
        parser.add_argument('--N',type=str,default=20,help='How many top results to analyze (default is 20)')
        parser.add_argument('--mirf',type=str,required=True,help='miRNA sequence fasta file')
        parser.add_argument('--fast',type=str,default=0,help='Use option if using FastLMM to solve the mixed linear models (default is False)')
        parser.add_argument('--out',type=str,required=False,default='mixMIR-out',help='Results output file basename (default is mixMIR-out)')

	args = parser.parse_args()

	seqf = args.seqf
	exprf = args.exprf
	k = int(args.k)
	Nkeep = int(args.N)
	mirf = args.mirf
	out = args.out

	# designate automatic files for PED/BED files
	global kinf,plinkf,resf
	odir = os.path.dirname(out)
	obase = os.path.basename(out)

	# if a kinship matrix is already specified, we can skip a few steps and save some time
	if (args.kinf):
		kinf = args.kinf
		if kinf[-9:] == '.fastlmmc':
			kinf = kinf[:-9]
		doKin = False
	else:
		kinf = os.path.join(odir,obase+'-kin.txt')
		doKin = True

	plinkf = out
	resf = os.path.join(odir,obase+'-mixMIR-results.txt')

	# if running fastLMM, file parsing is slightly different
        useFast = args.fast
	if useFast in ['False','F','false','f','0']:
		useFast = False
	elif useFast in ['True','T','true','t','1']:
		useFast = True

	if useFast == False:
		print('Solving MLM with GEMMA')
	elif useFast == True:
		print('Solving MLM with FastLMM')

	doAll(doKin=doKin)

