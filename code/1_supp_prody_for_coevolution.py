from prody import *
from pylab import *

# data and dirs

tfVector = ['Homeobox', 'HLH', 'zf-C4', 'Ets', 'Fork_head', 'HMG_box', 'zf-C2H2', 'Zn_clus', 'bZIP_1']
dat_dir = "/home/luanyizhao/project/pdi/data/pfam/pfam_full_prody/"
res_dir = "/home/luanyizhao/project/pdi/results/"

for tfFam in tfVector:
	print tfFam, '...\n'
	#input msa file
	msafile = dat_dir + tfFam + '.aln.trim.fasta'
	msa = parseMSA(msafile)
	
	# occupancy calculation for each column to see if there are any positions in the MSA file that have a log of gaps.
	#showMSAOccupancy(msa, occ='res')
	# entropy calculation
	entropy = calcShannonEntropy(msa)
	#showShannonEntropy(entropy) # plot entropy
	
	# mutual information
	mutinfo = buildMutinfoMatrix(msa)
	
	#: MI normalization
	mutinfo_norm = applyMutinfoNorm(mutinfo, entropy, norm='minent')
	#: Mi correction, MIp
	mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')
	#showMutinfoMatrix(mutinfo)
	#showMutinfoMatrix(mutinfo_norm)
	#showMutinfoMatrix(mutinfo_corr)
	#: OMES, another coevolution method.
	omes = buildOMESMatrix(msa)
	#showMutinfoMatrix(omes)
	#: DI 
	di = buildDirectInfoMatrix(msa) # return errors once submitting this command.
	#showDirectInfoMatrix(di)
	#: SCA
	sca = buildSCAMatrix(msa)
	#showSCAMatrix(sca) 
	# write the Mi and entropy arrays to local files.
	writeArray(res_dir + tfFam + "_MI.txt", mutinfo, format='%f', delimiter='\t')
	writeArray(res_dir + tfFam + "_MI_norm.txt", mutinfo_norm, format='%f', delimiter='\t')
	writeArray(res_dir + tfFam + "_MIp.txt", mutinfo_corr, format='%f', delimiter='\t')
	writeArray(res_dir + tfFam + "_OMES.txt", omes, format='%f', delimiter='\t')
	writeArray(res_dir + tfFam + "_SCA.txt", sca, format='%f', delimiter='\t')
	writeArray(res_dir + tfFam + "_DCA.txt", di, format='%f', delimiter='\t' )
print 'Done!\n'
