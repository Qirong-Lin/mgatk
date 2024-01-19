#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
import re
import os
import pysam
import random
import numpy as np
from collections import defaultdict

bamfile = sys.argv[1]
outpre = sys.argv[2]
mito_genome = sys.argv[3]
maxBP = sys.argv[4]
base_qual = float(sys.argv[5])
sample = sys.argv[6]
fasta_file = sys.argv[7]
alignment_quality = float(sys.argv[8])
emit_base_qualities = sys.argv[9]

# bamfile= "/mnt/cache/test_largest/AAACAGCAGTCTGCTA-1.bam"
# bamfile= "/mnt/cache/test_largest/GTAAGGGTTGTACGCA-1.bam"
# outpre = "/mnt/cache/test_largest/sumstatsBP_overlap"
# mito_genome = "/mnt/cache/test_largest/c2265_atac_sec1_mgatk_out_ho_ub_jm512_latencywait60/fasta/chrM.fasta"
# maxBP = 16569
# base_qual = 20
# sample = "test_largest"
# fasta_file = "/mnt/cache/test_largest/c2265_atac_sec1_mgatk_out_ho_ub_jm512_latencywait60/fasta/chrM.fasta"
# alignment_quality = 20
# emit_base_qualities = "True"


# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")

def writeSparseMatrix2(mid, vec1, vec2):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0 or vec2[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+"\n")

def writeSparseMatrix4(mid, vec1, vec2, vec3, vec4):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0 or vec3[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+","+str(vec3[i])+","+str(vec4[i])+"\n")

def findHighQualityBases(fwd_read, rev_read, overlap_start, overlap_end):
	# Get aligned pairs for forward and reverse reads
	fwd_pairs = fwd_read.get_aligned_pairs(matches_only=True)
	rev_pairs = rev_read.get_aligned_pairs(matches_only=True)

	# Find which read has a higher base quality at each overlapping ref position
	fwd_overlap_use_idx = []
	rev_overlap_use_idx = []

	for ref_pos in range(overlap_start, overlap_end):
		fwd_read_pos = next((qp for qp, rpos in fwd_pairs if rpos == ref_pos), None)
		rev_read_pos = next((qp for qp, rpos in rev_pairs if rpos == ref_pos), None)

		if fwd_read_pos is not None and rev_read_pos is not None:
			fwd_base_qual = fwd_read.query_qualities[fwd_read_pos]
			rev_base_qual = rev_read.query_qualities[rev_read_pos]

			# randomly choose one if the base qualities are the same
			# or choose the one with higher base quality
			if (
				(fwd_base_qual == rev_base_qual
				and random.choice([True, False]))
				or fwd_base_qual > rev_base_qual
			):
				fwd_overlap_use_idx.append(ref_pos)
			else:
				rev_overlap_use_idx.append(ref_pos)
		elif fwd_read_pos is None and rev_read_pos is not None:
			rev_overlap_use_idx.append(ref_pos)
		elif fwd_read_pos is not None and rev_read_pos is None:
			fwd_overlap_use_idx.append(ref_pos)

		return fwd_overlap_use_idx, rev_overlap_use_idx

n = int(maxBP)

# initialize with a pseudo count to avoid dividing by zero
countsA_fw = [0.00000001] * n 
countsC_fw = [0.00000001] * n 
countsG_fw = [0.00000001] * n 
countsT_fw = [0.00000001] * n 

qualA_fw = [0.0] * n
qualC_fw = [0.0] * n
qualG_fw = [0.0] * n
qualT_fw = [0.0] * n

countsA_rev = [0.00000001] * n 
countsC_rev = [0.00000001] * n 
countsG_rev = [0.00000001] * n 
countsT_rev = [0.00000001] * n 

qualA_rev = [0.0] * n
qualC_rev = [0.0] * n
qualG_rev = [0.0] * n
qualT_rev = [0.0] * n

# organize reads into a dict where key is readname
bam2 = [x for x in pysam.AlignmentFile(bamfile, "rb")]
ordered_bam2 = defaultdict(list)

for read in bam2:
	ordered_bam2[read.query_name].append(read)


for read_name in ordered_bam2:
	# disregard singlets and multiplets
	if len(ordered_bam2[read_name]) != 2:
		continue
	
	# identify fwd and rev in a pair
	read0, read1 = ordered_bam2[read_name]
	if read0.is_reverse and not read1.is_reverse:
		fwd_read, rev_read = read1, read0
	elif not read0.is_reverse and read1.is_reverse:
		fwd_read, rev_read = read0, read1
	else:
		# disregard a pair if both are the same strand
		continue
	
	# gather what we need
	fwd_seq, rev_seq = fwd_read.query_sequence, rev_read.query_sequence
	fwd_quality, rev_quality = np.array(fwd_read.query_qualities), np.array(rev_read.query_qualities)
	fwd_align_qual_read, rev_align_qual_read = fwd_read.mapping_quality, rev_read.mapping_quality
	#print(type(fwd_align_qual_read))  
	#sys.exit(0)
	
	# check alignment quality
	if fwd_align_qual_read > alignment_quality and rev_align_qual_read > alignment_quality:
		overlap_start = max(fwd_read.reference_start, rev_read.reference_start)
		overlap_end = min(fwd_read.reference_end, rev_read.reference_end)	

		# if there is no overlap, use all of fwd and rev
		if overlap_start >= overlap_end:
			fwd_use_idx = np.arange(fwd_read.reference_start, fwd_read.reference_end)
			rev_use_idx = np.arange(rev_read.reference_start, rev_read.reference_end)
		# if there is an overlap, use the high quality bases in the overlap region
		else:
			fwd_overlap_use_idx, rev_overlap_use_idx = findHighQualityBases(fwd_read, rev_read, overlap_start=overlap_start, overlap_end=overlap_end)
			
			# if reverse strand is included in the forward strand, partition the pair into fwd-only, overlap, and fwd-only
			if fwd_read.reference_start < rev_read.reference_start and fwd_read.reference_end > rev_read.reference_end:
				# merge the exclusive region and use idx in overlap region
				fwd_use_idx = np.concatenate([np.arange(fwd_read.reference_start, overlap_start), fwd_overlap_use_idx, np.arange(overlap_end, fwd_read.reference_end)])
				rev_use_idx = rev_overlap_use_idx

			# if forward strand is included in the reverse strand, partition the pair into rev-only, overlap, and rev-only
			elif fwd_read.reference_start > rev_read.reference_start and fwd_read.reference_end < rev_read.reference_end:
				fwd_use_idx = fwd_overlap_use_idx
				rev_use_idx = np.concatenate([np.arange(rev_read.reference_start, overlap_start), rev_overlap_use_idx, np.arange(overlap_end, rev_read.reference_end)])

			else:
			# partition the pair into fwd-only, overlap, and rev-only
			# merge the exclusive region and use idx in overlap region
				fwd_use_idx = np.concatenate([np.arange(fwd_read.reference_start, overlap_start), fwd_overlap_use_idx])
				rev_use_idx = np.concatenate([rev_overlap_use_idx, np.arange(overlap_end, rev_read.reference_end)])

	elif fwd_align_qual_read <= alignment_quality and rev_align_qual_read <= alignment_quality:
		# use none for either
		fwd_use_idx = np.array([])
		rev_use_idx = np.array([])
	
	elif fwd_align_qual_read > alignment_quality and rev_align_qual_read <= alignment_quality:
		# use none of rev and all of fwd
		fwd_use_idx = np.arange(fwd_read.reference_start, fwd_read.reference_end)
		rev_use_idx = np.array([])
	
	elif fwd_align_qual_read <= alignment_quality and rev_align_qual_read > alignment_quality:
		# use all of rev and none of fwd
		fwd_use_idx = np.array([])
		rev_use_idx = np.arange(rev_read.reference_start, rev_read.reference_end)


	# since mgatk will not handle indels anyway, we will use refpos instead of qpos in the fwd_use_idx and rev_use_idx.

	# handle fwd region, use refpos (pair[1]) instead of qpos (pair[0])
	fwd_aligned_pairs = fwd_read.get_aligned_pairs(True)
	fwd_region = [pair for pair in fwd_aligned_pairs if pair[1] in fwd_use_idx]
	for qpos, refpos in fwd_region:
		if refpos is not None and fwd_quality[qpos] > base_qual:
			if fwd_seq[qpos] == "A":
				qualA_fw[refpos] += fwd_quality[qpos]
				countsA_fw[refpos] += 1
			elif fwd_seq[qpos] == "C":
				qualC_fw[refpos] += fwd_quality[qpos]
				countsC_fw[refpos] += 1
			elif fwd_seq[qpos] == "G":
				qualG_fw[refpos] += fwd_quality[qpos]
				countsG_fw[refpos] += 1
			elif fwd_seq[qpos] == "T":
				qualT_fw[refpos] += fwd_quality[qpos]
				countsT_fw[refpos] += 1
	
	# handle rev region, use refpos (pair[1]) instead of qpos (pair[0])
	rev_aligned_pairs = rev_read.get_aligned_pairs(True)
	rev_region = [pair for pair in rev_aligned_pairs if pair[1] in rev_use_idx]
	for qpos, refpos in rev_region:
		if refpos is not None and rev_quality[qpos] > base_qual:
			if rev_seq[qpos] == "A":
				qualA_rev[refpos] += rev_quality[qpos]
				countsA_rev[refpos] += 1
			elif rev_seq[qpos] == "C":
				qualC_rev[refpos] += rev_quality[qpos]
				countsC_rev[refpos] += 1
			elif rev_seq[qpos] == "G":
				qualG_rev[refpos] += rev_quality[qpos]
				countsG_rev[refpos] += 1
			elif rev_seq[qpos] == "T":
				qualT_rev[refpos] += rev_quality[qpos]
				countsT_rev[refpos] += 1

meanQualA_fw = [round(x/y,1) for x, y in zip(qualA_fw, countsA_fw)]
meanQualC_fw = [round(x/y,1) for x, y in zip(qualC_fw, countsC_fw)]
meanQualG_fw = [round(x/y,1) for x, y in zip(qualG_fw, countsG_fw)]
meanQualT_fw = [round(x/y,1) for x, y in zip(qualT_fw, countsT_fw)]

countsA_fw = [ int(round(elem)) for elem in countsA_fw ]
countsC_fw = [ int(round(elem)) for elem in countsC_fw ]
countsG_fw = [ int(round(elem)) for elem in countsG_fw ]
countsT_fw = [ int(round(elem)) for elem in countsT_fw ]

meanQualA_rev = [round(x/y,1) for x, y in zip(qualA_rev, countsA_rev)]
meanQualC_rev = [round(x/y,1) for x, y in zip(qualC_rev, countsC_rev)]
meanQualG_rev = [round(x/y,1) for x, y in zip(qualG_rev, countsG_rev)]
meanQualT_rev = [round(x/y,1) for x, y in zip(qualT_rev, countsT_rev)]

countsA_rev = [ int(round(elem)) for elem in countsA_rev ]
countsC_rev = [ int(round(elem)) for elem in countsC_rev ]
countsG_rev = [ int(round(elem)) for elem in countsG_rev ]
countsT_rev = [ int(round(elem)) for elem in countsT_rev ]

# Allele Counts
bam = pysam.AlignmentFile(bamfile, "rb")

if(emit_base_qualities == "True"):
	writeSparseMatrix4("A", countsA_fw, meanQualA_fw, countsA_rev, meanQualA_rev)
	writeSparseMatrix4("C", countsC_fw, meanQualC_fw, countsC_rev, meanQualC_rev)
	writeSparseMatrix4("G", countsG_fw, meanQualG_fw, countsG_rev, meanQualG_rev)
	writeSparseMatrix4("T", countsT_fw, meanQualT_fw, countsT_rev, meanQualT_rev)
else:
	writeSparseMatrix2("A", countsA_fw, countsA_rev)
	writeSparseMatrix2("C", countsC_fw, countsC_rev)
	writeSparseMatrix2("G", countsG_fw, countsG_rev)
	writeSparseMatrix2("T", countsT_fw, countsT_rev)

zipped_list = zip(list(countsA_fw),list(countsC_fw),list(countsG_fw),list(countsT_fw), list(countsA_rev),list(countsC_rev),list(countsG_rev),list(countsT_rev))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)