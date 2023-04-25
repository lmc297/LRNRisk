#!/usr/bin/env python3

import argparse
import sys, os, glob

def get_gtdb(f):

	# get GTDB species
	# assumes there is one genome in the GTDB-Tk output file
	with open(f, "r") as infile:
		for line in infile:
			if "user_genome" not in line:
				gname = line.split("\t")[0].strip()
				tax = line.split("\t")[1].strip()
				tax = tax.split(";")[-1].strip()
				# split on GTDB species tag
				tax = tax.split("s__")[1].strip()
				if len(tax) == 0:
					tax = "(Unknown Species)"
	return(tax) 

def get_blast(f):

	# reads genes detected via BLAST
	# BLAST header is as follows:
	# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen
	d = {}
	with open(f, "r") as infile:
		for line in infile:
			gene = line.split("\t")[0].strip()
			contig = line.split("\t")[1].strip()
			pid = line.split("\t")[2].strip()
			alen = line.split("\t")[3].strip()
			e = line.split("\t")[-4].strip()
			qlen = line.split("\t")[-1].strip()
			# calculate query coverage by dividing alignment length by query length
			qcov = round(float(alen)/float(qlen)*100.0, 2)

			if gene not in d.keys():
				d[gene] = []
			d[gene].append(line.strip() + "\t" + str(qcov))

	return d
	
def get_blacklist(v, b):

	# identify high-risk isolates based on blacklisted genes
	# blacklisted genes file contains two columns:
	# column 0 = the gene name as it appears in the gene database
	# column 1 = the reason why the gene was blacklisted, which will be reported
	# e.g., "ANTHRAX TOXIN"
	bdict = {}
	with open(b, "r") as infile:
		for line in infile:
			gene = line.split("\t")[0].strip()
			val = line.split("\t")[1].strip()
			bdict[gene] = val

	blacklist_present = {}
	for key in v.keys():
		if key in bdict.keys():
			val = bdict[key]
			blacklist_present[key] = val

	return blacklist_present

def gene_dist(f, blast, gtdb):

	# get within-species prevalence of genes
	# for virulence factors (VFs): uses VFDB VFs detected via ABRicate's VFDB db
	# for AMR genes: uses AMR genes detected via ABRicate's ResFinder db
	# for VFs and AMR genes: genes were detected via ABRicate XXX
	# minimum nucleotide identity and coverage values >=80%
	# total of 61,161 genomes queried
	# takes VFDB or AMR gene distribution file as input (f)
	# BLAST file of VFDB or AMR genes (blast)
	# GTDB species (gtdb)
	
	# create dictionaries based on gene distribution
	d = {}
	annd = {}
	gtdbd = {}
	with open(f, "r") as infile:
		for line in infile:
			tax = line.split("\t")[0].strip()
			tax = tax.split("s__")[1].strip()
			if len(tax) == 0:
				tax = "(Unknown Species)"
			gene = line.split("\t")[1].strip()
			ann = line.split("\t")[-1].strip()
			denom = line.split("\t")[3].strip()
			d[tax + "___" + gene] = line.strip()
			annd[gene] = ann
			gtdbd[tax] = denom
	
	# parse BLAST results
	finallines = []
	for key in blast.keys():
		blastval = blast[key]
		for bv in blastval:
			testkey = gtdb + "___" + key
			if testkey in d.keys() and gtdb != "(Unknown Species)":
				taxval = d[testkey]
				tax = taxval.split("\t")[0].strip()
				tax = tax.split("s__")[1].strip()
				if len(tax) == 0:
					tax = "(Unknown Species)"
				gene = taxval.split("\t")[1].strip()
				pres = taxval.split("\t")[2].strip()
				denom = taxval.split("\t")[3].strip()
				perc = taxval.split("\t")[4].strip()
				perc = str(round(float(perc), 2))
				ann = taxval.split("\t")[-1].strip()
				freetext = "Gene {0} has been detected in {1}% of {2} genomes ({3} of {4} genomes queried)".format(gene, perc, tax, pres, denom)
			elif gtdb != "(Unknown Species)":
				ann = annd[key]
				denom = gtdbd[gtdb]
				freetext = "WARNING: Gene {0} ({1}) has never been detected in species {2} (n = {3} genomes queried)! Interpret with caution!".format(key, ann, gtdb, denom)
			else:
				ann = annd[key]
				freetext = "WARNING: Genome belongs to an undescribed species. Interpret with caution!"
			finalline = bv + "\t" + ann + "\t" + freetext
			finallines.append(finalline)
	return finallines

def print_output(blacklist, vfdist, amrdist, fname, gtdb):

	# print LRNRisk output files
	# takes detected blacklisted genes as input (blacklist)
	# takes distribution of virulence factors as input (vfdist)
	# takes distribution of AMR genes as input (amrdist)
	# takes path to output directory plus a prefix as input (fname)
	# takes GTDB species as input (GTDB)
	prefix = fname.split("/")[-1].strip()

	# blacklist results
	header_blacklist = ["Blacklisted Gene", "Reason", "Risk Category"]
	header_blacklist = "\t".join(header_blacklist).strip()
	with open(fname + "_lrnrisk_blacklist.tsv", "a") as outfile:
		
		print(header_blacklist, file = outfile)
		if len(blacklist.keys()) == 0:
			# print this if no blacklisted genes are detected
			print("\t".join(["(No blacklisted genes detected)", "NA", "Not high risk"]), file = outfile)
		else:
			# print this if blacklisted genes are detected
			# print a table with one row per detected blacklisted gene
			for key in blacklist.keys():
				val = blacklist[key]
				print("\t".join([key, val, "HIGH RISK"]), file = outfile)


	# VFDB results	
	header_vfdb = ["Gene", "Contig", "% Identity", "% Coverage", "E-Value", "Annotation", "Comparison to Publicly Available Genomes"]
	header_vfdb = "\t".join(header_vfdb).strip()
	with open(fname + "_lrnrisk_vfdb.tsv", "a") as outfile:
		print(header_vfdb, file = outfile)
		if len(vfdist) == 0:
			# print this if no VFs detected
			print("\t".join(["(No VFs Detected)"]*7), file = outfile)
		else:
			# print table of VFs if VFs detected
			for vline in vfdist:
				# blast_header = ["Gene", "Contig", "Percent (%) Nucleotide Identity", "Alignment Length", "Mismatches", "Gaps", "Query Start", "Query End", "Subject Start", "Subject End", "E-Value", "Bit Score",  "Identical Matches", "Query Length"]
				# lc_header = ["Query Coverage", "Annotation", "Comparison to Publicly Available Genomes"]
				vgene = vline.split("\t")[0].strip()
				vcontig = vline.split("\t")[1].strip()
				vid = vline.split("\t")[2].strip()
				vcov = vline.split("\t")[-3].strip()
				veval = vline.split("\t")[-7].strip()
				vann = vline.split("\t")[-2].strip()
				vnotes = vline.split("\t")[-1].strip()
				vfinal = [vgene, vcontig, vid, vcov, veval, vann, vnotes]
				vfinal = "\t".join(vfinal).strip()	
				print(vfinal, file = outfile)


	# AMR results
	with open(fname + "_lrnrisk_amr.tsv", "a") as outfile:
		print(header_vfdb, file = outfile)
		if len(amrdist) == 0:
			# print this if no AMR genes detected
			print("\t".join(["(No AMR Genes Detected)"]*7), file = outfile)
		else:
			# print this if AMR genes detected
			for aline in amrdist:
				# blast_header = ["Gene", "Contig", "Percent (%) Nucleotide Identity", "Alignment Length", "Mismatches", "Gaps", "Query Start", "Query End", "Subject Start", "Subject End", "E-Value", "Bit Score",  "Identical Matches", "Query Length"]
				# lc_header = ["Query Coverage", "Annotation", "Comparison to Publicly Available Genomes"]
				agene = aline.split("\t")[0].strip()
				acontig = aline.split("\t")[1].strip()
				aid = aline.split("\t")[2].strip()
				acov = aline.split("\t")[-3].strip()
				aeval = aline.split("\t")[-7].strip()
				aann = aline.split("\t")[-2].strip()
				anotes = aline.split("\t")[-1].strip()
				afinal = [agene, acontig, aid, acov, aeval, aann, anotes]
				afinal = "\t".join(afinal).strip()
				print(afinal, file = outfile)


def main():

	# lrnrisk_prototype arguments

	parser = argparse.ArgumentParser(prog = "lrnrisk_prototype.py", usage = "lrnrisk_prototype.py -v </path/to/vfdb_output.tab> -a </path/to/pima_amr_output.tab> -o </path/to/output/directory/prefix> [other options]")

	parser.add_argument("-g", "--gtdb", help = "Path to gtdbtk.bac120.summary.tsv, the tab-separated file containing GTDB-Tk results", nargs = 1, required = True)

	parser.add_argument("-v", "--virulence", help = "Path to tab-separated file containing VFDB virulence factors detected via BLAST", nargs = 1, required = True)

	parser.add_argument("-a", "--amr", help = "Path to tab-separated file containing Pima AMR determinants detected via BLAST", nargs = 1, required = True)

	parser.add_argument("-b", "--blacklist", help = "Path to tab-separated file containing blacklisted high-risk virulence factors", nargs = 1, required = True)

	parser.add_argument("-dv", "--vf_distribution", help = "Path to tab-separated file containing virulence factor distribution", nargs = 1, required = True)

	parser.add_argument("-da", "--amr_distribution", help = "Path to tab-separated file containing AMR determinant distribution", nargs = 1, required = True)

	parser.add_argument("-o", "--output", help = "Path to desired output directory, plus a prefix for the output files", nargs = 1, required = True)

	# parse arguments and run pipeline
	args = parser.parse_args()
	my_gtdb = get_gtdb(args.gtdb[0])
	my_vf = get_blast(args.virulence[0])
	my_amr = get_blast(args.amr[0])
	my_blacklist = get_blacklist(my_vf, args.blacklist[0])
	my_vfdist = gene_dist(args.vf_distribution[0], my_vf, my_gtdb)
	my_amrdist = gene_dist(args.amr_distribution[0], my_amr, my_gtdb)
	print_output(my_blacklist, my_vfdist, my_amrdist, args.output[0], my_gtdb)

if __name__ == "__main__":

	# run lrnrisk_prototype

	main()
