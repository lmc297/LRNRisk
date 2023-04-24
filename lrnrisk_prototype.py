#!/usr/bin/env python3

import argparse
import sys, os, glob

def get_gtdb(f):

	# get GTDB species
	with open(f, "r") as infile:
		for line in infile:
			if "user_genome" not in line:
				gname = line.split("\t")[0].strip()
				tax = line.split("\t")[1].strip()
				tax = tax.split(";")[-1].strip()
				tax = tax.split("s__")[1].strip()
				if len(tax) == 0:
					tax = "(Unknown Species)"
	return(tax) 

def get_blast(f):

	# reads genes detected via BLAST
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
			qcov = round(float(alen)/float(qlen)*100.0, 2)

			if gene not in d.keys():
				d[gene] = []
			d[gene].append(line.strip() + "\t" + str(qcov))

	return d
	
def get_blacklist(v, b):

	# identify high-risk isolates based on blacklisted genes

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
				freetext = "Gene <i>{0}</i> has been detected in {1}% of <i>{2}</i> genomes ({3} of {4} genomes queried)".format(gene, perc, tax, pres, denom)
			elif gtdb != "(Unknown Species)":
				ann = annd[key]
				denom = gtdbd[gtdb]
				freetext = "WARNING: Gene <i>{0}</i> ({1}) has never been detected in species <i>{2}</i> (<i>n</i> = {3} genomes queried)! Interpret with caution!".format(key, ann, gtdb, denom)
			else:
				ann = annd[key]
				freetext = "WARNING: Genome belongs to an undescribed species. Interpret with caution!"
			finalline = bv + "\t" + ann + "\t" + freetext
			finallines.append(finalline)
	return finallines

def print_html(blacklist, vfdist, amrdist, fname, gtdb):

	# print LRNRisk HTML report
	
	prefix = fname.split("/")[-1].strip()
	prefix = prefix.split(".html")[0].strip()

	html_header = """<!DOCTYPE html>
<html>
<head>
<style>
img {
  width: 25%; 
}
hr {
  height: 4px;
  background-color:#000000;
  opcaity: 0.5;
}
body {
  font-family: Arial;
}
#results {
  font-family: Arial, Helvetica, sans-serif;
  border-collapse: collapse;
  width: 100%;
}

#results td, #results th {
  border: 1px solid #ddd;
  padding: 8px;
}

#results tr:nth-child(even){background-color: #f2f2f2;}

#results tr:hover {background-color: #ddd;}

#results th {
  padding-top: 12px;
  padding-bottom: 12px;
  text-align: left;
  background-color: #00B0F0;
  color: white;
}

#blacklistresults {
  font-family: Arial, Helvetica, sans-serif;
  border-collapse: collapse;
  width: 100%;
}

#blacklistresults td, #blacklistresults th {
  border: 1px solid #ddd;
  padding: 8px;
}

#blacklistresults tr:nth-child(even){background-color: #f2f2f2;}

#blacklistresults tr:hover {background-color: #ddd;}

#blacklistresults th {
  padding-top: 12px;
  padding-bottom: 12px;
  text-align: left;
  background-color: #C43E96;
  color: white;
}

.sidebar {
  margin: 0;
  padding: 0;
  width: 200px;
  background-color: #f1f1f1;
  position: fixed;
  height: 100%;
  overflow: auto;
}

.sidebar a {
  display: block;
  color: black;
  padding: 16px;
  text-decoration: none;
}

.sidebar a.active {
  background-color: #00B0F0;
  color: white;
}

.sidebar a:hover:not(.active) {
  background-color: #555;
  color: white;
}

div.content {
  margin-left: 200px;
  padding: 1px 16px;
  height: 1000px;
}

@media screen and (max-width: 700px) {
  .sidebar {
    width: 100%;
    height: auto;
    position: relative;
  }
  .sidebar a {float: left;}
  div.content {margin-left: 0;}
}

@media screen and (max-width: 400px) {
  .sidebar a {
    text-align: center;
    float: none;
  }
}
</style>
</head>
</head>
<body>
<img src="lrnrisk_logo.jpg" alt="LRNRisk Logo">
<hr>

<!-- The sidebar -->
<div class="sidebar">
  <a class="active" href="#home">Home</a> 
  <a href="#blacklist">Blacklisted Genes</a>
  <a href="#gtdb">GTDB Species</a>
  <a href="#vfdb">VFDB Virulence Factors</a>
  <a href="#amr">Pima AMR Genes</a>
  <a href="#disclaimer">Disclaimer</a>
</div>

<!-- Page content -->
<div class="content">
"""

	html_title = '<h1 id="home"><i>LRNRisk</i> Results for {0}</h1>'.format(prefix)

	html_blacklist_true = """<h2 id="blacklist">Blacklisted Genes</h2>
<p align=justify style="background-color:#FDF0F6;
border-radius:4px;
font-size:36px;
padding:15px;
margin:5px;">
<b>ALERT! HIGH-RISK GENOME!</b> The following BLACKLISTED GENES were detected:
</p>
<table id="blacklistresults">
<tr>
<th>Blacklisted Gene</th>
<th>Reason</th>
</tr>
"""

	html_blacklist_true_boxes = """<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> LRNRisk identified one or more genes in the query genome, which have been "blacklisted" by the CDC and/or LRNRisk developers due to their public health risk. <b>Strains harboring one or more blacklisted genes should be treated as potentially dangerous, high-risk pathogens.</b>
</p>

<p align=justify style="background-color:#FDF0F6;
border-radius:4px;
font-size:36px;
padding:15px;
margin:5px;">
<b>ALERT! BLACKLISTED GENE(S) DETECTED! THIS IS A PREDICTED HIGH-RISK GENOME! PROCEED WITH EXTREME CAUTION!</b>
</p>
"""

	html_blacklist_false = """<h3>No blacklisted genes detected. Genome is not predicted to be high-risk.</h3>
<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> LRNRisk did not identify any genes in the query genome, which have been "blacklisted" by the CDC and/or LRNRisk developers due to their public health risk. <b>The query genome is not predicted to be high risk. Proceed with normal caution.</b>
</p>
"""

	gtdb_header = """<h2 id="gtdb">GTDB Species</h2>
<table id="results">
<tr>
<th>Genome</th>
<th>GTDB Species</th>
</tr>
"""

	gtdb_box = """<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> The query genome was assigned to the corresponding species using the Genome Taxonomy Database (GTDB) Toolkit (GTDB-Tk). GTDB-Tk assigns genomes to species within the GTDB taxonomy, a standardized, species-level taxonomy for prokaryotes.
<br>
<b>References:</b>
<br>
&#x2022; Chaumeil, et al. 2020. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7703759/">GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database.</a> <i>Bioinformatics</i> 36(6):1925-1927. doi: 10.1093/bioinformatics/btz848.
<br>
&#x2022; Parks, et al. 2022. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8728215/">GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy</a>. <i>Nucleic Acids Research</i> 50(D1): D785–D794. doi: 10.1093/nar/gkab776.
</p>

<p align=justify style="background-color:#FFFCF4;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>WARNING! Always interpret species names with care.</b> Pathogen misidentifications can result (and have resulted) from misinterpretation of GTDB species labels. <b>GTDB species labels should be considered holistically with other LRNRisk strain attributes to assess risk</b> (e.g., blacklisted genes, virulence factors).
<br>
<b>References:</b>
<br>
&#x2022; Carroll, et al. 2022. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9423903/">Laboratory Misidentifications Resulting from Taxonomic Changes to <i>Bacillus cereus</i> Group Species, 2018–2022.</a> <i>Emerging Infectious Diseases</i> 28(9): 1877–1881. doi: 10.3201/eid2809.220293.
</p>
"""

	gene_header = """<table id="results">
<tr>
<th>Gene</th>
<th>Contig</th>
<th>% Identity</th>
<th>% Coverage</th>
<th>E-Value</th>
<th>Annotation</th>
<th>Comparison to Publicly Available Genomes</th>
</tr>
"""

	
	vfdb_box_true = """<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> One or more genes within the Virulence Factor Database (VFDB) core database were detected in the query genome. According to <a href="http://www.mgc.ac.cn/VFs/main.htm">VFDB</a>, virulence factors are genes "that enable a microorganism to establish itself on or within a host of a particular species and enhance its potential to cause disease. Virulence factors include bacterial toxins, cell surface proteins that mediate bacterial attachment, cell surface carbohydrates and proteins that protect a bacterium, and hydrolytic enzymes that may contribute to the pathogenicity of the bacterium."
<br>
<b>References:</b>
<br>
&#x2022; Liu, et al. 2022. <a href="https://pubmed.ncbi.nlm.nih.gov/34850947/">VFDB 2022: a general classification scheme for bacterial virulence factors.</a> <i>Nucleic Acids Research</i> 50(D1):D912-D917. doi: 10.1093/nar/gkab1107.
<br>
&#x2022; Chen, et al. 2005. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC539962/">VFDB: a reference database for bacterial virulence factors.</a> <i>Nucleic Acids Research</i> 33(Database Issue): D325–D328. doi: 10.1093/nar/gki008.
</p>

<p align=justify style="background-color:#FFFCF4;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>WARNING! Always interpret VFDB results with caution.</b> The presence of VFDB virulence factors does not necessarily indicate increased risk (e.g., VFDB contains stress response genes, genes involved in nutrient metabolism); likewise, novel virulence factors may not be detected. To help assess risk, consider:<br>
1) <b>The function of the virulence factor</b> (see the "Annotation" column in the table above)<br>
2) <b>Other attributes reported by LRNRisk (e.g., blacklisted genes, GTDB species)</b>
"""

	vfdb_box_false = """<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> No genes within the Virulence Factor Database (VFDB) core database were detected in the query genome at the LRNRisk BLAST thresholds. According to <a href="http://www.mgc.ac.cn/VFs/main.htm">VFDB</a>, virulence factors are genes "that enable a microorganism to establish itself on or within a host of a particular species and enhance its potential to cause disease. Virulence factors include bacterial toxins, cell surface proteins that mediate bacterial attachment, cell surface carbohydrates and proteins that protect a bacterium, and hydrolytic enzymes that may contribute to the pathogenicity of the bacterium."
<br>
<b>References:</b>
<br>
&#x2022; Liu, et al. 2022. <a href="https://pubmed.ncbi.nlm.nih.gov/34850947/">VFDB 2022: a general classification scheme for bacterial virulence factors.</a> <i>Nucleic Acids Research</i> 50(D1):D912-D917. doi: 10.1093/nar/gkab1107.
<br>
&#x2022; Chen, et al. 2005. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC539962/">VFDB: a reference database for bacterial virulence factors.</a> <i>Nucleic Acids Research</i> 33(Database Issue): D325–D328. doi: 10.1093/nar/gki008.
</p>
"""

	amr_box_true = """<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> One or more antimicrobial resistance (AMR) determinants were detected within the query genome via the CDC's Pima pipeline. Pima is a bioinformatic pipeline, which was developed by the CDC "to perform rapid <i>de novo</i> genome and plasmid assemblies, to detect AMR variants and genes, and to identify variants (mismatches and indels)."
<br>
<b>References:</b>
<br>
&#x2022; Gargis, et al. 2019. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6751186/">Rapid Detection of Genetic Engineering, Structural Variation, and Antimicrobial Resistance Markers in Bacterial Biothreat Pathogens by Nanopore Sequencing</a> <i>Scientific Reports</i> 9: 13501. doi: 10.1038/s41598-019-49700-1.

<p align=justify style="background-color:#FFFCF4;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>WARNING! Always interpret Pima AMR results with caution.</b> Detected AMR determinants may not confer phenotypic AMR; likewise, novel AMR determinants may not be detected.
"""

	amr_box_false= """<p align=justify style="background-color:#EEF4FB;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>What does this mean?</b> No antimicrobial resistance (AMR) determinants were detected within the query genome via the CDC's Pima pipeline. Pima is a bioinformatic pipeline, which was developed by the CDC "to perform rapid <i>de novo</i> genome and plasmid assemblies, to detect AMR variants and genes, and to identify variants (mismatches and indels)."
<br>
<b>References:</b>
<br>
&#x2022; Gargis, et al. 2019. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6751186/">Rapid Detection of Genetic Engineering, Structural Variation, and Antimicrobial Resistance Markers in Bacterial Biothreat Pathogens by Nanopore Sequencing</a> <i>Scientific Reports</i> 9: 13501. doi: 10.1038/s41598-019-49700-1.
"""

	disclaimer = """<p align=justify style="background-color:#F0F8F5;
border-radius:4px;
font-size:16px;
padding:15px;
margin:5px;">
<b>Disclaimer:</b> LRNRisk aims to rapidly identify potential bacterial pathogens and provide insight into their potential public health risk. However, no tool is perfect, and LRNRisk cannot definitively prove whether an isolate is pathogenic or not. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic potential, and/or misinterpretations (biological, statistical, or otherwise) of LRNRisk results.
"""

	with open(fname, "a") as outfile:
		print(html_header, file = outfile)
		print(html_title, file = outfile)
		
		# blacklist results	
		if len(blacklist.keys()) == 0:
			print(html_blacklist_false, file = outfile)
		else:
			print(html_blacklist_true, file = outfile)
			for key in blacklist.keys():
				val = blacklist[key]
				print("<tr>", file = outfile)
				print("<td><i>" + key + "</i></td>", file = outfile)
				print("<td>" + val + "</td>", file = outfile)
				print("</tr>", file = outfile)
			print("</table>", file = outfile)
			print(html_blacklist_true_boxes, file = outfile)

		# GTDB results
		print(gtdb_header, file = outfile)
		print("<tr>", file = outfile)
		print("<td>{0}</td>".format(prefix), file = outfile)
		print("<td><i>{0}</i></td>".format(gtdb), file = outfile)
		print("</tr>", file = outfile)
		print("</table>", file = outfile)
		print(gtdb_box, file = outfile)

		# VFDB results
		print('<h2 id="vfdb">VFDB Virulence Factors</h2>', file = outfile)
		if len(vfdist) == 0:
			print("<h3>No VFDB virulence factors detected.</h3>", file = outfile)
			print(vfdb_box_false, file = outfile)
		else:
			print(gene_header, file = outfile)
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
				vfinal = ["<i>" + vgene + "</i>", vcontig, vid, vcov, veval, vann, vnotes]
				print("<tr>", file = outfile)
				for vfi in vfinal:
					print("<td>" + vfi.strip() + "</td>", file = outfile)
				print("</tr>", file = outfile)
			print("</table>", file = outfile)
			print(vfdb_box_true, file = outfile)


		# AMR results
		print('<h2 id="amr">Pima AMR Genes</h2>', file = outfile)
		if len(amrdist) == 0:
			print("<h3>No Pima AMR genes detected.</h3>", file = outfile)
			print(amr_box_false, file = outfile)
		else:
			print(gene_header, file = outfile)
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
				afinal = ["<i>" + agene + "</i>", acontig, aid, acov, aeval, aann, anotes]
				print("<tr>", file = outfile)
				for afi in afinal:
					print("<td>" + afi.strip() + "</td>", file = outfile)
				print("</tr>", file = outfile)
			print("</table>", file = outfile)
			print(amr_box_true, file = outfile)	

		print('<h2 id="disclaimer">Disclaimer</h2>', file = outfile)	
		print(disclaimer, file = outfile)		
		print("</div>", file = outfile)
		print("</body>", file = outfile)
		print("</html>", file = outfile)

def main():

	# lrnrisk_prototype arguments

	parser = argparse.ArgumentParser(prog = "lrnrisk_prototype.py", usage = "lrnrisk_prototype.py -v </path/to/vfdb_output.tab> -a </path/to/pima_amr_output.tab> -o </path/to/output/directory/> [other options]")

	parser.add_argument("-g", "--gtdb", help = "Path to gtdbtk.bac120.summary.tsv, the tab-separated file containing GTDB-Tk results", nargs = 1, required = True)

	parser.add_argument("-v", "--virulence", help = "Path to tab-separated file containing VFDB virulence factors detected via BLAST", nargs = 1, required = True)

	parser.add_argument("-a", "--amr", help = "Path to tab-separated file containing Pima AMR determinants detected via BLAST", nargs = 1, required = True)

	parser.add_argument("-b", "--blacklist", help = "Path to tab-separated file containing blacklisted high-risk virulence factors", nargs = 1, required = True)

	parser.add_argument("-dv", "--vf_distribution", help = "Path to tab-separated file containing virulence factor distribution", nargs = 1, required = True)

	parser.add_argument("-da", "--amr_distribution", help = "Path to tab-separated file containing AMR determinant distribution", nargs = 1, required = True)

	parser.add_argument("-o", "--output", help = "Path to desired output HTML file", nargs = 1, required = True)

	args = parser.parse_args()
	my_gtdb = get_gtdb(args.gtdb[0])
	my_vf = get_blast(args.virulence[0])
	my_amr = get_blast(args.amr[0])
	my_blacklist = get_blacklist(my_vf, args.blacklist[0])
	my_vfdist = gene_dist(args.vf_distribution[0], my_vf, my_gtdb)
	my_amrdist = gene_dist(args.amr_distribution[0], my_amr, my_gtdb)
	print_html(my_blacklist, my_vfdist, my_amrdist, args.output[0], my_gtdb)

if __name__ == "__main__":

	# run lrnrisk_prototype

	main()
