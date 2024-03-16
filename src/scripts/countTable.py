from collections import defaultdict
import re
import argparse 

def main():
	options = setup()
	file_list = options.file_list.split(",")
	print(file_list)
	gene_dict, sample_list = count_reads(file_list)
	save_output(gene_dict, sample_list, options.output_file)

#########
# Setup #
#########

def setup():
	"""
	Gets command line arguments and returns a Namespace object
	"""

	parser = argparse.ArgumentParser()

    # House keeping to read in arguments from the command line
	parser.add_argument("-f", "--file_list", dest = "file_list",
		help = "The list of input files that came from featureCounts, separate individual files with a comma and no space.",
		default = "none",
		action = "store",
		metavar = "\b")

	parser.add_argument("-o", "--output", dest = "output_file",
		help = "The path to the output file",
		default = "counts.txt",
		action = "store",
		metavar = "\b")

	args = parser.parse_args()

	return(args)

def count_reads(file_list):
	"""
	Opens each feture count file and adds to a dictionary the number of reads from the sample that was counted
	"""

	gene_dict = defaultdict(list)
	sample_list = list()

	file_path = "featureCounts/"

	for i in file_list:
		sample_name = re.sub(rf'{file_path}', '', i)
		sample_name = re.sub(r'_countsOutput', '', sample_name)
		sample_list.append(sample_name)
		with open(i, "r") as countFile:
			for line in countFile:
				line = line.strip().split("\t")
				if "ENS" in line[0]:
					ens_id = line[0]
					gene = line[6]
					count = line[8]
					gene_ens = gene + "_" + ens_id
					gene_dict[gene_ens].append(count)

	return(gene_dict, sample_list)

def save_output(gene_dict, sample_list, output_file):
	"""
	Creates an output file where columns are samples and rows are genes
	"""

	with open(output_file, "w") as count_file:
		count_file.write("gene" + "\t" + "\t".join(sample_list) + "\n")
		for gene in gene_dict:
			count_file.write(gene + "\t" + "\t".join(gene_dict[gene]) + "\n")

if __name__ == "__main__":
	main()