#!/usr/bin/python
import sys, os, re, getopt, glob, subprocess, os.path, numpy as np
import cPickle as pickle
import timeit
import itertools

now = timeit.default_timer()

Usage = """
Usage:  ./Metagenome_Assembly.py -o ASSEMBLED -i READ_FILES_LIST.tsv -S 1,2 -A NexteraPE -F Y -M IDBA

REQUIRED ARGUMENTS:
		-o	output directory
	
                -i	your input file names in the form of a .tsv file with three columns : "SampleID", "Forward PE" and "Reverse PE"
			(If reads are already interleaved, do not provide a third column) 
			(files can be compressed as either ".gz" or ".bz2")

		-S	Specify which processing steps to perform. If you only perform assembly, make sure input .tsv has three columns: "SampleID", "PE Reads", "SE Reads"
			(1: Quality Control, 2: Assembly)

		-p	Specify Number of Processors to Use When Possible

DEPENDENCIES:
		-FastX Toolkit (built with v.0.7)
		-khmer (built with v.1.0)
		-Any of : 	
				-IDBA (built with v.1.1.1)
				-RAY-meta (built with v.v2.3.1)
				-MEGAHIT (built with v1.0.4-beta)
				-metaSPAdes (built with v.3.10.1)

OPTIONAL SOFTWARE:
		-Trimmomatic (built with v.0.36)
		-FLASH (built with v.2.0)
                -EA-UTILS
	
OPTIONAL:
		## Related To Adapter Removal
		-A 	Specify whether to Remove Illumina Adapters using Trimmomatic
			(Must specify which adapter library to use, i.e. "NexteraPE","TruSeq2" or "TruSeq3")

                -W <Y>  Specify QC filtering method (either FastX Toolkig or EA-UTILS)
                        (DEFAULT: FastX)

		## IF FASTX TOOLKIT: Additional Specification Related to Quality Trimming
		-Q	Specify Phred Format  [default: Phred33]

		-q	Specify Minimum Quality Score for Quality Trimming [default: 30]

		-l	Specify Percentage of Read That Must Be Above Quality Score [default: 50]

		## Related to Improving Assembly
		-F <Y>	Use FLASH to extend paired end reads

		## Related to Assembly
		-M	Specify which assembler to use by either name or number
			(1: IDBA, 2: RAYmeta, 3: MEGAHIT, 4: metaSPAdes, 5: SOAPdenovo2)

UTILITY:
This script will perform the basic quality control (i.e. remove artefacts from Illumina adapters and quality trim)and assemble metagenomes using a variety of methods."

NOTES:
- The QC Pipeline was appropriated from the khmer workflow described here: <https://khmer-protocols.readthedocs.org/en/latest/metagenomics/1-quality.html>
- File extensions accurately denote what kind of analysis has been done to the file:

.pe - Interleaved Paired end reads file
.se - Unpaired reads file
.fq - Fastq format
.qc - Quality trimmed
.tr - Adapter trimmed
.fl - Flash merged

- Logging is set ON to by default. Editorializing is logged to "run.log", while all commands executed to your files in "command.log". This can be shut off by changing the section of script "LOG = 'ON'" to 'OFF'


Usage:  ./Metagenome_Assembly.py -o ASSEMBLED -i READ_FILES_LIST.tsv -S 1,2 -A NexteraPE -F Y -M IDBA


Note: This script will need minor modification for how input is read if you're providing three files (i.e. two paired read files + orphaned reads)
      You will also need to specify the specific location of the Trimmomatic jar file and adapter
"""


if len(sys.argv)<2:
        print Usage
        sys.exit()

# Initialize all containers for input 
LIST=''
PROCESSES=[]
ADAPTER=''
QUAL_FORMAT=''
MIN_QUAL=''
MIN_LENGTH=''
FLASH=''
ASSEMBLER=''
OUTPUT=''
EA_UTILS=''
PROCESSORS=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"o:i:A:Q:q:l:F:S:M:W:p:L:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-o':
        OUTPUT= a
    if o == '-i':
        LIST= a
    if o == '-A':
        ADAPTER= a
    if o == '-Q':
        QUAL_FORMAT= a
    if o == '-q':
        MIN_QUAL= a
    if o == '-l':
        MIN_LENGTH= a
    if o == '-F':
        FLASH= a
    if o == '-S':
        PROCESSES= a
    if o == '-M':
        ASSEMBLER= a
    if o == '-W':
        EA= a
    if o == '-p':
        PROCESSORS= a
    if o == '-L':
        LOG= a

## Set Defaults unless user specified otherwise
if len(QUAL_FORMAT) == 0:
	QUAL_FORMAT = "33"

if len(MIN_QUAL) == 0:
	MIN_QUAL = "30"

if len(MIN_LENGTH) == 0:
	MIN_LENGTH = "50"

try:
	len(LOG)
	LOG = 'OFF'
except:
	LOG = 'ON'

if LOG == 'ON':
        print "\n-- You are Logging a Description of Analyses and All Commands --"
	log = open("run.log", "w")
	command = open("command.log", "w")

if len(OUTPUT)>0:
        if os.path.exists('./' + OUTPUT):
                print "\n-- Output Folder Exists - Caution: Files May Be Over-written --"
        else:
                os.mkdir(OUTPUT)


if LOG == 'ON':
	log.write("\n\nThe following defaults were used:\n\nFastX Quality Trim\nPhred Format: "+QUAL_FORMAT+"\nMinimum Quality Score: "+MIN_QUAL+"\nMinimum Length: "+MIN_LENGTH+"\n")

###################
## Define Functions
###################

def UnZip(file):
	#Unzip if necessary
	if re.search(".bz2", file):
		os.system(' '.join([
			"tar -xjfv",
			file
		]))

		# Rename file
        	file = re.sub(".bz2","", file)

	if re.search(".gz", file):
		os.system(' '.join([
			"pigz -d -p",
			PROCESSORS,
			file,
		]))

		# Remove zip extension from name values
		file = re.sub(".gz","", file)

	return file

def Zip(file):
	os.system(' '.join([
		"pigz -p",
		PROCESSORS,
		file,
	]))
	

def Interleave(sampleID, file):
	os.system(' '.join([
		"interleave-reads.py",
		file,
		interleave,
		">",
		"./"+OUTPUT+"/"+sampleID+".pe.fq"
		]))

	if LOG == 'ON':
		command.write(' '.join([
        	        "\ninterleave-reads.py",
			file,
       	        	interleave,
                	">",
	                "./"+OUTPUT+"/"+sampleID+".pe.fq"
               	]))
			
	file = "./"+OUTPUT+"/"+sampleID+".pe.fq"

	return file


def Convert_FQ(x):
	os.system(' '.join([
		"fastq_to_fasta -n -i",
		x,
		"-o",
		re.sub("fq|fastq","fasta",x)
	]))
	
	if LOG == 'ON':
		command.write(' '.join([
			"\nfastq_to_fasta -n -i",
			x,
                        "-o",
       	                re.sub("fq|fastq","fasta",x)
                ]))
			
	x = re.sub("fq|fastq","fasta",x)

	return x

		
def FastX(IN,OUT):
        os.system(' '.join([
                "fastq_quality_filter",
                "-Q"+QUAL_FORMAT,
                "-q"+MIN_QUAL,
                "-p"+MIN_LENGTH,
                "-i",
                IN,
                ">",
                OUT
        ]))

        if LOG == 'ON':
                command.write(' '.join([
                        "\nfastq_quality_filter",
                        "-Q"+QUAL_FORMAT,
                        "-q"+MIN_QUAL,
                        "-p"+MIN_LENGTH,
                        "-i",
                        IN,
                        ">",
                        OUT
                ]))


def EAUTILS(IN,OUT):

        dummy = open("adapter-file.fq", "w")
        dummy.close()

        os.system(' '.join([
                "fastq-mcf",
                "adapter-file.fq",
                IN,
                "-o",
                OUT,
                ">>",
                "./" + OUTPUT + "/" + sampleID + "_fastq-mcf.log",
                "2>&1"
        ]))

        if LOG == 'ON':
                command.write(' '.join([
                        "\nfastq-mcf",
                        "adapter-file.fq",
                        IN,
                        "-o",
                        OUT,
                        ">>",
                        "./" + OUTPUT + "/" + sampleID + "_fastq-mcf.log",
                        "2>&1"
                ]))

        os.system('rm ./adapter-file.fq')


def split_paired(file):
	os.system(' '.join([
		"extract-paired-reads.py",
		file
	]))

	if LOG == 'ON':
		command.write("\nextract-paired-reads.py "+file)

	name = file.split("/")[len(file.split("/"))-1]
	return name+".pe", name+".se"


def flash(file):
	if LOG == 'ON':
		log.write("\nYou Selected To Merge Paired Reads Using FLASH.\n")

	os.system(' '.join([
		"flash --interleaved",
		file,
	]))

	if LOG == 'ON':
		command.write("\nflash --interleaved -z "+file)

	os.system(' '.join([
		"rm *.hist*",
		file,
	]))


	return "out.extendedFrags.fastq","out.notCombined.fastq"


def trimmomatic(file1, file2):
	print "-- You Have Selected To Trim Illumina Adapters Using the Following Adapter File: "+ADAPTER

	if LOG == 'ON':
		log.write("\nYou performed adapter trimming with"+ADAPTER+"as the adapter type.\n")

	## Run Trimmomatic for Adapter Removal
	os.system(' '.join([
		"java -jar /home/roli/Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE",
		file1,
		file2,
		"./TRIM/s1_pe ./TRIM/s1_se ./TRIM/s2_pe ./TRIM/s2_se",
		"ILLUMINACLIP:/home/roli/Software/Trimmomatic-0.36/adapters/"+ADAPTER+"-PE.fa:2:30:10"
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"java -jar /home/roli/Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE",
			file1,
			file2,
			"./TRIM/s1_pe ./TRIM/s1_se ./TRIM/s2_pe ./TRIM/s2_se",
			"ILLUMINACLIP:/home/roli/Software/Trimmomatic-0.36/adapters/"+ADAPTER+"-PE.fa:2:30:10"
		]))

	## Interleave Reads from the two paired ends files
	os.system(' '.join([
		"interleave-reads.py ./TRIM/s?_pe",
		"> ./TRIM/combined.pe.fq"
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"\ninterleave-reads.py ./TRIM/s?_pe",
			"> ./TRIM/combined.pe.fq"
		]))

	## Combined Orphaned Reads
	os.system(' '.join([
		"cat",
		"./TRIM/s1_se.trim",
		"./TRIM/s2_se.trim",
		"> ./TRIM/combined.se.fq"
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"cat",
			"./TRIM/s1_se.trim",
			"./TRIM/s2_se.trim",
			"> ./TRIM/combined.se.fq"
		]))

	pe = "./TRIM/combined.pe.fq"
	se = "./TRIM/combined.se.fq"

	return pe, se


def idba(pe, sampleID):
	print "\n\n-----Using IDBA-----\n"

	if LOG == 'ON':
		log.write("\nNOTE: IDBA-UD does NOT make use of unpaired reads.")

	os.system(' '.join([
		"idba_ud",
		"-o",
		"./"+OUTPUT+"/ASSEMBLY/IDBA/"+sampleID+"/",
		"-r",
		pe,
		"--pre_correction"
	]))		

	if LOG == 'ON':
		command.write(' '.join(["\nidba_ud","-o","./"+OUTPUT+"/ASSEMBLY/IDBA/"+sampleID+"/","-r",pe,"--pre_correction"])+"\n")		


def ray(pe, se, sampleID):
	os.system(' '.join([
		"mpiexec",
		"-n",
		PROCESSORS,
		"Ray",
		"-o",
		"./"+OUTPUT+"/ASSEMBLY/RAYMETA/"+sampleID,
		"-i",
		pe,
		"-s",
		se,
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"\nmpiexec",
			"-n",
			PROCESSORS,
			"Ray",
			"-o",
			"./"+OUTPUT+"/ASSEMBLY/RAYMETA/"+sampleID,
			"-i",
			pe,
			"-s",
			se,
		]))

	
def megahit(pe, se, sampleID):
	os.system(' '.join([
		"megahit",
		"--12",
		pe,
		"-r",
		se,	
		"-o",
		"./"+OUTPUT+"/ASSEMBLY/MEGAHIT/"+sampleID,
		"--num-cpu-threads",
		PROCESSORS
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"\nmegahit",
			"--12",
			pe,
			"-r",
			se,	
			"-o",
			"./"+OUTPUT+"/ASSEMBLY/MEGAHIT/"+sampleID,
			"--num-cpu-threads",
			PROCESSORS
		]))

	
def metaSPAdes(pe, se, sampleID):
	os.system(' '.join([
		"spades.py",
		"--pe1-12",
		pe,
		"-s",
		se,	
		"--meta",
		"-o",
		"./"+OUTPUT+"/ASSEMBLY/METASPADES/"+sampleID,
		"--threads",
		PROCESSORS
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"spades.py",
			"--pe1-12",
			pe,
			"-s",
			se,	
			"--meta",
			"-o",
			"./"+OUTPUT+"/ASSEMBLY/METASPADES/"+sampleID,
			"--threads",
			PROCESSORS
		]))
	

####################################
# IMPORT FILE LIST AND HOUSE KEEPING
####################################

## IMPORT FILE NAMES into DICTIONARY
FILE_DICT={}

with open(LIST) as f:
	next(f)

	for line in f:
		line = line.strip("\r\n")
		line = line.split("\t")

		try:
			FILE_DICT[line[0]] = [line[1], line[2]]

		except:
			FILE_DICT[line[0]] = [line[1], 'interleaved']

if LOG == 'ON':
	log.write("Your Input Files Were:\n\n")
	log.write(str(FILE_DICT))
	command.write("Your Input Files Were:\n\n")
	command.write(str(FILE_DICT))

#########################
##QUALITY FILTERING STEPS
#########################

if re.search("1", PROCESSES):
	print "\n-- Performing Quality Trimming --\n\n"

	#Make Temporary Trim Folder
       	if os.path.exists("TRIM"):
		pass
	else:
                os.mkdir("TRIM")

	for sampleID, value in FILE_DICT.iteritems():
		file1 = value[0]
		file2 = value[1]

		## Unzip
		if file2 == 'interleaved':
			file1 = UnZip(file1)

		else:
			file1 = UnZip(file1)
			file2 = UnZip(file2)

		## Run TRIMMOMATIC
		if ADAPTER:
			if file2 == 'interleaved':
				file1_mod,file2_mod = split_paired(file1)
				pe, se = trimmomatic(file1_mod,file2_mod)
			else:
				pe, se = trimmomatic(file1,file2)

			## Run Quality Filtering
		        try:
				EA
				EAUTILS(pe,"./TRIM/combined-trim.fq")
				EAUTILS(se,"./TRIM/se.trim")

		        except NameError: #Or DEFAULT: FastX Toolkit
				FastX(pe,"./TRIM/combined-trim.fq")
				FastX(se,"./TRIM/se.trim")
		else:
			## Run Quality Filtering
		        try:
				EA
				EAUTILS(file1,"./TRIM/combined-trim.fq")

		        except NameError: #Or DEFAULT: FastX Toolkit
				FastX(file1,"./TRIM/combined-trim.fq")

			pe, se = split_paired("./TRIM/combined-trim.fq")

			#Re-locate Files for consistent processing
			os.system(' '.join([
				"mv",
				pe,
				"./TRIM/combined-trim.fq"
			]))

			os.system(' '.join([
				"mv",
				se,
				"./TRIM/se.trim.fq"
			]))

			# QC Second Raw Read file
		        if file2 != 'interleaved':
			        if EA:
					EAUTILS(file2,"./FOO.trim")

				else:
					FastX(file2,"./FOO.trim")

				pe, se = split_paired("./FOO.trim.fq")

				os.system(' '.join([
					"cat",
					pe,
					"./TRIM/combined-trim.fq",
					">>",
					"./FOO.fq"
				]))

				os.system(' '.join([
					"mv",
					"./FOO.fq",
					"./TRIM/combined-trim.fq"
				]))

				os.system(' '.join([
					"cat",
					se,
					"./TRIM/se.trim.fq",
					">>",
					"./FOO.fq"
				]))

				os.system(' '.join([
					"mv",
					"./FOO.fq",
					"./TRIM/se.trim.fq"
				]))

				os.system(' '.join([
					"rm",
					"./FOO.trim.fq"
				]))

		## Zip source files
		Zip(file1)
		if file2 != "interleaved":
			Zip(file2)

		## Merge Paired Reads using FLASH if user specified
		if FLASH:
			se, pe = flash("./TRIM/combined-trim.fq")

			## Merge orphaned single read files with newly extended single reads from FLASH
			os.system(' '.join([
				"cat",
				se,
				"./TRIM/se.trim.fq",
				">>",
				"./"+OUTPUT+"/"+sampleID+".se.qc.fl.fq"
			]))

			os.system(' '.join([
				"mv",
				pe,
				"./"+OUTPUT+"/"+sampleID+".pe.qc.fl.fq"
			]))

			os.system(' '.join([
				"pigz -p",
				PROCESSORS,
				"./"+OUTPUT+"/"+sampleID+".pe.qc.fl.fq",
				"./"+OUTPUT+"/"+sampleID+".se.qc.fl.fq"
			]))

			## Save output file names
			OUTPUT_DICT[sampleID] = ["./"+OUTPUT+"/"+sampleID+".pe.qc.fl.fq.gz","./"+OUTPUT+"/"+sampleID+".se.qc.fl.fq.gz"]

		else:
			## Compress Files and Move to Output Folder
			os.system(' '.join([
				"mv",
				"./TRIM/combined-trim.fq",
				"./"+OUTPUT+"/"+sampleID+".pe.qc.fq"
			]))

			os.system(' '.join([
				"mv",
				"./TRIM/se.trim.fq",
				"./"+OUTPUT+"/"+sampleID+".se.qc.fq"
			]))

			os.system(' '.join([
				"pigz -p",
				PROCESSORS,
				"./"+OUTPUT+"/"+sampleID+".pe.qc.fq",
				"./"+OUTPUT+"/"+sampleID+".se.qc.fq"
			]))

			## Save output file names
			OUTPUT_DICT[sampleID] = ["./"+OUTPUT+"/"+sampleID+".pe.qc.fq.gz","./"+OUTPUT+"/"+sampleID+".se.qc.fq.gz"]

		## Remove Temp Folder "TRIM" and clean up
		os.system(' '.join([
			"rm -fr TRIM",
			"out.extendedFrags.fastq",
			"*pe",
			"*.se"
		]))

	if LOG == 'ON':
		log.write("\n-- Quality Trimming --\n Following QC your File Names Are Now: \n"+str(OUTPUT_DICT))


###########
# ASSEMBLY
###########

if re.search("2", PROCESSES):
	print "\n\n-- Performing Assembly --"

	# Check for existence of OUTPUT_DICT
	try:
		ASSEMBLE_DICT = OUTPUT_DICT
		
	except NameError:
		ASSEMBLE_DICT = FILE_DICT

	for sampleID, value in ASSEMBLE_DICT.iteritems():
		# Unzip Files
		pe = UnZip(value[0])
		se = UnZip(value[1])

		## IDBA
		if re.search("1|IDBA", ASSEMBLER):
			print "\n\n-----Using IDBA_UD-----\n"
       	                os.makedirs('./' + OUTPUT + "/ASSEMBLY/IDBA/" + sampleID)

			##If necessary convert from .fastq to .fasta
			pe_fa = Convert_FQ(pe)

			## Run IDBA
			idba(pe_fa, sampleID)

			## Zip FASTA
			Zip(pe_fa)

			print "\n\n--- Assembly Completed --- \n\nYou've Completed Assembly with IDBA.The files you are most interested in are entitled \"contig.fa\".\n" 

		## RAYmeta
		if re.search("2|RAYmeta", ASSEMBLER):
			print "\n\n-----Using RAYmeta-----\n"

			# Select KMER length
			KMER = '39'
	
			# Convert from .fastq to .fasta if necessary
			pe_fa = Convert_FQ(pe)
			se_fa = Convert_FQ(se)

			## Run Ray-meta
			ray(pe, se, sampleID)

			## Zip FASTA
			Zip(pe_fa)
			Zip(se_fa)

			print "\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with RAY-meta. Be sure to look at the RAY-meta documentation for other tools in the package for community ecology of metagenomic data.\n" 
			print "\nThe files you are most interested in are entitled \"Contigs.fasta\" and \"Scaffolds.fasta\"."
			print "\nThe assembly statistics can be found in the \"OutputNumbers.txt\"."

		## MEGAHIT
		if re.search("3|MEGAHIT", ASSEMBLER):
			print "\n\n-----Using MEGAHIT-----\n"
       	                os.makedirs('./' + OUTPUT + "/ASSEMBLY/MEGAHIT/" + sampleID)

			# Run MEGAHIT  (super simple b/c it is flexible with file formats)
			megahit(pe, se, sampleID)

		## metaSPAdes
		if re.search("4|metaSPAdes", ASSEMBLER):
			print "\n\n-----Using metaSPAdes-----\n"
       	                os.makedirs('./' + OUTPUT + "/ASSEMBLY/METASPADES/" + sampleID)

			# Run metaSPAdes 
			metaSPAdes(pe, se, sampleID)


		## Zip FASTQ
		Zip(pe)
		Zip(se)

log.close()
command.close()
end = timeit.default_timer()

print end - now
