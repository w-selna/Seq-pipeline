#!/usr/bin/env python3

import os
from os import path
import shutil
import re
import rpy2
import rpy2.robjects as robjects
import sys

#---------Static Variables---------#

#all paths to applications, will build paths to files later.
prin = "/home/linh/Documents/Pipeline_files/prinseq-lite" 
fastqc = "/home/linh/Documents/Pipeline_files/FastQC"
cutadapt = "/home/linh/Documents/Pipeline_files/cutadapt-1.2.1"
bwa = "/home/linh/Documents/Pipeline_files/bwa.kit"
pepr = "/home/linh/Documents/Pipeline_files/PePr-master/PePr"

geno = "/media/linh/OS/ref_genomes/genomes"
# list of applications
apps = ["cutadapt ", "perl prinseq-lite.pl", "bwa index ","bwa aln", "bwa samse", "bwa sampe", "bedtools bamtobed", "python3 ~/Documents/Pipeline_files/PePr-master/PePr/PePr.py"]

# lists of all options that will be used. for make_cmd function
fastqc_ops = [] 							# useless since only running with files, use file names, or not make_cmd function
prin_se_ops = ["-fastq", "-out_good", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-derep", "-derep_min", "-lc_method", "-lc_threshold", "-trim_to_len","-graph_data", "-noniupac","-phred64"]  							#last 3 dont take values 
prin_se_ops_non64 = ["-fastq", "-out_good", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-derep", "-derep_min", "-lc_method", "-lc_threshold", "-trim_to_len","-graph_data", "-noniupac"]
prin_pe_ops = ["-fastq", "-fastq2", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-derep", "-derep_min", "-lc_method", "-lc_threshold", "-trim_to_len","-graph_data", "-noniupac", "-phred64"]  						#last 3 dont take values,
prin_pe_ops_non64 = ["-fastq", "-fastq2", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-derep", "-derep_min", "-lc_method", "-lc_threshold", "-trim_to_len","-graph_data", "-noniupac"]
cut_ops =  ["-e", "-a", "-g"]							# first one may have multiple values, cant use make_cmd
bwa_aln_ops = []
bwa_samse_ops = []
bwa_sampe_ops = []
bedtools_ops = ["-i", ">"]
pepr1_ops = ["-c","-i","-n","-f","--peaktype=", "2>"]  			# run pepr twice, once for narrow peaks, once for broad peaks.

#quasi static variables
postcut_read_len = []


#----------Functions!----------#


def out_name_fastq(pth):	# takes in a path to a fastq file and returns an appropriate output file name, designed with FastQC in mind, works there, maybe not every where.
				# pth: path to a file
	base = os.path.basename(pth)
	base = base.replace(".", "_")					# do this specifically because of hof FASTQC outputs files
	base = base + "_out"						# add a designation that this in an output file
	return(base)

def isPE(pth):			# tests if a fastq file is PE of not based on read ID's containing a '#' or not. '#' indicated PE, this is how Illumina does it at least.
				# pth: path to a file
	f_open = open(pth)
	line = f_open.readline()					# read the first line, is a descriptive one about the first read in file
	if "#" in line:							# testing for '#'
		f_open.close()
		return (True)
	else:
		f_open.close()
		return (False)
	f_open.close()



def make_filepath_list(path): 	#makes a list of paths to files kept in a folder, used for fastqc
				#path: a path to the folder with files in it. 
	files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]	# Extracts path of a file in 'path' folder
	sfiles = sorted(files)											# sorts list, why alphabetical ordering of files in set up is so important.
	return sfiles

def get_trim_len(pth): 		# Trims reads that were not trimmed down below 90% of the longest read in that file, we trim them to 90% of longest because 
				# all the other trimming creates and unnatural GC bias in the tail region, and this is eliminating that. 
				# pth: a path to the txt file created by FASTQC. 
	f_open = open(pth)
	line = f_open.readline()
	while (re.search(r'Sequence length',line) == None):		# finds Sequence length section of FastQC .txt file report
		line = f_open.readline()
	Seq,length,ex =  line.split("\t")
	if ("-" in length):						# checks if there is a range of length, in txt file a range of lengths is reported as  # - #.
		l = list(length)
		for i in l:
			if i != "-":					# selects larger of the 2 values by "eating up" the string til the - has passed
				l = l[1:]
			else:
				l = l[1:]
				tot = int("".join(l))			# makes the number and integer
				trim = int(tot- tot*.1)			# maths to get what 90% of total is, off by 1 sometimes due to interger calculations
				return (trim)		
	else:
		trim = int(int(length)-int(length)*.1)			# maths to get what 90% of total is, off by 1 sometimes due to interger calculations
		return (trim)
	f_open.close()



def get_trim_qual(pth): # returns an integer corresponding to the phred score of the last peak in the "Per sequence quality scores" graph.
			# an important behind the scenes selection process. Simply, this method was arrived at by testing different ways of trimming down reads to ensure high qaulity. other methods tested, 1. Trimming all reads to the length of when the average phred score for that position dropped below 25. ended up only slightly improving average read quality, while also trimming good bases too, 2. Trimming to a specific quality, by manually decided. This worked much better, but produces reads of variable length and required a manual input, the automation of this input was the motivation to make this function. I found that the sequence quality scores usually show a nice increace up to some peak around 35-39 phred, as this was the mode, then it would be easy to get most reads to be of that quality or higher. the algo to find this peak is a backward search from phred 40 down to 0, checking if the n+1 phred number has a higher or lower number of reads with that average phred score than the nth position. If the n+1st is higher than we are past the highest local mode and we accept n+1st phred value as teh quaility to trim to. 
			# pth is a path to the actual txt file.
	f_open = open(pth)
	line = f_open.readline()
	master = []
	while (re.search(r'>>Per sequence quality scores',line) == None):   # iterates til get to section
		line = f_open.readline()
	line = f_open.readline()
	line = f_open.readline()
	while (re.search(r'>>END_MODULE',line) == None):		    # fills array
		line = line.strip()		
		lst =  line.split("\t")
		master.append(lst)
		line = f_open.readline()
	post = float(master[len(master)-1][1])				    # reverse order checking to find last max point.
	for c in range(1,(len(master)-1)):
		pre = float(master[len(master)-1-c][1])
		if (pre<= post):
			return(float(master[len(master)-c][0]))		   # we use this number because this number indicates a score most sequences, and therefore bases have that is still of high 										   # quality. Cutting above this would cut too little, and cutting below this ok, but just harder to determine, computationally, 										   # when to stop and assign the cut value. 
		else:
			post = pre
	return(20)							   # safeguard number
	f_open.close()



def get_overrep_seq (pth): 	# a function to get any adaptors that may be present in the data set, useful especially if the user doesnt know the adaptors used. 
				# pth is a path to the FastQC data txt file.
	adap_list = []
	f_open = open(pth)
	line = f_open.readline()
	while (re.search(r'>>Overrepresented sequences',line) == None): 	# finds start of Overrepresented sequences section
		line = f_open.readline()
	line = f_open.readline()
	if (re.search(r'>>END_MODULE',line) != None):				# exits if no sequences
		return (adap_list)
	line = f_open.readline()
	while (re.search(r'>>END_MODULE',line) == None):
		Sequence,Count,Percentage,Possible_Source = line.split("\t") 	# finds adaptors reported by FastQC
		adap_list.append(Sequence)					# adds adaptor to list
		line = f_open.readline()	
	return (adap_list)
	f_open.close()

def make_cmd(app, options, parameters): #command line command generator for easy to make commands
					# app is the application you want to use i.e prin-seq, ect.
					# options is a predefined list of options like prin_pe_ops.
					# parameters is a parallel list to optinons of parameter values for that app. generated in the code chunks for that application.         
	cmd = app
	for i in range(len(options)):
		if (len(parameters) > i):					# 'parameters' is always the same length or shorter than options, if there are still paramnters to put in, then we do so
			cmd = cmd + ' ' + options[i] + ' ' + parameters[i]
		else:								# some options do not take parameters, and are added onto the end.
			cmd = cmd + ' ' + options[i]
	return(cmd)

def qualencode(pth):	# A function to extract the quality encoding of the reads. this function may break if not using illumina reads
			# pth is a path to the FastQC data txt file.
	f_open = open(pth)
	line = f_open.readline()
	while (re.search(r'Encoding',line) == None): 				# finds Encoding line
		line = f_open.readline()
	title,num,end = line.split("\t")
	try:
		val = float(num[-3:]) 						# try to get a number from the line, done in case of mess ups. This try should be able to handel it.
	except ValueError:
		val = 1.2 							# default to not phred64
	if ((val >= 1.3) & (val<= 1.7)):					# if worked normally then test in what range the number falls, inside range here
		return(1)
	elif ((val < 1.3) | (val > 1.7)):					# outside range here
		return(-1)
	else:
		return(-1)							# default to not phred64

#tested, may not work only because of size of genome, testing on a new microbial genome, or again on ecoli would be a good idea.
def make_refgenome(genome,txt,bwa_pth):	# a function to make new reference genomes and add their names to a txt file so they may be used at a later date as well.
					# genome is a path to a fasta file, the name of the file should be descriptive of the organism and will be used in the text file to denote what reference genomes 
					# 	are avaliable (I am assuming only bacterial genomes right now so it wont be too bad to store all of them locally.)
					# txt is the text file genomes in ref_genomes folder that lists the genomes in there. variable "geno"
					# bwa_pth is the path to bwa, will use the index function from there to make the appropriate index of the fasta file.

	basename = os.path.basename(genome)
	print("Your reference fasta file is being moved to the folder /media/linh/OS/ref_genomes")		
	os.rename(genome, "/media/linh/OS/ref_genomes/"+basename)				#moves to ref genome folder
	genome = "/media/linh/OS/ref_genomes/"+basename
	os.chdir(bwa_pth)
	os.system(apps[2] + " /media/linh/OS/ref_genomes/"+basename)				#index command
	list_txt = open(txt,'a')
	list_txt.write(basename)								#adds its name to the text file	
	list_txt.close()
	return(genome)

sam_to_bam = robjects.r('''
	sam_to_bam <- function(lst_o_pths, lst_o_outs){				# A function to convert SAM file types to BAM file types with Rsamtools in R. 
	  k = 1									# lst_o_pths: a list containing paths to the SAM files to convert
	  suppressMessages(library(Rsamtools));					# lst_o_outs: a list of paths that are locations to store the BAM files. 
	  x = vector("list", length(lst_o_pths) )
	  for (i in 1:length(lst_o_pths))
	  {
	    x[[k]]= asBam(lst_o_pths[[i]],lst_o_outs[[i]])			# do the conversion
	    k  = k + 1
	  }
	  x									#return the vectir of BAMs
	}
''')

########### V work to do here!! V ###########    not saving pdf correctly or making 
htseqtools = robjects.r('''
	htseqtools <- function(lst_o_pths, lst_o_fldrs,CI){			# a function to execute htSeqTools functions on the the bed files it is given. creates plots and checks for
		suppressMessages(library(rtracklayer))				# immunopercipitation success.
		suppressMessages(library(htSeqTools))				# lst_o_pths: a list containing paths to the BED files to analyze
		my.rangeddatalst = vector("list", length(lst_o_pths))		# lst_o_fldrs: a list containing paths to folders, these will be the locations to store plots made in the function
		for (i in 1:length(lst_o_pths))					# CI: a list of integers, positive and negative, that is parallel to lst_o_pths that denote ChIP vs Control and number 
		{								# in that specific section.
			my.rangeddatalst[[i]] <- RangedData(import.bed(con=lst_o_pths[[i]]))			# create a ranged data list of the bed files
		}
		x = RangedDataList(my.rangeddatalst)
		nme = vector("list", length = length(lst_o_pths))						# assign names to each list element
		ch_plce = 1
		co_plce = 1
		for(i in 1:length(lst_o_pths))									#Generate the names
		{
			if(CI[[i]] == -1)
			{			
				nme[[i]] = paste("Control_rep",co_plce,sep = "")				# control names with increasing number designation
				co_plce = co_plce +1
			}else{
				nme[[i]] = paste("IP_rep",ch_plce,sep = "")
				ch_plce = ch_plce +1
			}
		}
		names(x) = nme

		pdf(paste(lst_o_fldrs[[1]],'/','mds_plots.pdf', sep = ''))					# create an MDS plot of the data, right idea, but polt commands dont work...
		#file.create(paste(lst_o_fldrs[[1]],'/','mds_plots.pdf', sep = ''))
		MDS = cmds(x,k=2)
		plot(MDS)
		#dev.copy(pdf, paste(lst_o_fldrs[[1]],'/','mds_plots.pdf', sep = ''))
		dev.off()

		ssd = ssdCoverage(x)										# Calculate the Standard deviation in coverage of the different files
		write.csv(ssd, paste(lst_o_fldrs[[1]],'/ssd.csv',sep = ''), row.names = lst_o_pths)		# write values to a CSV file, but npot working as it did in Rstudio

		for (i in 1:length(x))										# loop to graph lorenze curves and calculate gini coefficients of files
		{	
			pdf(paste(lst_o_fldrs[[i]],'/','gini.pdf', sep = ''))					# again graphing os off, ikd why, need to test it in a smaller program. 
			#file.create(paste(lst_o_fldrs[[i]],'/','gini.pdf', sep = ''))
			giniCoverage(x[[i]],mc.cores=1,mk.plot=TRUE)		
			#dev.copy(pdf, paste(lst_o_fldrs[[i]],'/','gini.pdf', sep = ''))
			dev.off()
		}
	}
''')

    #----------Code Chunks---------#
#----------Perliminary Questions----------#

geno_open = open(geno)     							# get the reference genomes already on file.
ref_list = []
new = -1
for lin in geno_open:
	ref_list.append(lin)
geno_open.close()

print("Hello,\n\nWelcome to the Chip-Seq/ChIP-exo Pipeline writen by Weston Selna\n\nThis is for use by the IT Lab, any and all reference documents for this Pipeline can be found at /home/linh/Documents/Pipeline_files/pipeline_docs.pdf") # intro
pth1 = input('Please enter the path to the folder containing the .fastq files that you want analyzed\n\nPlease make sure you have reviwed and kept to the file name requierments\n\nFailure to do so will most likely result in an error and this pipeline exiting, or it will somehow work but give you meaningless results.\n\n') # most important location! The Main folder.

skip = input("skip ahead?(y/n)")						# for debugging and testing purposes, will be removed.
if (re.search(r"[y,Y]",skip) == None): 
	print("The avaliable reference genomes to map to are:")			# this section is for determination of the organism under study.
	for acc in ref_list:
		print ("%s" % acc)

	s_new = input("Will you be using one of these? (y/n)\n")
	if (re.search(r"[y,Y]",s_new)!=None):
		new = 0
		accession = input("Which then? (Please type the name in exactly as seen in the prior print out)\n")  
		pthgl = [r for r in os.listdir("/media/linh/OS/ref_genomes") if (r == accession)]
		pthg = os.path.join("/media/linh/OS/ref_genomes",pthgl[0])
	else:
		pthg = input('Enter the path to the file containing the .fasta files you would like to set as a genome.\nThis will be saved for later use so please make sure the file in named appropriately so other may use it too.\n')					# a new organism to align too.
		new = 1

	if (new == 1):
		print("Making new index of reference genome.")	
		pthg = make_refgenome(pthg,geno,bwa)				# creation of a new indexed genome

	#----------Part 1: Preprocessing----------#
	fplist = make_filepath_list(pth1)					#paths to fastq files
	# checking that each path has word chip or control in it and that there are at least 4 files present. Exits if conditions are not met
	if (len(fplist)<4):
		sys.exit("There are not enough files. The minimum number of files allowed is 4, there should be at least 2 replicates for both the case and control states, and correspondingly at least 1 single end read fastq file each totaling 4. More files are allowed above this for extra replicates and PE data") 
	for pths in fplist:
		if((re.search(r'control',pths) == None)&(re.search(r'chip',pths) == None)):
			sys.exit("At least one file is named incorrectly, please review the required naming schema. All files should contain the word 'control' or 'chip' in them to designate their group.")
	print("beginning FastQC analysis")

	cmd1 = 'fastqc'								# builds a factqc command.
	for f in fplist:
		cmd1 = cmd1 + " " + f

	os.chdir(fastqc)							# changes to teh appropriate directory
	os.system(cmd1)  							# executes command
		
	print("Beginning automated read analysis: filtering and trimming")
	txtplist = []
	files = [f for f in os.listdir(pth1) if os.path.isfile(os.path.join(pth1, f))]
	for a in files:  							# getting/making paths to appropriate fastqc_data.txt files in different folders 
		aa = a.replace(".","_")
		aa_p = aa.split()
		aa_p.append("c")						# ths for loop is very specific to exactly wht FastQC outputs.
		aaa = ''.join(aa_p)
		if os.path.isdir(os.path.join(pth1, aaa)):		
			txtplist.append(os.path.join(os.path.join(pth1, aaa),"fastqc_data.txt"))

	txtplist = sorted(txtplist)      					# to maintain the alphabetical order
	PE_ToF = []
	trim_qual_params = []		
	trim_len_params = []							# these are important factors from the FastQC report and fastq file. 
	adapter_params = []
	qual_encoding = []

	for b in txtplist:							# building the parameter lists 
		trim_qual_params.append(int(get_trim_qual(b)))    
		trim_len_params.append(int(get_trim_len(b)))
		postcut_read_len.append(int(get_trim_len(b)))
		adapter_params.append(get_overrep_seq(b))
		qual_encoding.append(qualencode(b))
	for f in fplist:
		PE_ToF.append(isPE(f))

	switch = 0
	save = None								#these are place saver variables for the main preprocessing 'for' loop.
	save_name = None

	for f in fplist:
		name_out = out_name_fastq(f) 
		#cutadapt
		if (adapter_params[fplist.index(f)] != []):	
			cut_params = []
			cut_params.append(adapter_params[fplist.index(f)])
			cut_params.append(.1)							# error allowed in adaptor match. default is .1
			cut_params.append(f+" > "+os.path.join(pth1,name_out+"_cut"+".fastq")) 	# output names of files,
			cut_params.append('2>'+ os.path.join(pth1,name_out+"_cutadapt"+".log"))	# suppressing standard error output from cutadapt. send it to a log file.
			cmd2 = apps[0]
			for e in cut_params:
				if (type(e) == list):
					for g in e:
						cmd2 = cmd2 + ' ' + cut_ops[cut_params.index(e)] + ' ' + g + '$' + ' ' + cut_ops[cut_params.index(e)+1] + ' ' + '^' + g # adds all adapters to command.
				elif(type(e)!= float):										# if not -e option
					cmd2 = cmd2 + ' ' + e
				else:												#for -e option
					cmd2 = cmd2 + ' ' + cut_ops[cut_params.index(e)] + ' ' + str(e)
			os.chdir(cutadapt)
			os.system(cmd2)
			fplist[fplist.index(f)] = os.path.join(pth1,name_out+"_cut"+".fastq")
			f = (os.path.join(pth1,name_out+"_cut"+".fastq"))  #new path to new adapterless reads	

		if (PE_ToF[fplist.index(f)] == True):
		#start of prin-seq
			if (switch == 1):								# enter here is on the second of the 2 PE fastq files	
				prin_params = []
				prin_params.extend((f,save,"null","ld,gc,qd,ns,pt,ts,aq,de,sc,dn",'20',str(trim_qual_params[fplist.index(f)]),'0','1','12',"entropy",'70',str(trim_len_params[fplist.index(f)])))
				if (qual_encoding[fplist.index(f)] == 1):				#different options depending on encoding of quality
					cmd3 = make_cmd(apps[1],prin_pe_ops,prin_params)
				else:
					cmd3 = make_cmd(apps[1],prin_pe_ops_non64,prin_params)
				cmd3 = cmd3 + " 2> " + os.path.join(pth1,name_out+"_prin"+".log")	#suppressing output to a log file
				os.chdir(prin)
				os.system(cmd3)

				# finding the fastq files in  pth1 folder, moving, renaming and getting their new paths.
				outs = [f for f in os.listdir(pth1) if (re.search(r"prinseq_good",os.path.basename(f)) != None)]
				outs[0] = outs[0].split("_")
				if (re.search(outs[0][0],f)):
					outs[0] = "_".join(outs[0])
					outs[0],outs[1] = outs[1],outs[0]				# this if-else maintains order of files.
				else:
					outs[0] = "_".join(outs[0])
				outs_p = [None,None]
				for g in outs:
					if (outs.index(g) == 1):
						outs_p[outs.index(g)] = os.path.join(pth1,g)
						os.rename(outs_p[outs.index(g)],os.path.join(pth1,name_out+"_prin_PE.fastq"))		#the paths to the new fastq files. 
						fplist[fplist.index(f)] = os.path.join(pth1,name_out+"_prin_PE.fastq")
					else:
						outs_p[outs.index(g)] = os.path.join(pth1,g)
						os.rename(outs_p[outs.index(g)],os.path.join(pth1,save_name+"_prin_PE.fastq"))
						fplist[fplist.index(f)-1] = os.path.join(pth1,save_name+"_prin_PE.fastq")
				switch = 0
				save = None
				save_name = None
			else:
				switch = 1								# encountered the first of the pair of PE fastq files.
				save = f								# save associated files
				save_name = name_out
				print("PE pair analysis.")
		else:											# for SE data
			prin_params = []
			prin_params.extend((f,os.path.join(pth1,name_out+"_prin"),"null","ld,gc,qd,ns,pt,ts,aq,de,sc,dn",'20',str(trim_qual_params[fplist.index(f)]),'0','1','12',"entropy",'70',str(trim_len_params[fplist.index(f)])))
			if (qual_encoding[fplist.index(f)] == 1):					#different options depending on encoding of quality
				cmd3 = make_cmd(apps[1],prin_se_ops,prin_params)
			else:
				cmd3 = make_cmd(apps[1],prin_se_ops_non64,prin_params)
			cmd3 = cmd3 + " 2> " + os.path.join(pth1,name_out+"_prin"+".log")		#suppressing output to a log file
			os.chdir(prin)
			os.system(cmd3)
			fplist[fplist.index(f)] = os.path.join(pth1,name_out+"_prin.fastq")
	
	print("End of Preprocessing\nYou can find all files made in this process in %s, each with descriptive filenames" %pth1)

	#----------Part 2: Alignment----------#

	print("Starting Alignment process")
	samplist = []
	print("using BWA-ALN, best for reads less than 70bp in length\nThe longest reads for the given files after trimming are as follows:")
	print(postcut_read_len)
	for f in fplist:										# display informative information about the reads, post cutting. relevant to alignment process.
		print ("for file %s ----------- %d" %(os.path.basename(f), postcut_read_len[fplist.index(f)]))


	both = 0
	save_sai = None			# place holder variables for alignment for loop.
	save_reads = None
	for f in fplist:
		# pe vs se selector/matcher
		if (PE_ToF[fplist.index(f)] == False):
			#SE
			if (qual_encoding[fplist.index(f)] == 1): 								# different parameters based on quality encoding 
				bwa_params = [apps[3],"-I",pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")] 	#-14 because descriptive names wer getting too long.
			else:
				bwa_params = [apps[3],pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
			cmd4 = " ".join(bwa_params)										# bwa aln (-I) ref.fa reads.fq > reads.sai
			os.chdir(bwa)
			os.system(cmd4)
			sesai = os.path.join(pth1,os.path.basename(f)[:-14]+".sai")
			bwa_params_samse = [apps[4],pthg,sesai,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_se.sam")]
			cmd5 = " ".join(bwa_params_samse)									# bwa samse ref.fa reads.sai reads.fq > aln-se.sam
			os.chdir(bwa)
			os.system(cmd5)
			samplist.append(os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_se.sam"))  # will not make a perfectly parallel list cause and PE data will only out 1 sam file for 2 fastq files.
		else:	#PE
			if (both == 1):
				if (qual_encoding[fplist.index(f)] == 1):								# different parameters based on quality encoding
					bwa_params_2 = [apps[3],"-I",pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				else:
					bwa_params_2 = [apps[3],pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				cmd7 = " ".join(bwa_params_2)										# bwa aln (-I) ref.fa reads.fq > reads.sai
				os.chdir(bwa)
				os.system(cmd7)
				pesai = os.path.join(pth1,os.path.basename(f)[:-14]+".sai")
				both = 0
				bwa_params_sampe = [apps[5],pthg,save_sai,pesai,save_reads,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_pe.sam")]
				cmd8 = " ".join(bwa_params_sampe)									# bwa sampe ref.fa reads.sai reads.fq > aln-pe.sam
				os.chdir(bwa)
				os.system(cmd8)
				samplist.append(os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_pe.sam"))# will not make a perfectly parallel list cause and PE data will only out 1 sam file for 2 fastq files. the sam file will have the name of the second pair that goes into this fcn.
			else:
				if (qual_encoding[fplist.index(f)] == 1):								# different parameters based on quality encoding
					bwa_params_1 = [apps[3],"-I",pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				else:
					bwa_params_1 = [apps[3],pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				cmd6 = " ".join(bwa_params_1)										# bwa aln (-I) ref.fa reads.fq > reads.sai
				os.chdir(bwa)
				os.system(cmd6)
				save_sai = os.path.join(pth1,os.path.basename(f)[:-14]+".sai")
				save_reads = f
				both = 1

#--------Part 3: Post QC--------# 
print("Starting post alignment, pre peak calling quality control")
print("Converting to other file types for analysis") 

os.chdir(pth1)
samoutnames = []
for sam in samplist:									# make sam out names.
	bse = os.path.basename(sam)
	if (re.search(r"sorted",bse) != None):     					# difference between PE and SE.
		b_name = bse[:-17] 							# truncate names to make them shorter, not losing any of the personal descriptive name tho
	else:
		b_name = bse[:-10]
	samoutnames.append(os.path.join(pth1,b_name+"asbam")) #(revisit documentation of sam_to_bam functions in R if needed)    
 
bamplist_0 = sam_to_bam(samplist,samoutnames) 						#make bams, are R objects here #not the best implementation since dont use the bam file for anything.
bamplist = []
for i in range(0, len(samoutnames)):
	bamplist.append(bamplist_0[i][0]) 						#get paths to bams, as strings here (not r objects)

bedplist = []
bamoutnames = []
for bm in bamplist:									#convert to bed for htSEQtools. makecmd for bedtools
	bm_out = os.path.basename(bm)							# generation of the names
	bm_out = bm_out[:-9]+"as_bed.bed"
	bm_out = os.path.join(pth1,bm_out)
	bedplist.append(bm_out)
	bedtool_params = [bm,bm_out]
	os.system(make_cmd(apps[6],bedtools_ops,bedtool_params))  			#makes beds

print("Corresponding BED files have been made from the alignments. These .bed files can be viewed in IGV or on the UCSC genome browser. See the additional help text for how to use IGV to locally view the alignments\n")


#htSEQtools
rplotfolderlist = []
for flder in bedplist:   						# makes folders for R plots, will be named so differentiation between samples is easy. 
	fldr_nme = os.path.basename(flder)[:-18]			# again reduce names to make managable
	os.mkdir(pth1+"/Rplots_"+fldr_nme)				# makes folder
	rplotfolderlist.append(pth1+"/Rplots_"+fldr_nme)
	

CoI = []								# Holder variables
ctrl_pths = []
chip_pths = []

for bd in bedplist:							#seprate beds into ChIP and ctrl
	if(re.search(r'control',os.path.basename(bd)) != None):		#control and chip are required portions of file names. looks for that and assigns a 1 or -1 and distributes them into seprate lists. 										#Parallelness may be lost, but still recoverable if needed. should still be parallel wrt PE being next to each other and what not.
		CoI.append(-1)
		ctrl_pths.append(bd)
	elif(re.search(r'chip',os.path.basename(bd)) != None):
		CoI.append(1)
		chip_pths.append(bd)
	else:
		print("ERROR......\nABORT.....\nFILES NAMED INCORRECTLY")# Double check correct naming
		exit


print("Creating MDS plots of coverage\nCalculating Gini index for all samples to see if immunopercipitation worked\nCalculating the coverage standard deviaton\n")

htseqtools(bedplist,rplotfolderlist,CoI)				#Statistical htSeqTools analysis

print("Coverage and immunopercipitation statistics and plots are avaliable to view in the newly crated XXXX_Rplots folder in the main folder holding the original fastq files.")

#####-------Part 4: Peak Calling---------#####
#PePr
#build PePr params
pepr_params_sharp = [','.join(chip_pths), ','.join(ctrl_pths), 'sharppeaks', 'bed', 'sharp', "pepr_sharp.log"]# We run this twice because we are not sure what we are looking for, may be a good question to ask in 
pepr_params_broad = [','.join(chip_pths), ','.join(ctrl_pths), 'broadpeaks', 'bed', 'broad', "pepr_broad.log"]# beginning, but it is quick so isnt really a big deal, just hard drive space is the problem. 


cmd9 = make_cmd(apps[7],pepr1_ops,pepr_params_sharp)
cmd10 = make_cmd(apps[7],pepr1_ops,pepr_params_broad)
# must remove space after =    
cmd9 = cmd9.replace("= ","=")
cmd19 = cmd10.replace("= ","=")

os.chdir(pth1)
print('Starting Peak Calling with Peak Prioritization softwear.')
os.system(cmd9)	 
os.system(cmd19)
print("Peak calling done. You will be able to find multiple output files for this process.\n")

print("Thank You for using the IT lab ChIP-seq/ChIP-exo Pipeline!\nYour files should all be located in the original folder you indicated at the beginning of the process.\n\nThis pipeline has created a wealth of useful files for you to review at your lesiure.\n")
	
