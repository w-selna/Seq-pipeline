#!/usr/bin/env python3
#####This is the Exo version of the pipline

#Differences from Seq version:
# 1: does not remove seq copies in prin-seq since them seem to make up most of the Peaks in exo,
# 2: MACE instead of PePr. 
# 3: dealing with non-uniquely mapped reads? (this should be done for both)
# file alphabetical order should still be in there, but removed the chip control differentiation since PePr only works on sets of replicates. 

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
geno = "/media/linh/OS/ref_genomes/genomes"
mace = "/home/linh/Documents/Pipeline_files/MACE-1.2/bin"

# list of applications
apps = ["cutadapt ", "perl prinseq-lite.pl", "bwa index","bwa aln", "bwa samse", "bwa sampe", "bedtools bamtobed", "samtools idxstats", "cut", "python /home/linh/Documents/Pipeline_files/MACE-1.2/bin/preprocessor.py","/home/linh/Documents/Pipeline_files/wigToBigWig","python /home/linh/Documents/Pipeline_files/MACE-1.2/bin/mace.py"]

# lists of all options that will be used. for make_cmd function
fastqc_ops = [] # useless since only running with files, use file names, or not make_cmd function
prin_se_ops = ["-fastq", "-out_good", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-lc_method", "-lc_threshold","-graph_data", "-noniupac","-phred64"]  #last 3 dont need values 
prin_se_ops_non64 = ["-fastq", "-out_good", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-lc_method", "-lc_threshold", "-graph_data", "-noniupac"]	     #last 2 dont need values 
prin_pe_ops = ["-fastq", "-fastq2", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-lc_method", "-lc_threshold", "-graph_data", "-noniupac", "-phred64"]  #last 3 dont need values,
prin_pe_ops_non64 = ["-fastq", "-fastq2", "-out_bad", "-graph_stats", "-min_len", "-trim_qual_right", "-ns_max_n", "-lc_method", "-lc_threshold", "-graph_data", "-noniupac"] 	     #last 2 dont need values 
cut_ops = ["-e", "-a", "-g"] # first one may have multiple values, can't use make_cmd
bwa_aln_ops = []
bwa_samse_ops = []
bwa_sampe_ops = []
bedtools_ops = ["-i", ">"]
#smtools_idx_ops = []
ct_ops = ["-f",">"]
prepros_ops = ["-i", "-r", "-o", "-w", "2>"]
mace_ops = ["-s", "-f", "-r", "-o","-p","2>"]	

#quasi static variables
postcut_read_len = []
seq_len = 0

#----------Functions!----------#


def out_name_fastq(pth):# takes in a path the a fastq file and returns an appropriate output file name, designed with FastQC in mind, works there and maybe not every where.
	base = os.path.basename(pth)
	base = base.replace(".", "_")
	base = base + "_out"
	return(base)

# tests if a fastq file is PE of not based on read ID's containing a '#' or not. '#' indicated PE, this is how Illumina does it at least.
def isPE(pth): 
	f_open = open(pth)
	line = f_open.readline()
	if "#" in line:
		f_open.close()
		return (True)
	else:
		f_open.close()
		return (False)
	f_open.close()


#makes a list of paths to files kept in a folder, used for fastqc (so far) messes up order of shit, not parallel to what is in folder. fix?    #### has been fixed with sorting.
def make_filepath_list(path): #path: a path to the folder with files in it. 
	files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
	sfiles = sorted(files)
	return sfiles

def get_trim_len(pth): # pth is a path to the actual txt file. for reads that were not trimmed down below 90% of the longest read in that file, we trim them to 90% of longest because all the other trimming creates and unnatural GC bias in the tail region, and this is eliminating that. 
	f_open = open(pth)
	line = f_open.readline()
	while (re.search(r'Sequence length',line) == None):
		line = f_open.readline()
	Seq,length,ex =  line.split("\t")
	if ("-" in length):
		l = list(length)
		for i in l:
			if i != "-":
				l = l[1:]
			else:
				l = l[1:]
				tot = int("".join(l))
				trim = int(tot- tot*.1)
				return (trim)		
	else:
		trim = int(int(length)-int(length)*.1)
		return (trim)
	f_open.close()


#returns an integer corresponding to the phred score of the last peak in the "Per sequence quality scores" graph.
def get_trim_qual(pth): # pth is a path to the actual txt file, 
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
			return(float(master[len(master)-c][0]))		   # we use this number because this number indicates a score more sequences, and therefore bases have that is still of high quality. Cutting above this wouls cutt too much, and cutting below this ok, but just harder to determine, computationally what the best trim quality score is. 
		else:
			post = pre
	return(20)							   # safeguard number
	f_open.close()


# a function to get any adaptors that may be present in the data set, useful especially if the user doesnt know the adaptors used. 
def get_overrep_seq (pth): # pth is a path to the actual txt file,
	adap_list = []
	f_open = open(pth)
	line = f_open.readline()
	while (re.search(r'>>Overrepresented sequences',line) == None): # finds start of Overrepresented sequences section
		line = f_open.readline()
	line = f_open.readline()
	if (re.search(r'>>END_MODULE',line) != None):			# exits if no sequences
		return (adap_list)
	line = f_open.readline()
	while (re.search(r'>>END_MODULE',line) == None):
		Sequence,Count,Percentage,Possible_Source = line.split("\t") 
		adap_list.append(Sequence)				# adds adaptor to list
		line = f_open.readline()	
	return (adap_list)
	f_open.close()

def make_cmd(app, options, parameters): #command line command generator for easy to make commands: app is the app u want to use: prin-seq, ect   /   options is a predefined list of options like prin_pe_ops /parameters is a parallel list to optinons of parameter values for that app.         
	cmd = app
	for i in range(len(options)):
		if (len(parameters) > i):
			cmd = cmd + ' ' + options[i] + ' ' + parameters[i]
		else:
			cmd = cmd + ' ' + options[i]
	return(cmd)

def qualencode(pth):# this function may break if not using 
	f_open = open(pth)
	line = f_open.readline()
	while (re.search(r'Encoding',line) == None): # finds Encoding line
		line = f_open.readline()
	title,num,end = line.split("\t")
	try:
		val = float(num[-3:]) #hopefully this wont mess up, it could be characters if not illumina stuff. this try should be able to handel it, default to not phred64
	except ValueError:
		val = 1.2 
	if ((val >= 1.3) & (val<= 1.7)):
		return(1)
	elif ((val < 1.3) | (val > 1.7)):
		return(-1)
	else:
		return(-1)

def make_refgenome(genome,txt,bwa_pth):	# a function to make new reference genomes and add their names to a txt file so they may be used at a later date as well.
					# genome is a path to a fasta file, the name of the file should be descriptive of the organism and will be used in the text file to denote what reference genomes 
					# are avaliable (I am assuming only bacterial genomes right now so it wont be too bad to store all of them locally.)
					# txt is the text file genomes in ref_genomes folder that lists the genomes in there. variable "geno"
					# bwa is the path to bwa, will use the index function from there to make the appropriate index of the fasta file.

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


###################################################################################################

sam_to_bam = robjects.r('''
	sam_to_bam <- function(lst_o_pths, lst_o_outs){
	  k = 1
	  suppressMessages(library(Rsamtools));
	  x = vector("list", length(lst_o_pths) )
	  for (i in 1:length(lst_o_pths))
	  {
	    x[[k]]= asBam(lst_o_pths[[i]],lst_o_outs[[i]])
	    k  = k + 1
	  }
	  x
	}
''')

########### V work to do here!! V ###########
htseqtools = robjects.r('''
	htseqtools <- function(lst_o_pths, lst_o_fldrs){
		suppressMessages(library(rtracklayer))
		suppressMessages(library(htSeqTools))
		my.rangeddatalst = vector("list", length(lst_o_pths))
		for (i in 1:length(lst_o_pths))
		{
			my.rangeddatalst[[i]] <- RangedData(import.bed(con=lst_o_pths[[i]]))
		}
		x = RangedDataList(my.rangeddatalst)
		nme = vector("list", length = length(lst_o_pths))
		for(i in 1:length(lst_o_pths))
		{			
			nme[[i]] = paste("Exo_rep",i,sep = "")
		}
		names(x) = nme

		ssd = ssdCoverage(x)
		write.csv(ssd, paste(lst_o_fldrs[[1]],'/ssd.csv',sep = ''), row.names = lst_o_pths)

		for (i in 1:length(x))
		{
			pdf(paste(lst_o_fldrs[[i]],'/','gini.pdf', sep = ''))
			#file.create(paste(lst_o_fldrs[[i]],'/','gini.pdf', sep = ''))
			giniCoverage(x[[i]],mc.cores=1,mk.plot=TRUE)
			#dev.copy(pdf, paste(lst_o_fldrs[[i]],'/','gini.pdf', sep = ''))
			dev.off()
		}
	}
''')
########################################3

    #----------Code Chunks---------#
#----------Perliminary Questions----------#

geno_open = open(geno)     #this bit gets the reference genomes already on file.
ref_list = []
new = -1
for lin in geno_open:
	ref_list.append(lin)
geno_open.close()
print("Hello,\n\nWelcome to the Chip-Seq/ChIP-exo Pipeline writen by Weston Selna\n\nThis is for use by the IT Lab, any and all reference documents for this Pipeline can be found at /home/linh/Documents/Pipeline_files/pipeline_docs.pdf")
pth1 = input('Please enter the path to the folder containing the .fastq files that you want analyzed\n\nPlease make sure you have reviwed and kept to the file name requierments\n\nFailure to do so will most likely result in an error and this pipeline exiting, or it will somehow work but give you meaningless results.\n\n')

print("The avaliable reference genomes to map to are:")			# this section is for determination of the organism under study.
for acc in ref_list:
	print ("%s" % acc)
s_new = input("Will you be using one of these? (y/n)\n")
if (re.search(r"[y,Y]",s_new)!=None):
	new = 0
	accession = input("Which then? (please type the name in exactly as seen in the prior print out) ") 
	pthgl = [r for r in os.listdir("/media/linh/OS/ref_genomes") if (r == accession)]
	pthg = os.path.join("/media/linh/OS/ref_genomes",pthgl[0])
	print(pthg)
else:
	pthg = input('Enter the Path to the flie containing the fasta files you would like to set as a genome\nThis will be saved for later use so please make sure the file in named appropriately so other may use it too ')# a new organism to align too.
	new = 1						
	if (new == 1):
		print("Making new index of reference genome.")
		pthg = make_refgenome(pthg,geno,bwa)
	
skip = input("skip ahead?(y/n)")
if (re.search(r"[y,Y]",skip) == None):   #does this really do what I think it does?

	#----------Part 1: Preprocessing----------#
	fplist = make_filepath_list(pth1)		#paths to fastq files
	# checking that each path has word chip or control in it. an that there are at least 4 files present. 
	if (len(fplist)<2):
		sys.exit("There are not enough files. The minimum number of files allowed is 2, there should be at least 2 replicates, and correspondingly at least 1 single end read fastq file each totaling 2. More files are allowed above this for extra replicates and PE data") 
	print("beginning FastQC analysis")

	cmd1 = 'fastqc'					# builds a factqc command. but rearranges the names, bad!!!! need a loop it looks like
	for f in fplist:
		cmd1 = cmd1 + " " + f

	os.chdir(fastqc)
	os.system(cmd1)  			

	print("Beginning automated read analysis: filtering and trimming")
	txtplist = []
	files = [f for f in os.listdir(pth1) if os.path.isfile(os.path.join(pth1, f))]
	for a in files:  # getting/making paths to appropriate fastqc_data.txt files in different folders 
		aa = a.replace(".","_")
		aa_p = aa.split()
		aa_p.append("c")
		aaa = ''.join(aa_p)
		if os.path.isdir(os.path.join(pth1, aaa)):		
			txtplist.append(os.path.join(os.path.join(pth1, aaa),"fastqc_data.txt")) #lists still parallel? no they are not!!!

	txtplist = sorted(txtplist)      

	PE_ToF = []
	trim_qual_params = []		# these are important factors from the FastQC report and fastq file.
	#trim_len_params = []		# removed for exo      
	adapter_params = []		
	qual_encoding = []

	for b in txtplist:
		trim_qual_params.append(int(get_trim_qual(b)))
		#trim_len_params.append(int(get_trim_len(b))) unneeded since read ar already so short, will see if need, but may interfer with integrity if used.
		postcut_read_len.append(int(get_trim_len(b)))
		adapter_params.append(get_overrep_seq(b))
		qual_encoding.append(qualencode(b))
	for f in fplist:
		PE_ToF.append(isPE(f))

	#build params lists and run programs.
	switch = 0
	save = None				#these are place saver variables for the main preprocessing 'for' loop.
	save_name = None

	for f in fplist:
		name_out = out_name_fastq(f) 
		if (adapter_params[fplist.index(f)] != []):	# cutadapt
			cut_params = []
			cut_params.append(.1) 							# error allowed in adaptor match. tunable in expert mode, default is .1
			cut_params.append(adapter_params[fplist.index(f)])
			cut_params.append(f+" > "+os.path.join(pth1,name_out+"_cut"+".fastq")) 	# output names of files, need .fastq on end

	
			cut_params.append('2>'+ os.path.join(pth1,name_out+"_cutadapt"+".log"))	# suppressing standard error output from cutadapt. send it to a log file.
	

			cmd2 = apps[0]
			for e in cut_params:
				if (type(e) == list):
					for g in e:
						cmd2 = cmd2 + ' ' + cut_ops[cut_params.index(e)] + ' ' + g + '$' + ' ' + cut_ops[cut_params.index(e)+1] + ' ' + '^' + g # adds all adapters to command.
				elif(type(e)!= float):						# if not -e
					cmd2 = cmd2 + ' ' + e
				else:								#for -e
					cmd2 = cmd2 + ' ' + cut_ops[cut_params.index(e)] + ' ' + str(e)
			os.chdir(cutadapt)
			os.system(cmd2)
			fplist[fplist.index(f)] = os.path.join(pth1,name_out+"_cut"+".fastq")
			f = (os.path.join(pth1,name_out+"_cut"+".fastq"))  							#new path to new adapterless reads	

		if (PE_ToF[fplist.index(f)] == True):
			#start of prin-seq
			if (switch == 1):		
				prin_params = []
				prin_params.extend((f,save,"null","ld,gc,qd,ns,pt,ts,aq,de,sc,dn",'13',str(trim_qual_params[fplist.index(f)]),'0',"entropy",'70'))
				if (qual_encoding[fplist.index(f)] == 1):				#different options depending on encoding of quality
					cmd3 = make_cmd(apps[1],prin_pe_ops,prin_params)
				else:
					cmd3 = make_cmd(apps[1],prin_pe_ops_non64,prin_params)
				cmd3 = cmd3 + " 2> " + os.path.join(pth1,name_out+"_prin"+".log")	#suppressing output to a log file
				os.chdir(prin)
				os.system(cmd3)
				# finding the fastq files in prinseq folder(actually in pth1 folder), moving, renaming and getting their new paths.
				outs = [f for f in os.listdir(pth1) if (re.search(r"prinseq_good",os.path.basename(f)) != None)]  	# got to keep them straight.
				outs[0] = outs[0].split("_")
				if (re.search(outs[0][0],f)):
					outs[0] = "_".join(outs[0])
					outs[0],outs[1] = outs[1],outs[0]		# this if-else maintains order of files.
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
				switch = 1		# encountered the first of the pair of PE fastq files.
				save = f
				save_name = name_out
		else:
			prin_params = []
			prin_params.extend((f,os.path.join(pth1,name_out+"_prin"),"null","ld,gc,qd,ns,pt,ts,aq,de,sc,dn",'13',str(trim_qual_params[fplist.index(f)]),'0',"entropy",'70'))
			if (qual_encoding[fplist.index(f)] == 1):
				cmd3 = make_cmd(apps[1],prin_se_ops,prin_params)
			else:
				cmd3 = make_cmd(apps[1],prin_se_ops_non64,prin_params)
			cmd3 = cmd3 + " 2> " + os.path.join(pth1,name_out+"_prin"+".log")
			os.chdir(prin)
			os.system(cmd3)
			fplist[fplist.index(f)] = os.path.join(pth1,name_out+"_prin.fastq")
	
	print("End of Preprocessing\nYou can find all files and folders made in this process in %s, each with descriptive filenames" %pth1)

	#----------Part 2: Alignment----------#

	print("Starting Alignment process")
	samplist = []
	
	print("using BWA-ALN, best for reads less than 70bp in length\nThe longest reads for the given files after trimming are as follows:")
	for f in fplist:
		print ("for file %s ----------- %d" %(os.path.basename(f), postcut_read_len[fplist.index(f)]))

	both = 0
	save_sai = None			# place holder variables for alignment for loop.
	save_reads = None
	for f in fplist:
		# pe vs se selector/matcher
		if (PE_ToF[fplist.index(f)] == False):
			#SE
			if (qual_encoding[fplist.index(f)] == 1): # different parameters based on quality encoding 
				bwa_params = [apps[3],"-I",pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")] ## 14 is good, yes? Check! ## -14 because descriptive names wer getting too long.
			else:
				bwa_params = [apps[3],pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
			cmd4 = " ".join(bwa_params)					# bwa aln (-I) ref.fa reads.fq > reads.sai
			os.chdir(bwa)
			os.system(cmd4)
			sesai = os.path.join(pth1,os.path.basename(f)[:-14]+".sai")
			bwa_params_samse = [apps[4],pthg,sesai,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_se.sam")]
			cmd5 = " ".join(bwa_params_samse)				# bwa samse ref.fa reads.sai reads.fq > aln-se.sam
			os.chdir(bwa)
			os.system(cmd5)
			samplist.append(os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_se.sam"))  # will not make a perfectly parallel list cause and PE data will only out 1 sam file for 2 fastq files.
		else:	#PE
			if (both == 1):
				if (qual_encoding[fplist.index(f)] == 1):
					bwa_params_2 = [apps[3],"-I",pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				else:
					bwa_params_2 = [apps[3],pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				cmd7 = " ".join(bwa_params_2)				# bwa aln (-I) ref.fa reads.fq > reads.sai
				os.chdir(bwa)
				os.system(cmd7)
				pesai = os.path.join(pth1,os.path.basename(f)[:-14]+".sai")
				both = 0
				bwa_params_sampe = [apps[5],pthg,save_sai,pesai,save_reads,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_pe.sam")]
				cmd8 = " ".join(bwa_params_sampe)			# bwa sampe ref.fa reads.sai reads.fq > aln-pe.sam
				os.chdir(bwa)
				os.system(cmd8)
				samplist.append(os.path.join(pth1,os.path.basename(f)[:-14]+"_aln_pe.sam"))# will not make a perfectly parallel list cause and PE data will only out 1 sam file for 2 fastq files. the sam file will have the name of the second pair that goes into this fcn.
			else:
				if (qual_encoding[fplist.index(f)] == 1):
					bwa_params_1 = [apps[3],"-I",pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				else:
					bwa_params_1 = [apps[3],pthg,f,">",os.path.join(pth1,os.path.basename(f)[:-14]+".sai")]
				cmd6 = " ".join(bwa_params_1)				# bwa aln (-I) ref.fa reads.fq > reads.sai
				os.chdir(bwa)
				os.system(cmd6)
				save_sai = os.path.join(pth1,os.path.basename(f)[:-14]+".sai")
				save_reads = f
				both = 1

#--------Part 3: Post QC--------#

print("Starting post alignment, pre peak calling quality control")
print("converting to other file types for analysis") 

# need sam out names.
os.chdir(pth1)
samoutnames = []
for sam in samplist:					# fix this so naming is appropriate.
	bse = os.path.basename(sam)
	if (re.search(r"sorted",bse) != None):     # difference between PE and SE.
		b_name = bse[:-17] 
	else:
		b_name = bse[:-10]
	samoutnames.append(os.path.join(pth1,b_name+"asbam")) #still parallel?   what about .bam?(seemed to be fine, revisit documentation of sam_to_bam functions in R if needed)    

bamplist_0 = sam_to_bam(samplist,samoutnames) #make bams, are R objects here
bamplist = []
for i in range(0, len(samoutnames)):
	bamplist.append(bamplist_0[i][0]) #get paths to bams, as strings here (not r objects), add .bai to end to get sorted/indexed versions

bedplist = []
bamoutnames = []
#convert to bed for htSEQtools. makecmd for bedtools				
for bm in bamplist:
	bm_out = os.path.basename(bm)			
	bm_out = bm_out[:-9]+"as_bed.bed"
	bm_out = os.path.join(pth1,bm_out)
	bedplist.append(bm_out)

	bedtool_params = [bm,bm_out]
	os.system(make_cmd(apps[6],bedtools_ops,bedtool_params))  #makes beds (tell them they can view them in IGV)

print("Corresponding BED files have been made from the alignments. These .bed files can be viewed in IGV or on the UCSC genome browser. See the additional help text for how to use IGV to locally view the alignments\n")

#htSEQtools R stuff: make a folder for pictures from the stuff, generate a list of folder paths to put plots in.
rplotfolderlist = []
for flder in bedplist:   # makes folders for R plots, will be named so differentiation between samples is easy. 
	fldr_nme = os.path.basename(flder)[:-18]			
	os.mkdir(pth1+"/Rplots_"+fldr_nme)
	rplotfolderlist.append(pth1+"/Rplots_"+fldr_nme)	

print("Calculating Gini index for all samples to see if immunopercipitation worked\nCalculating the coverage standard deviaton\n")

htseqtools(bedplist,rplotfolderlist)			#Statistical stuff

print("Coverage and immunopercipitation statistics and plots are avaliable to view in the newly crated Rplots folders in the Main folder holding the original fastq files.")

#####-------Part 4: Peak Calling---------#####
os.chdir(pth1)

#samtools method to get names and lengths of chromosomes
smtools_idx_params = [bamplist[0]]
sizes_file = os.path.basename(pthg)+".sizes"
ct_params = ["1,2", sizes_file]
cmd10 = apps[7] + " " + smtools_idx_params[0] # as a string?
cmd11 = make_cmd(apps[8],ct_ops,ct_params)

#print(cmd10, end ='')
#print(" | ", end ='')
#print(cmd11)
cmd1011 = cmd10+" | "+ cmd11				# this should work, but sometimes it fucks up.
os.system(cmd1011)
															
#MACE preprocessing 
prefix = os.path.basename(pth1)+"_MACE"
wigs = [os.path.join(pth1,prefix+"_Forward.wig"),os.path.join(pth1,prefix+"_Reverse.wig")]
prepros_params = [",".join(bamplist), sizes_file, prefix, "0", os.path.join(pth1,"preprocessing.log")]
cmd12 = make_cmd(apps[9],prepros_ops,prepros_params)

os.system(cmd12)

#bigwig conversion    
os.chdir(pth1)
bigwigs = [wigs[0][:-3]+"bw", wigs[1][:-3]+"bw"]
if ("Reverse" in bigwigs[0]):
	bigwigs[0],bigwigs[1] = bigwigs[1],bigwigs[0] #error checking since order is imperative, Reverse must always be Second in the list of 2.

cmd13 = apps[10] + " " + wigs[0] + " " + sizes_file+ " " +bigwigs[0]
cmd14 = apps[10] + " " + wigs[1] + " " + sizes_file+ " " +bigwigs[1]
os.system(cmd13)
os.system(cmd14)

#MACE call
#example: mace.py -s hg19.chrom.sizes -f CTCF_MACE_Forward.bw -r CTCF_MACE_Reverse.bw -o CTCF_MACE -p .05
mace_params = [sizes_file,bigwigs[0],bigwigs[1], prefix, ".05", "mace.log"]
cmd15 = make_cmd(apps[11], mace_ops, mace_params)
os.system(cmd15)


print("Thank You for using the IT lab ChIP-seq/ChIP-exo Pipeline!\nYour files should all be located in the original Main folder you indicated at the beginning of the process.\n\nThis pipeline has created a wealth of useful files for you to review at your lesiure.\nConfer with the documentation and Usage Manual for specifics on each file.\nGood Bye\n")
	
