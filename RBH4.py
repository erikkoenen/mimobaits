"""reciprocal top blast hit for four fasta files, gathering the data on %identity and length of match for each sucessful hit """

'usage: python RBH4.py <cdhit_threshold> <taxon_A> <taxon_B> <taxon_C> <taxon_D> <number_of_cpus>'

import sys
import subprocess
import fasta_stuff_rbh4

# specify cdhit threshold

cdhit_threshold = sys.argv[1]

# taxa

taxon_A = sys.argv[2]
taxon_B = sys.argv[3]
taxon_C = sys.argv[4]
taxon_D = sys.argv[5]

# cpus

cpus = sys.argv[6]

# define blast functions

def blastp_A():
	my_command = "blastp -query to_blast.fasta -db " + taxon_A + ".pep." + cdhit_threshold + ".SC -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length evalue sstart send' -out process_blast"
	output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)

def blastp_B():
	my_command = "blastp -query to_blast.fasta -db " + taxon_B + ".pep." + cdhit_threshold + ".SC -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length evalue sstart send' -out process_blast"
	output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)

def blastp_C():
	my_command = "blastp -query to_blast.fasta -db " + taxon_C + ".pep." + cdhit_threshold + ".SC -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length evalue sstart send' -out process_blast"
	output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)

def blastp_D():
	my_command = "blastp -query to_blast.fasta -db " + taxon_D + ".pep." + cdhit_threshold + ".SC -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length evalue sstart send' -out process_blast"
	output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)

def blast_parse(blast_output):
	file = open(blast_output)
	split_line = "no hit"
	if len(blast_output) != 0:
		for line in file:
			line = line.rstrip('\n')
			split_line = line.split('\t')
	return split_line



#Gets a fasta file for a specfic name from a dict

def get_fasta(name, db):
        if db == 'file1':
                seq = file1_dict[name]
                fasta = ">" +name + "\n" + seq

        if db == 'file2':
                seq = file2_dict[name]
                fasta = ">" +name + "\n" + seq

        if db == 'file3':
                seq = file3_dict[name]
                fasta = ">" +name + "\n" + seq

        if db == 'file4':
                seq = file4_dict[name]
                fasta = ">" +name + "\n" + seq

        return fasta

# sets up a blast search for a named fasta sequence to a particular database

def Find_hits(fasta, db):
	fasta_stuff_rbh4.write_fasta('to_blast.fasta', fasta)
	if db == 'A':
		blastp_A()
	if db == 'B':
		blastp_B()  
	if db == 'C':
		blastp_C()
	if db == 'D':
		blastp_D()
	hit = blast_parse('process_blast')

	return hit

# command for alignment with MAFFT linsi option (suitable algorithm for aligning sequences with alignable domains flanked by non-alignable regions, change "8" to the number of cores on your machine)
def mafft(fasta_input, output_name):
	my_command = "mafft --maxiterate 1000 --globalpair --thread " + cpus + " " +str(fasta_input) +" > " +str(output_name)
	output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)

# trim alignments to keep only reliably algined regions with BMGE
def bmge(output_name, output2_name):
	my_command = "java -jar ~/Software/BMGE-1.1/BMGE.jar -i "+str(output_name) +" -t AA -h 1 -g 0.51:1 -o "+str(output2_name)
	output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)

# run RAxML to build the trees (change "8 to the number of cores on your machine)
def raxml(output2_name,outputID):
	try:
		my_command = "raxmlHPC-PTHREADS-SSE3 -T" + cpus + " -f a -m PROTCATLGF -p 12345 -x 12345 -# 100 -s "+str(output2_name)+" -n "+str(outputID)
		output = subprocess.check_output(my_command, shell=True, stderr=subprocess.STDOUT)
	except:
		print ("Problems with RAxML")

#*******************START SCRIPT *********************

fastafile1 = taxon_A + ".pep." + cdhit_threshold + ".SC"
fastafile2 = taxon_B + ".pep." + cdhit_threshold + ".SC"
fastafile3 = taxon_C + ".pep." + cdhit_threshold + ".SC"
fastafile4 = taxon_D + ".pep." + cdhit_threshold + ".SC"

#Make dicts of the fasta sequences in each fiel wiht he fasta-line as the key

file1_dict = fasta_stuff_rbh4.fasta_dict(fastafile1)
file2_dict = fasta_stuff_rbh4.fasta_dict(fastafile2)
file3_dict = fasta_stuff_rbh4.fasta_dict(fastafile3)
file4_dict = fasta_stuff_rbh4.fasta_dict(fastafile4)

# Prepare output files

output = ['rbh	A	B	C	D	A_B_pident	A_B_length	A_C_pident	A_C_length	B_C_pident	B_C_length	A_D_pident	A_D_length	B_D_pident	B_D_length	C_D_pident	C_D_length']
output2 = ['rbh	A	B	C	A_B_pident	A_B_length	A_C_pident	A_C_length	B_C_pident	B_C_length']
reciprocal_count = 0
lines_processed = 0

# Aim is to generate a file with names of reciprocal blast hits, percentage IDs and lengths of those hits

# Start with the genes in the A list

for name in file1_dict:
	RBH3 = "no"	
	lines_processed += 1
	A_fasta = ">" + name + "\n" + file1_dict[name]

# 1: blast the A to B      
	AB_hit = Find_hits(A_fasta, 'B')
	if AB_hit != "no hit":
		B_hit_name = AB_hit[1]
		B_hit_fasta = get_fasta(B_hit_name, 'file2')
		AB_recip_hit = Find_hits(B_hit_fasta, 'A')
		if AB_recip_hit[1] == AB_hit[0]:
			print("Found AB match!")
			
			
# 2: blast the B to C
			BC_hit = Find_hits(B_hit_fasta, 'C')
			if BC_hit != "no hit":
				C_hit_name = BC_hit[1]
				C_hit_fasta = get_fasta(C_hit_name, 'file3')
				BC_recip_hit = Find_hits(C_hit_fasta, 'B')
				if BC_recip_hit[1] == BC_hit[0]:
					print("Found BC match!")

					
# 3: blast the A to C and check C match is same as A to B
					AC_hit = Find_hits(A_fasta, 'C')
					if AC_hit != "no hit":
						C_hit_name = AC_hit[1]
						C_hit_fasta = get_fasta(C_hit_name, 'file3')
						AC_recip_hit = Find_hits(C_hit_fasta, 'A')
						if AC_recip_hit[1] == AC_hit[0] and AC_hit[1] == BC_hit[1]:
							print("Found ABC match!")
							RBH3 = "yes"
							reciprocal_count += 1
# 4: blast the C to D
							CD_hit = Find_hits(C_hit_fasta, 'D')
							if CD_hit != "no hit":
								D_hit_name = CD_hit[1]
								D_hit_fasta = get_fasta(D_hit_name, 'file4')
								CD_recip_hit = Find_hits(D_hit_fasta, 'C')
								if CD_recip_hit[1] == CD_hit[0]:
									print ("Found CD match!")
											
# 5: blast the A to D and check D match is same as C to D
									AD_hit = Find_hits(A_fasta, 'D')
									if AD_hit != "no hit":
										D_hit_name = AD_hit[1]
										D_hit_fasta = get_fasta(D_hit_name, 'file4')
										AD_recip_hit = Find_hits(D_hit_fasta, 'A')
										if AD_recip_hit[1] == AD_hit[0] and AD_hit[1] == CD_hit[1]:
											print ("Found AD match!")

# 6: blast the B to D and check D match is same as C to D and A to D
											BD_hit = Find_hits(B_hit_fasta, 'D')
											if BD_hit != "no hit":
												D_hit_name = BD_hit[1]
												D_hit_fasta = get_fasta(D_hit_name, 'file4')
												BD_recip_hit = Find_hits(D_hit_fasta, 'B')
												if BD_recip_hit[1] == BD_hit[0] and BD_hit[1] == CD_hit[1]:
													print ("Found ABCD match!")

													RBH3 = "no"

# Now we have a set of four rbhs, get the info together for the output

													A = ">" + taxon_A + "\n" + file1_dict[name]
													B = ">" + taxon_B + "\n" + file2_dict[B_hit_name]
													C = ">" + taxon_C + "\n" + file3_dict[C_hit_name]
													D = ">" + taxon_D + "\n" + file4_dict[D_hit_name]
													multi_fasta = A + "\n" + B + "\n" + C + "\n" + D
													file_name = "rbh" + str(reciprocal_count)+ ".fasta"
													path = "fasta_RBH4/"+file_name
													fasta_stuff_rbh4.write_fasta(path, multi_fasta)

													print ("aligning rbh"+str(reciprocal_count))
													fasta_input = "fasta_RBH4/rbh"+str(reciprocal_count)+".fasta"
													output_name = "mafft_RBH4/rbh"+str(reciprocal_count)+"_aln.fa"
													mafft(fasta_input, output_name)
													
													print ("trimming alignment of rbh"+str(reciprocal_count))
													output2_name = "bmge_RBH4/rbh" + str(reciprocal_count) + ".phy"
													bmge(output_name, output2_name)
													
													print ("building rapid bootstrap tree for rbh"+str(reciprocal_count))
													outputID = "rbh"+str(reciprocal_count)
													raxml(output2_name,outputID)


# set up outputfile
													output_line = "rbh" + str(reciprocal_count) + "\t" + str(AB_hit[0]) + "\t" + str(AB_hit[1]) + "\t" + str(AC_hit[1]) +  "\t" + str(AD_hit[1]) + "\t" + str(AB_hit[2]) + "\t" + str (AB_hit[3]) +  "\t" + str(AC_hit[2]) + "\t" + str (AC_hit[3])+  "\t" + str(BC_hit[2]) + "\t" + str(BC_hit[3]) + "\t" + str(AD_hit[2]) + "\t" + str(AD_hit[3]) + "\t" + str(BD_hit[2]) + "\t" + str (BD_hit[3])  + "\t" + str(CD_hit[2]) + "\t" + str(CD_hit[3])
													output.append(output_line)

# if we find only 3 rbhs (so not for D), then we print this to an output file as well
							if RBH3 == "yes":
								output_line2 = "rbh3_" + str(reciprocal_count) + "\t" + str(AB_hit[0]) + "\t" + str(AB_hit[1]) + "\t" + str(AC_hit[1]) +  "\t" + str(AB_hit[2]) + "\t" + str (AB_hit[3]) +  "\t" + str(AC_hit[2]) + "\t" + str (AC_hit[3])+  "\t" + str(BC_hit[2]) + "\t" + str (BC_hit[3])
								output2.append(output_line2)


print("sequences with reciprocol best hits = " + str(reciprocal_count))
print("Lines processed = " + str(lines_processed))

outfile = open('RBH4_output', "w")
outfile.write("\n".join(output))
outfile.close()
outfile = open('RBH3_output', "w")
outfile.write("\n".join(output2))
outfile.close()
        
