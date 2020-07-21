#fasta assocationed stuff

def get_fasta(name, db):
	if db == 'file1':	
		seq = file1_dict[name]
		fasta = ">" + name + "\n" + seq

	if db == 'file2':	
		seq = file2_dict[name]
		fasta = ">" + name + "\n" + seq

	if db == 'file3':	
		seq = file3_dict[name]
		fasta = ">" + name + "\n" + seq

	if db == 'file4':	
		seq = file4_dict[name]
		fasta = ">" + name + "\n" + seq


	return fasta

def write_fasta(file_name,contents):
	outfile = open(file_name, "w")
	outfile.write(contents)
	outfile.close()

def fasta_dict(fastafile):
	name2seq = {}
	file = open(fastafile)

	for line in file:
		if line.startswith(">"):
			seq = ''
			split_line = line.split(' ')
			name = split_line[0]
			name = name.rstrip()
			name = name.lstrip('>')
		else:
			seq = seq + str(line)
			seq = seq.rstrip()
			seq = seq.upper()

		name2seq[name] = seq
	
	return name2seq

def fasta_list(fastafile):
	file = open(fastafile)
	fasta_list = []

	for line in file:
		if line.startswith(">"):
			seq = ''
			split_line = line.split(' ')
			name = split_line[0]
			name = name.rstrip()
			name = name.lstrip('>')
		else:
			seq = seq + str(line)
			seq = seq.rstrip()
			seq = seq.upper()

		fasta_line = name + "\n" + seq

		fasta_list.append(fasta_line)
	
	return fasta_list
