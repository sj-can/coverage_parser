import argparse, os, subprocess, fnmatch

#initialises command line arguments and calls handler function that ensures all arguments have been passed
class Argument_handler():
    def __init__(self):
        parser = argparse.ArgumentParser()
	#parser.add_argument(alamut_file_path)
        parser.add_argument('-s', action='store', dest='sample_names', help='path to file containing exome samples to be analysed')
        parser.add_argument('-g', action='store', dest='gene_names', help='path to file containing the genes to be analysed')
        arguments = parser.parse_args()
        self.arg_dict = vars(arguments)
        print 'initialised exome coverage parser with: ' + str(self.arg_dict)
	self.handler()

#handle if no sample path is given nonetype error
    def handler(self):
	switch = 0
	if 'sample_names' in self.arg_dict.keys() and os.path.isfile(self.arg_dict['sample_names']):
	    print 'samples file location confirmed. Ensure one sample per line.'
	    switch += 1
        elif not os.path.isfile(self.arg_dict['sample_names']):
	    print 'Cannot find smaples file check path or ass -s flag before path'
	if 'gene_names' in self.arg_dict.keys() and os.path.isfile(self.arg_dict['gene_names']):
	    print 'genes file location confirmed. Ensure one sample per line.'
	    switch += 1
	elif not os.path.isfile(self.arg_dict['sample_names']):
            print 'Cannot find smaples file check path or add -g flag before path'
        if switch == 2:
#	    call rest of program
	    instance = Coverage_parser(self.arg_dict['sample_names'], self.arg_dict['gene_names'], '/mnt/Data1/resources/alamut-genes/grch37_2016-05-10.txt')
            instance.exome_coverage_finder()
	    instance.find_gene_intervals()
	    instance.longest_transcript()
	    
	else:	
	    print 'ERROR - Input criteria not satisfied.'

#Models a transcript using data from the Alamut genes file 
class Transcript():

    def __init__(self, alamut_line):
        self.alamut_line = alamut_line
	self.gene_symbol = self.alamut_line[1]
        print 'gene_symbol = ' + self.gene_symbol
        self.transcript_id = self.alamut_line[7]
        print 'trans_id = ' + self.transcript_id
        self.transcript_start = self.alamut_line[8]
        self.transcript_end = self.alamut_line[9]
        self.transcript_length = int(self.transcript_end) - int(self.transcript_start) + 1
        self.chromosome = self.alamut_line[3]
        self.gene_start = self.alamut_line[4]
        self.gene_stop = self.alamut_line[5]
        self.exons = {self.alamut_line[12] : [self.alamut_line[13], self.alamut_line[14]]}

class Coverage_parser(Transcript):

    transcript_instances = []
    exome_coverage_files = ''

    def __init__(self, sample_name_file, gene_name_file, alamut_file):
        self.sample_list = self.list_from_file(sample_name_file)
        self.gene_list = self.list_from_file(gene_name_file)
	self.alamut_file = alamut_file

    def line_strip_split(self, line):
        line = line.strip('\n')
        split_line = line.split()
        return split_line

    def list_from_file(self, input_file):
        new_list = []
        with open(input_file, 'r') as inputfile:
            for line in inputfile:
                line = line.strip()
                new_list.append(line)
        return new_list

    def exome_coverage_finder(self):
	exome_coverage_files = []
	print 'finding coverage files for : '
        print (self.sample_list)
	sample_iteration = 0
	for exome_sample in self.sample_list:
            print 'searching for ' + exome_sample + ' coverage file'
	    match_count = 0
            for root, directories, filenames in os.walk('/mnt/Data5/exome_sequencing/'):
                for filename in filenames:
                    match_pattern = '*' + exome_sample + '*'
                    if fnmatch.fnmatch(str(filename), match_pattern) and match_count == 0:
                        match_count += 1
                        if len(self.exome_coverage_files) == sample_iteration:
                            sample_iteration += 1
                            if 'r01_metrics' in root:
                                print '------match------'
                                coverage_dict = {}
                                coverage_dict[exome_sample] = root + '/Coverage'
                                print coverage_dict
                                exome_coverage_files.append(coverage_dict)
                                print
                            else:
                                match_count -= 1
                                sample_iteration -= 1
			else:
			    pass
        print 'exome_coverage_files: ' + str(exome_coverage_files) + '\n'
        self.exome_coverage_files = exome_coverage_files
	return exome_coverage_files

    def find_gene_intervals(self):
        lines_parsed = 0
        print 'calculating exon intervals from alamut file \n'
        with open(self.alamut_file, 'r') as alamut_file:
            for line in alamut_file:
                for gene in self.gene_list:
                    #if its the first line
                    if gene in line and len(self.transcript_instances) == 0:
                        line = line.split()
                        initial_transcript_instance = Transcript(line)  
			self.transcript_instances.append(initial_transcript_instance)
		    elif gene in line and len(self.transcript_instances) > 0:
		        line = line.split()
 			if self.transcript_instances[-1].transcript_id == line[7] and line [12] not in self.transcript_instances[-1].exons:
			    self.transcript_instances[-1].exons[line[12]] = [line[13], line [14]]
			if self.transcript_instances[-1].transcript_id != line[7]:
			    additional_transcript_instance = Transcript(line)
			    self.transcript_instances.append(additional_transcript_instance)

#need to account for it they are equal
    def longest_transcript(self):
        print 'finding longest transcript'
        longest_transcripts = []
	current_longest = []
        #iterate through the gene list provided
        print 'identifying longest transcripts'
        for gene in self.sample_list:
            print gene
            this_gene = []
            #current_longest = []
            #iterate through the all transcript list of objects
            for item in self.transcript_instances:
                if gene == item.__dict__['gene_symbol']:
                    if len(current_longest) == 0:
                        this_gene.append(item)
                        current_longest.append(item)
                    elif len(current_longest) == 1:
                        this_gene.append(item)
            for item in this_gene:
                current = current_longest[-1].__dict__['transcript_length']
                challenger = item.__dict__['transcript_length']
                if int(challenger) >= int(current):
                    current_longest[-1] = item
            longest_transcripts.append(current_longest[-1])
        print 'the longest transcripts are:'
        for item in longest_transcripts:
            print str(item.__dict__['gene_symbol']) + ' ' + str(item.__dict__['transcript_id']) + '\n'
        return longest_transcripts



Argument_handler()
	
