import argparse, os, fnmatch, subprocess

#argument parser
class Exome_coverage_parser():

    exon_objects = []

    def __init__(self, transcript_class):
        parser = argparse.ArgumentParser()
        parser.add_argument('-s', action='store', dest='sample_names', help='path to file containing exome samples to be analysed')
        parser.add_argument('-g', action='store', dest='gene_names', help='path to file containing the genes to be analysed')
        arguments = parser.parse_args()
        self.arg_dict = vars(arguments)
        self.alamut_file_instance = transcript_class

    def exome_coverage_finder(self):
        exome_coverage_files = []
	sample_iteration = 0
        exome_sample_list = self.list_from_file(self.arg_dict['sample_names'])
        print 'exome_sample_list: '
        for exome_sample in exome_sample_list:
	    print 'searching for ' + exome_sample + 'coverage file'
            match_count = 0
            for root, directories, filenames in os.walk('/mnt/Data5/exome_sequencing/'):
		for filename in filenames:
		    match_pattern = '*' + exome_sample + '*'
	            if fnmatch.fnmatch(str(filename), match_pattern) and match_count == 0:
			match_count += 1
			if len(exome_coverage_files) == sample_iteration:
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
        print str(exome_coverage_files) + '\n'
        return exome_coverage_files

    def line_strip_split(self, line):
        line = line.strip('\n')
        split_line = line.split()
        return split_line

    #opens a file and turns each line into a lsit object
    def list_from_file(self, input_file):
        new_list = []
        with open(input_file, 'r') as inputfile:
            for line in inputfile:
                line = line.strip()
                new_list.append(line)
        return new_list

    #parse the transcript start_stop and exon_intervals from the
    #/mnt/Data1/resources/alamut-genes/
    #return { gene_symbol : '', transcript: '', transcript_coords : [start, stop], exons :{ transcript : '', [start, stop, CDS_boolean]}}
    #create the two separately then add them together to pass to the plotting function
    def gene_interval_parser(self, list_of_genes):
	print 'calculating exon intervals from alamut coverage file \n'
        with open('/mnt/Data1/resources/alamut-genes/grch37_2016-05-10.txt', 'r') as alamut_genes:
            for line in alamut_genes:
                for gene in list_of_genes:
                    #if its the first line
                    if gene in line and len(self.exon_objects) == 0:
                        line = line.split()
                        initial_transcript_instance = self.create_new_instance(line)                       
                        self.exon_objects.append(initial_transcript_instance)
                    #if its not the first line
                    elif gene in line and len(self.exon_objects) > 0:
                        line = line.split()
                        if self.exon_objects[-1].transcript_id == line[7] and line[12] not in self.exon_objects[-1].exons:
                            self.exon_objects[-1].exons[line[12]] = [line[13], line[14]]
                        if self.exon_objects[-1].transcript_id != line[7]:
                            additional_transcript_instance = self.create_new_instance(line)
                            self.exon_objects.append(additional_transcript_instance)

    #creates a new instance of a transcript from a line in the alamut file
    def create_new_instance(self, line):
        instance = self.alamut_file_instance(line[1], line[7], line[8], line[9], line[3], line[4], line[5])
        instance.exons[line[12]] = [line[13], line[14]]
	return instance

    #identifies the transcript object with the lingest transcript
    def longest_transcript(self):
        longest_transcripts = []
        current_longest = []
        #iterate through the gene list provided
        print 'identifying longest transcripts'
        for gene in self.list_from_file(self.arg_dict['gene_names']):
            this_gene = []
            current_longest = []
            #iterate through the all transcript list of objects
            for item in self.exon_objects:
                if gene == item.__dict__['gene_symbol']:
                    if len(current_longest) == 0:
                        this_gene.append(item)
                        current_longest.append(item)
                    elif len(current_longest) == 1:
                        this_gene.append(item)                     
            for item in this_gene:
                current = current_longest[-1].__dict__['transcript_length']
                challenger = item.__dict__['transcript_length']
                if int(challenger) > int(current):
                    current_longest[-1] = item
            longest_transcripts.append(current_longest[-1])
        for item in longest_transcripts:
            print str(item.__dict__) + '\n'
        return longest_transcripts

    #create a file containing the transcripts to be plotted
    def exon_interval_file_creator(self, transcript_instance):
	print 'creating exon interval files'
        output_name = transcript_instance.__dict__['transcript_id']
        sorted_name = transcript_instance.__dict__['transcript_id'] + '.intervals'
        with open(output_name, 'a') as output_file:
            exon_dictionary = transcript_instance.__dict__['exons']
            for k,v in exon_dictionary.iteritems():
                output = k + ',' + v[0] + ',' + v[1] + '\n'
                output_file.write(output)

    def sorter(self, file):
        output = 'transcript_intervals/' + file +'.intervals'
        with open(output, 'a') as sorted_file:
            command = ["sort", "-V", file]
            process = subprocess.Popen(command, stdout=subprocess.PIPE)
            sorted_output = process.communicate()[0]
            sorted_file.write(sorted_output)
            remove_command = ['rm', file]
            subprocess.call(remove_command)

    def parse_coverage_file(self, list_of_transcript_instances, exome_id, coverage_file_location):
        for item in list_of_transcript_instances:
            cds_length = self.exon_counter(item)
            with open(coverage_file_location, 'r') as coverage_data:
                header_list = coverage_data.readline().split()
                print header_list
                print coverage_file_location
                match_sample_column = self.match_sample_to_column(header_list, exome_id)
                print match_sample_column
                locus_count = 0
                locus_array = []
                print str(item) + 'iteration'
                instance_file_name = exome_id + '_' + item.__dict__['gene_symbol']
                print instance_file_name
                with open(instance_file_name, 'a') as output:
                    print 'opened_output'
                    output.write('>' + instance_file_name + '\n')
                    for line in coverage_data:
                        if line.startswith(item.__dict__['chromosome']):
                            split_line = line.split()
                            locus = split_line[0].split(':')
                            if int(locus[1]) >= int(item.__dict__['transcript_start']) and int(locus[1]) <= int(item.__dict__['transcript_end']):
                                output.write(locus[1] + ',' + str(split_line[match_sample_column[0]] + '\n')) 
                                locus_count += 1
                                locus_array.append(split_line[match_sample_column[0]])
                                print exome_id
                                print header_list
                                print split_line
                                print 'depth: ' + str(split_line[match_sample_column[0]])
                                print match_sample_column[1]
			        print str(locus_count)
				print len(locus_array)
				print cds_length
				print str(item.__dict__['transcript_length'])
                            else:
                                pass


    def parse_coverage_filezzz(self, list_of_transcript_instances, exome_id, coverage_file_location):
        with open(coverage_file_location, 'r') as coverage_data:
            header_list = coverage_data.readline().split()
	    print header_list
	    print coverage_file_location
	    match_sample_column = self.match_sample_to_column(header_list, exome_id)
            print match_sample_column
            for item in list_of_transcript_instances:
                locus_count = 0
                locus_array = []
                print str(item) + 'iteration'
                instance_file_name = exome_id + '_' + item.__dict__['gene_symbol']
                print instance_file_name
                with open(instance_file_name, 'a') as output:
                    print 'opened_output'
                    for line in coverage_data:
	                if line.startswith(item.__dict__['chromosome']):
			    split_line = line.split()
			    locus = split_line[0].split(':')
                            if int(locus[1]) >= int(item.__dict__['transcript_start']) and int(locus[1]) <= int(item.__dict__['transcript_end']):
                                locus_count += 1
                                locus_array.append(split_line[match_sample_column[0]])
                                print exome_id
				print split_line[match_sample_column[0]]
                                print match_sample_column[1]
				print header_list
		        else:
	                    pass
		    

    def match_sample_to_column(self, header_list, exome_identifier):
        print header_list
        print exome_identifier
        y = 'Depth_for_' + exome_identifier        
        for i, x in enumerate(header_list):  
	    if x == y:
		result = [i,x]
	return result

    def exon_counter(self, instance):
        exon_total = 0
        for k,v in instance.__dict__['exons'].iteritems():
            distance = int(v[1]) - int(v[0])
	    exon_total + distance
	return exon_total
     
class Transcript_instance():

    def __init__(self, gene_symbol, transcript_id, transcript_start, transcript_end, chromosome, g_start, g_stop):
        self.gene_symbol = gene_symbol
        self.transcript_id = transcript_id
        self.transcript_start = transcript_start
        self.transcript_end = transcript_end
        self.transcript_length = int(transcript_end) - int(transcript_start)
        self.chromosome = chromosome
        self.gene_start = g_start
        self.gene_stop = g_stop
        self.exons = {}

#don't plot the genomic coordinate - include a buffer of 50 either side

instance = Exome_coverage_parser(Transcript_instance)
list_of_coverage_file_dicts = instance.exome_coverage_finder()
gene_intervals = instance.gene_interval_parser(instance.list_from_file(instance.arg_dict['gene_names']))
list_of_longest_transcripts = instance.longest_transcript()
for item in list_of_longest_transcripts:
    instance.exon_interval_file_creator(item)
    instance.sorter(item.__dict__['transcript_id'])
for item in list_of_coverage_file_dicts:    
    for k,v in item.iteritems():
        instance.parse_coverage_file(list_of_longest_transcripts, k, v)
    






