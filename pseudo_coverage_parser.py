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

    #give it a list of exomes to be analysed, do this step once and return a list of matches
    def exome_coverage_finder(self):
        exome_coverage_files = []
        gene_names_list = self.list_from_file(self.arg_dict['gene_names'])
        exome_sample_list = self.list_from_file(self.arg_dict['sample_names'])
        #locate the metrics directory  within a batch directory
        for root, directories, filenames in os.walk('/mnt/Data5/exome_sequencing/'):
            for filename in filenames:
                for exome_sample in exome_sample_list:
                    match_pattern = '*' + exome_sample + '*'
   		    if fnmatch.fnmatch(str(filename), match_pattern):
                        if 'r01_metrics' in directories:
                            coverage_dict = {}
                            #print root
                            coverage_dict[exome_sample] = root + '/r01_metrics/Coverage'
                            exome_coverage_files.append(coverage_dict)
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
        instance = self.alamut_file_instance(line[1], line[7], line[8], line[9])
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
        return longest_transcripts

    #now I know which transcript I need to use by deafult,
    #I need to feed the transcript info (exons too into) the data to be plotted on the graph for each exome sample
        
    #create a file containing the transcripts to be plotted
    def exon_interval_file_creator(self, transcript_instance):
        output_name = transcript_instance.__dict__['transcript_id']
        sorted_name = transcript_instance.__dict__['transcript_id'] + '.intervals'
        with open(output_name, 'a') as output_file:
            exon_dictionary = transcript_instance.__dict__['exons']
            print exon_dictionary
            print
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
       
class Transcript_instance():

    def __init__(self, gene_symbol, transcript_id, transcript_start, transcript_end):
        self.gene_symbol = gene_symbol
        self.transcript_id = transcript_id
        self.transcript_start = transcript_start
        self.transcript_end = transcript_end
        self.transcript_length = int(transcript_end) - int(transcript_start)
        self.exons = {}

instance = Exome_coverage_parser(Transcript_instance)
#instance.exome_coverage_finder()
instance.gene_interval_parser(instance.list_from_file(instance.arg_dict['gene_names']))
list_of_longest_transcripts = instance.longest_transcript()
for item in list_of_longest_transcripts:
    instance.exon_interval_file_creator(item)
    instance.sorter(item.__dict__['transcript_id'])






