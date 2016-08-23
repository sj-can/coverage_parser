#creates integers from transscript intervals - i.e. genomic coordinates of every locus

ranges = []
range_array = []
in_file = []
not_in_file = []
pos_iter_array = []
neg_iter_array = []

def generate_transcript_range(input_file):
    with open(input_file, 'r') as input:
	for line in input:
            split_line = line.split(',')
            ranger = (int(split_line[4]) - int(split_line[3])) + 1
            ranges.append(ranger)
            for i in range(int(split_line[3]), int(split_line[4]) + 1, 1):
	        i = ['11', str(i)]
                range_array.append(i)
    return range_array

#array of targets in format [[c, g],[c,g]] and coverge file in order of chromosome
def binary_search_coverage(locus_array, coverage_file):
    with open(coverage_file, 'r') as coverage_file:
	items_checked = 0
	this_item = 0
	for item in locus_array:
	    this_item += 1
	    last_target_locus = []
	    end_array = []
	    begin_array = [0]
	    coverage_file.seek(0,2)
	    end_array.append(coverage_file.tell())
	    while items_checked < this_item:
	            this_begin = begin_array[-1]
	            this_end = end_array[-1]
	            coverage_file.seek(((this_end + this_begin)/2),0)
	            partial_line = coverage_file.readline()
	            after_partial = coverage_file.tell()
	            coverage_file.seek(after_partial)
	            whole_line = coverage_file.readline()
	            split_line = whole_line.split()
	            if ':' in split_line[0]:
	                split_locus = split_line[0].split(':')
	                chrom = int(item[0])
	                target_chrom = int(split_locus[0])
        	        locus = int(item[1])
                	target_locus = int(split_locus[1])
			if chrom == target_chrom:
			    last_target_locus = [target_locus]
        	            chrom_end_array = [end_array[-1]]
	                    chrom_begin_array = [begin_array[-1]]
			    neg_iter = 0
			    pos_iter = 0
			    while locus != last_target_locus[-1]:
                                coverage_file.seek(((chrom_end_array[-1] + chrom_begin_array[-1])/2),0)
	                        partial_line2 = coverage_file.readline()
	                        end_partial2 = coverage_file.tell()
	                        coverage_file.seek(end_partial2)
	                        whole_line2 = coverage_file.readline()
	                        split_line2 = whole_line2.split()
        	                split_locus2 = split_line2[0].split(':')
                	        target_chrom2 = int(split_locus2[0])
                        	target_locus2 = int(split_locus2[1])
				if chrom == target_chrom2:
				    if locus == target_locus2:
					in_file.append(locus)
				        last_target_locus.append(target_locus2)
					items_checked += 1
					if neg_iter > 1:
					    neg_iter_array.append(neg_iter)
					if pos_iter > 1:
					    pos_iter_array.append(pos_iter)
				    elif locus < target_locus2:
					chrom_end_array.append(coverage_file.tell())
					if chrom_end_array[-1] == chrom_end_array[-2]:
					    chrom_begin_array.append(chrom_begin_array[-1] - 1)
					    neg_iter += 1
					    if neg_iter >= 200:
					        not_in_file.append(locus)
						last_target_locus.append(locus)
						items_checked += 1
						neg_iter_array.append(neg_iter)
				    elif locus > target_locus2:
					chrom_begin_array.append(coverage_file.tell())
					if chrom_begin_array[-1] == chrom_begin_array[-2]:
					    chrom_end_array.append(chrom_end_array[-1] + 1)
					    if pos_inter >= 200:
						not_in_file.append(locus)
						last_target_locus.append(locus)
						items_checked += 1
						pos_iter_array.append(pos_iter)
				elif chrom > target_chrom2:
				    #move beginning to current_location	
				    chrom_begin_array.append(coverage_file.tell())
				    #reset the end to the one before it stars repeating
				    if chrom_end_array[-2] == chrom_end_array[-1] - 1:
					end_iter1 = 2
					end_iter2 = 1
					while chrom_end_array[-iter1] == chrom_end_array[-iter2] - 1:
					    end_iter1 += 1
					    end_iter2 += 2
				    	chrom_end_array.append(end_iter1)
				    else:
					del chrom_end_array[-1]
				elif chrom < target_chrom2:
				    #move end to current location
				    chrom_end_array.append(coverage_file.tell())
				    #reset the begining to before it starts repeating
				    if chrom_begin_array[-2] == chrom_begin_array[-1] - 1:
                                        begin_iter1 = 2
                                        begin_iter2 = 1
                                        while chrom_beign_array[-iter1] == chrom_begin_array[-iter2] - 1:
                                            begin_iter1 += 1
                                            begin_iter2 += 2
                                        chrom_begin_array.append(begin_iter1)
                                    else:
                                        del chrom_begin_array[-1]

			elif chrom > target_chrom:
			    #move begining to current_location
			    begin_array.append(coverage_file.tell())
			elif chrom <target_chrom:
			    #move end to current_location
			    end_array.append(coverage_file.tell())

generate_transcript_range("transcript_intervals/NM_000352.3.intervals")
search = binary_search_coverage(range_array, "/mnt/Data5/exome_sequencing/WE0343-WE0350/r01_metrics/Coverage")
print len(in_file)
print len(not_in_file)
print not_in_file
