import subprocess

#neeed to incorporate_chromosome from object to remove hard coded while loop 

ranges = []
range_array = [[11,17498126]]
#locus_range_array = []
#with open("transcript_intervals/NM_000352.3.intervals", 'r') as input:
#    for line in input:
#        split_line = line.split(',')
#        ranger = (int(split_line[4]) - int(split_line[3])) + 1
#        ranges.append(ranger)
#	for i in range(int(split_line[3]), int(split_line[4]) + 1, 1):
#            range_array.append(i)
    #for i in range_array:
    #	i = '11:' + str(i)
    #	locus_range_array.append(i)

    #print len(range_array)
    #print range_array[0]
    #print range_array[-1]

#binary search function
in_file = []
not_in_file = []
with open('/mnt/Data5/exome_sequencing/WE0343-WE0350/r01_metrics/Coverage', 'r') as coverage_file, open('not_in_coverage.txt', 'w') as out_file:
    for item in range_array:
	end_array = []
	begin_array = [0]
        coverage_file.seek(0,2)
	end_array.append(coverage_file.tell())
        while begin_array[-1] < end_array[-1]:
	    this_begin = begin_array[-1]
            #print 'begin_array ' +str(begin_array)
            #print 'this begin ' + str(begin_array[-1])

	    this_end = end_array[-1]
            #print 'end_array ' + str(end_array)

	    print 'this_begin = ' +str(begin_array[-1]) + ' this end = ' +str(end_array[-1]) 
	    #half the current_distance
            coverage_file.seek(((this_end + this_begin)/2),0)
            partial_line = coverage_file.readline()
	    current = coverage_file.tell()
	    coverage_file.seek(current)
	    whole_line = coverage_file.readline()
	    split_line = whole_line.split()	

	    try:
	        if ':' in split_line[0]:
	            split_locus = split_line[0].split(':')
	            if int(item[1]) == int(split_locus[1]) and int(item[0]) == int(split_locus[0]):
                       in_file.append(item)
		       print 'in_file : ' + str(item)
		       print
	            elif int(item) > int(split_locus[1]):
		        print 'item > line in coverage_file'
		        new_begin = coverage_file.tell()
			print new_begin
			begin_array.append(new_begin)
			print
	            elif int(item) < int(split_locus[1]):
			print 'item < line in coverage_file'
		        new_end = coverage_file.tell()
			end_array.append(new_end)
			print end_array
			print 
	    except:
	        for i in range(100):
		    print
	    #try:
	    #    split_locus = split_line[0].split(':')
	    #    print split_locus
            #    if split_locus[1] == item:
	   #	    in_file.append(item)
	#	    #print in_file
	 #   	elif int(item) > int(split_locus[1]):
	#	    begin = coverage_file.tell()
	 #       elif int(item) < int(split_locus[1]):
	  #          end = coverage_file.tell()
	  #  except:
	#	pass
	        #print 'couldnt split_locus'
	        #print item


#for item in range_array:
#    if item not in in_file:
#	out_file.write(item)
	

	
    




'''
with open('/mnt/Data5/exome_sequencing/WE0343-WE0350/r01_metrics/Coverage', 'r') as coverage_file, open('not_in_coverage.txt', 'w') as out_file:
    range_array[0] = "not_in_file"
    for item in range_array:
	item = '11:' + str(item)
        command = ["grep", item, "/mnt/Data5/exome_sequencing/WE0343-WE0350/r01_metrics/Coverage"]
	process = subprocess.Popen(command, stdout=subprocess.PIPE)
        x = process.communicate()
	if x[0] == '':
	    out_file.write(str(item) + '\n')
'''

