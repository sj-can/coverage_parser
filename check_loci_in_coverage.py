

import subprocess

#neeed to incorporate_chromosome from object to remove hard coded while loop 

ranges = []
range_array = []
locus_range_array = [['11', '17498305'],['11', '17498129']]
'''
with open("transcript_intervals/NM_000352.3.intervals", 'r') as input:
    for line in input:
        split_line = line.split(',')
        ranger = (int(split_line[4]) - int(split_line[3])) + 1
        ranges.append(ranger)
	for i in range(int(split_line[3]), int(split_line[4]) + 1, 1):
            range_array.append(i)
    for i in range_array:
    	i = ['11',str(i)]
    	locus_range_array.append(i)

    print len(range_array)
    print locus_range_array
    print len(locus_range_array)
'''

#binary search function
in_file = []
not_in_file = []
with open('/mnt/Data5/exome_sequencing/WE0343-WE0350/r01_metrics/Coverage', 'r') as coverage_file, open('not_in_coverage.txt', 'w') as out_file:
    items_checked = 0
    this_item = 0
    for item in locus_range_array:
	print item
	this_item += 1
	end_array = []
	begin_array = [0]
	last_target_locus = []
        coverage_file.seek(0,2)
	end_array.append(coverage_file.tell())
	while items_checked < this_item:
	    this_begin = begin_array[-1]
	    print begin_array
	    this_end = end_array[-1]
	    print end_array
	    print 'this_begin = ' +str(begin_array[-1]) + ' this end = ' +str(end_array[-1]) 
            coverage_file.seek(((this_end + this_begin)/2),0)
            partial_line = coverage_file.readline()
	    after_partial = coverage_file.tell()
	    coverage_file.seek(after_partial)
	    whole_line = coverage_file.readline()
	    split_line = whole_line.split()
	    print split_line
	    print 'in_file = ' + str(in_file)
	    if ':' in split_line[0]:
                split_locus = split_line[0].split(':')
		print item[1] + ' ' + split_locus[1]
		print item[0] + ' ' + split_locus[0]
		chrom = int(item[0])
		target_chrom = int(split_locus[0])
                locus = int(item[1])
                target_locus = int(split_locus[1])
	        if chrom > target_chrom:
		    print 'chrom bigger than target'
		    new_begin = coverage_file.tell()
		    begin_array.append(new_begin)
		elif chrom < target_chrom:
		    print 'chrom smaller than target'
		    new_end = coverage_file.tell()
		    end_array.append(new_end)
		elif chrom == target_chrom:
		    print 'found chrom #####################################################################'
		    last_target_locus = [target_locus]
		    chrom_end_array = [end_array[-1]]
                    chrom_begin_array = [begin_array[-1]]
		    #items_checked = items_checked + 1
		    offset = 0
   		    while locus != last_target_locus[-1]:
			print locus
			#print last_target_locus
		    	#print chrom_begin_array
			#print chrom_end_array
			coverage_file.seek(((chrom_end_array[-1] + chrom_begin_array[-1])/2 - offset),0)
			partial_line = coverage_file.readline()
			print partial_line
			end_partial2 = coverage_file.tell()
			coverage_file.seek(end_partial2)
			whole_line2 = coverage_file.readline()
			split_line2 = whole_line2.split()
                        split_locus2 = split_line2[0].split(':')
                        target_chrom2 = int(split_locus2[0])
                       	target_locus2 = int(split_locus2[1])
			last_target_locus.append(target_locus2)
			print split_line2
			print target_chrom2
			print target_locus2
			print chrom_begin_array 
			print chrom_end_array
			print
			#find the start of the chromosome and make that begin	    
			if chrom == target_chrom2:
				print 'infile:'
				print in_file
				if locus > target_locus2:				    
				    #bring the start to present
				    chrom_begin_array.append(coverage_file.tell())
				    #leave the end where it is
				    if chrom_begin_array[-1] == chrom_begin_array[-2]:
					while chrom_begin_array[-1] == chrom_end_array[-2]:
						print '#####################################################################################' 
				elif locus < target_locus2:
				    #bring the end closer
				    chrom_end_array.append(coverage_file.tell())
				    #leave the start where it is
				    if chrom_end_array[-1] == chrom_end_array[-2]:
					print 'end_array the same'
					print 'partial_line'
					print partial_line
					chrom_begin_array.append(chrom_begin_array[-1] - 1)
					print chrom_begin_array
					print chrom_begin_array[-1] 
				elif locus == target_locus2:
					last_target_locus.append(target_locus2)
					in_file.append(locus)
					items_checked += 1
					print 'this_item : ' + str(this_item)
					print 'items_checked : ' + str(items_checked)
					print in_file
					
			elif chrom > target_chrom2:
			    #move begining to where you are now
			    chrom_begin_array.append(coverage_file.tell())
			    #reset the end to the one before the last
			    del chrom_end_array[-1]
			elif chrom < target_chrom2:
				for i in range(1000):
				    print 'error chrom < target_chrom2'

			#dive the 


'''
                	target_locus2 = int(split_locus[1])
			if chrom2 < target_chrom2:
			    chrom_begin_array.append(coverage_file.tell())
			    chrom_end_array.append(end_array[0])
		        if locus == target_locus and chrom2 == target_chrom2:
                            in_file.append(item)
                            print 'Found locus in file on correct chromosome' + str(item)
			    items_checked += 1
			    last_target_locus.append(target_locus)
		        elif locus > target_locus:
			    print 'locus is bigger moving beginning to current position'
			    last_target_locus.append(target_locus)
			    begin_array.append(coverage_file.tell())
		        elif locus < target_locus:
			    chrom_end_array.append(coverage_file.tell())
                            last_target_locus.append(target_locus)
    	                    print 'this_begin = ' +str(chrom_begin_array[-1]) + ' this end = ' +str(chrom_end_array[-1])
		            print 'locus is smaller moving end to current position and checking last locus'
			    #if target_locus == last_target_locus[-1]:
			    #	print 'this target locus is equal to the last target locus - negative offset requried'
			    
'''
'''
		elif int(item[1]) == int(split_locus[1]) and int(item[0]) > int(split_locus[0]):
		    print 'found locus on a preceeding chromosome to the target'
		    new_begin = coverage_file.tell()
		    begin_array.append(new_begin)

		elif int(item[1]) == int(split_locus[1]) and int(item[0]) < int(split_locus[0]):
		    print 'found locus on a subsequent chromosome to the target'
		    new_end = coverage_file.tell()
		    end_array.append(new_end)
                elif int(item[1]) > int(split_locus[1]) and int(item[0]) == int(split_locus[0]):
                    print 'coord is higer coverage_file and chromosome the same'
                    new_begin = coverage_file.tell()
                    begin_array.append(new_begin)
                    print
                elif int(item[1]) > int(split_locus[1]) and int(item[0]) > int(split_locus[0]):
		    print 'coord is higher and chromosome is higher than targets'
                    new_begin = coverage_file.tell()
                    begin_array.append(new_begin)
                elif int(item[1]) > int(split_locus[1]) and int(item[0]) < int(split_locus[0]):
		    print 'locus is higher and chromosome is lower than targets'
                    new_end = coverage_file.tell()
                    end_array.append(new_end)

                elif int(item[1]) < int(split_locus[1]) and int(item[0]) == int(split_locus[0]):
                    print 'item < line in coverage_file correct chromosome'
                    new_begin = coverage_file.tell()
                    begin_array.append(new_begin)

                elif int(item[1]) < int(split_locus[1]) and int(item[0]) > int(split_locus[0]):
                    print 'item < line in coverage_file greater than chromosome'
                    new_begin = coverage_file.tell()
                    begin_array.append(new_begin)	
                elif int(item[1]) < int(split_locus[1]) and int(item[0]) < int(split_locus[0]):
		    print 'locus is lower and chromosome is lower'
                    new_end = coverage_file.tell()
                    end_array.append(new_end)
		else:
		    print '###########################################################################################################################################################################################################'#
		    not_in_file.append(item)
		    break		    
'''	    


