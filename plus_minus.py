#need to check that the difference between exons is > 100
#need to check if the sequence needs reversing or not

interval_array = []
with open('transcript_intervals/NM_000352.3.intervals', 'r') as intervals:
    for line in intervals:
        line = line.split(',')
        interval_array.append([line[3], line[4]])

#with open('WE0345_ABCC8', 'r') as input, open('WE0345_ABCC8_ext', 'w') as output:
#    completed_intervals = []
#    for line in input:
#        split_line = line.split(',')
#        print split_line
#        for interval in interval_array:
#            if int(split_line[0]) >= int(interval[0]) and int(split_line[0]) <= int(interval[1]):
#                output.write(line)
#                if interval not in completed_intervals:
#		    completed_intervals.append(interval)

with open('WE0345_ABCC8_ext_gaps', 'w') as output:
    completed_intervals = []
    for interval in reversed(interval_array):
        with open('WE0345_ABCC8', 'r') as input:
            print completed_intervals
            print interval
            for line in input:
	        split_line = line.split(',')
                interval_1 = interval[1].strip('\n')
	        if int(split_line[0]) >= int(interval[0]) and int(split_line[0]) <= int(interval_1):                
                    output.write(line)
	    if interval not in completed_intervals:
	        completed_intervals.append(interval)
	        output.write('\n\n\n')
