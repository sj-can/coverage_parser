#need to check that the difference between exons is > 100
#need to check if the sequence needs reversing or not

#takes the interval file and sample coverae data for a gene to outputs a subset that only includes coverage
#plus or minus 50 from the ends of the exons

interval_array = []
with open('transcript_intervals/NM_000352.3.intervals', 'r') as intervals:
    for line in intervals:
        line = line.split(',')
        interval_array.append([line[3], line[4]])

with open('WE0345_ABCC8_pm50_1', 'w') as output:
    completed_intervals = []
    print len(completed_intervals)
    for interval in reversed(interval_array):
        with open('WE0345_ABCC8', 'r') as input:
            #print len(completed_intervals)
            print interval
            for line in input:
	        split_line = line.split(',')
                interval_1 = interval[1].strip('\n')
	        if int(split_line[0]) >= int(interval[0]) and int(split_line[0]) <= int(interval_1):                
                    output.write(line)
	    if interval not in completed_intervals:
	        completed_intervals.append(interval)
	        output.write('\n\n\n')
            else:
                print 'error is in here!'

print len(interval_array)
print len(completed_intervals)


