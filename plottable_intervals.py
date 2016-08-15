with open('transcript_intervals/NM_000352.3.intervals') as interval_file, open('NM_000352.3_exons', 'a') as exons, open('NM_000352.3_extended', 'a') as extended:
    for line in reversed(interval_file.readlines()):
        split_line = line.split(',')
        #print split_line
        exon_range = int(split_line[2])-int(split_line[1])
        extension = 50
        for i in range(extension):
            exons.write('\n')
        for i in range(exon_range):
            exons.write('5\n')
        for i in range(extension):
            exons.write('\n')
        for i in range(extension):
            extended.write('10\n')
        for i in range(exon_range):
            extended.write('\n')
        for i in range(extension):
            extended.write('10\n')




    
    
