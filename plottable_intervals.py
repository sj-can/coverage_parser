#creates two files that are used to represent the exons and +/- 50 for each. The two files are then plotted on the same graph.
#adds lines ot data taht are plotted as non and therefore no line. 

#file 1: nnnnnnnnnnn<exon1>nnnnn<exon2>
#file 2: <single in>nnnnnnn<si >nnnnnnn

with open('transcript_intervals/NM_000352.3.intervals') as interval_file, open('NM_000352.3_exons', 'a') as exons, open('NM_000352.3_extended', 'a') as extended, open('master_intervals', 'a') as master:
    for line in reversed(interval_file.readlines()):
        split_line = line.split(',')
        #print split_line
        exon_range = int(split_line[2])-int(split_line[1]) + 1
        extension = 50
        for i in range(extension):
            exons.write('\n')
        for i in range(exon_range):
            exons.write('-0.05\n')
        for i in range(extension):
            exons.write('\n')
        for i in range(extension):
            extended.write('-0.05\n')
        for i in range(exon_range):
            extended.write('\n')
        for i in range(extension):
            extended.write('-0.05\n')

#        for i in range(extension):
#            master.write()
#        for i in range(exon_range):
#            master.write()
#        for i in range(extension):
#            master.write()

        




    
    
