#adds line number to an outputfile. takes first column of input

with open("WE0345_ABCC8_pm50_1", 'r') as input, open("WE0345_ABCC8_pm50_enum_1", 'w') as output:
    for i,line in enumerate(input):
        if line != '\n':
            line = line.split(',')
            if '\n' not in line[0]:
                output.write(str(i) + ',' + line[0] + '\n')
            else:
                output.write(str(i) + ',' + line[0])
        elif line == '\n':
            output.write('\n')
        
