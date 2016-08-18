#adds line number to an outputfile. takes first column of input

with open("checker_no_blanks.txt", 'r') as input, open("checker_no_blanks_enum.txt", 'w') as output:
    for i,line in enumerate(input):
        if line != '\n':
            line = line.split(',')
            if '\n' not in line[0]:
                output.write(str(i) + ',' + line[0] + '\n')
            else:
                output.write(str(i) + ',' + line[0])
        elif line == '\n':
            output.write('\n')
        
