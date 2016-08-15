with open("NM_000352.3_extended", 'r') as input, open("NM_000352.3_extended_2col", 'w') as output:
    for i,line in enumerate(input):
        if line != '\n':
            line = line.split(',')
            output.write(str(i) + ',' + line[0])
        elif line == '\n':
            output.write('\n')
        
