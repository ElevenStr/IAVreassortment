fread = open('influenza.cds','r')
fwrite = open('influenza1Line.cds','w')
# fread = open('influenza.faa','r')
# fwrite = open('influenza1Line.faa','w')
lines = fread.readlines()
fread.close()
for i in range(len(lines)):
    line = lines[i].rstrip()
    if '>' in line:
        if i > 0:
            fwrite.write(seq + '\n')
        fwrite.write(line + '\n')
        seq = ''
        seqInfo = line
    else:
        seq = seq + line
fwrite.write(seq + '\n')
fwrite.close()