f = open('lukowski_embo2019_raw_count_matrix.csv', 'r')
outf = open('cleaned_matrix.tsv', 'w')
counter = 0
for line in f:
    if counter > 0:
        newline = line.replace(' ','').replace(',,',',0,').replace(',','\t')
        outf.write(newline)
    counter += 1
f.close()
outf.close()
