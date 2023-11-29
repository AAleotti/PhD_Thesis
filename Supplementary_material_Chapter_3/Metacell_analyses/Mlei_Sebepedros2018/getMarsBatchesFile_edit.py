import os

fNames = sorted([fName[:fName.find('.')] for fName in os.listdir('raw_umi_tables')])

f = open('MARS_Batches.txt', 'w')

header = ['Amp.Batch.ID', 'Seq.Batch.ID', 'Batch.Set.ID', 'dataset', 'color']

delimiter = '\t'

f = open('MARS_Batches.txt', 'w')
f.write(delimiter.join(header)+'\n')

for fName in fNames:
    f.write(delimiter.join(3*[fName]+['Adult', 'black'])+'\n')

f.close()
print('Done. Please see MARS_Batches.txt.')
