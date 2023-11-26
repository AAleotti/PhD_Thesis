import os
import xlwt
import xlrd
 
stageDB_file_name = '1-s2.0-S0092867418305968-mmc1.xlsx'
dbFilePath = os.getcwd()+'/'+stageDB_file_name
umiFolderName = 'raw_umi_tables' 
adultOnly = True

if adultOnly:
    outFileName = 'MARS_Batches_Adult_Only.txt'
else:
    outFileName = 'MARS_Batches_All_Stages.txt'



def readDB(filePath):
    db = dict()
    insheets = xlrd.open_workbook(filePath).sheets()
    for i in range(len(insheets)):
        insheet = insheets[i]
        for rowID in range(1,insheet.nrows):
            db[insheet.cell_value(rowID,0)] = insheet.cell_value(rowID,1)
    return db

def writeFile(outFileName, umiFolderName, stageDB):
    fNames = sorted([fName[:fName.find('.')] for fName in os.listdir(umiFolderName)])

    header = ['Amp.Batch.ID', 'Seq.Batch.ID', 'Batch.Set.ID', 'dataset', 'color']

    delimiter = '\t'

    f = open(outFileName, 'w')
    f.write(delimiter.join(header)+'\n')

    for fName in fNames:
        if adultOnly:
           if 'Adult' in stageDB[fName] or 'Juvenile' in stageDB[fName]:        
            f.write(delimiter.join(3*[fName]+[stageDB[fName], 'black'])+'\n')
        else:
            f.write(delimiter.join(3*[fName]+[stageDB[fName], 'black'])+'\n')
    
    f.close()
    print('Done. Please see '+outFileName+'.')

def main():
    stageDB = readDB(dbFilePath)
    writeFile(outFileName, umiFolderName, stageDB)

main()
