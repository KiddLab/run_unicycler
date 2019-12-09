import sys
import os
import gzip
import shutil

# python3 codes to run and process BAC assembly with unicycler
# needs Illumina + long reads

###############################################################################
def add_breaks_to_line(seq,n=50):
    myList = []
    myList = [i for i in seq]
    newList = []
    c = 0
    for i in myList:
        newList.append(i)
        c += 1
        if c % n == 0 and c != (len(myList)):
            newList.append('\n')
    myStr = ''.join(newList)
    return myStr    
###############################################################################
# Helper function to run commands,
# doesn't check return or print to log.  Use for 'grep' so that doesn't
# fail if no items found
def runCMDNoFail(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
###############################################################################
def read_fasta_file_to_list(fastaFile):
    myDict = {}
    inFile = open(fastaFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print('ERROR, FILE DOESNNOT START WITH >')
        sys.exit()
    myName = line[1:]
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    mySeq = ''
    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            mySeq = ''
            continue
        mySeq += line
    inFile.close()
    return myDict
###############################################################################
# poor name choice on function when first started.  Keep to save compatability
def read_fasta_file_to_dict(fastaFile):
    return read_fasta_file_to_list(fastaFile)
###############################################################################
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):    
    print('Checking bwa...')
    if shutil.which('bwa') is None:
        print('bwa not found in path! please fix (module load?)')
        sys.exit()

    print('Checking samtools...')
    if shutil.which('samtools') is None:
        print('samtools not found in path! please fix (module load?)')
        sys.exit()

    print('Checking blat...')
    if shutil.which('blat') is None:
        print('blat not found in path! please fix (module load?)')
        sys.exit()

    print('Checking RepeatMasker...')
    if shutil.which('RepeatMasker') is None:
        print('RepeatMasker not found in path! please fix (module load?)')
        sys.exit()
#####################################################################







