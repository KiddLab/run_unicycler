import sys
import os
import gzip
import shutil
import subprocess

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
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
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
# read top line of blat PSL file
def read_top_blat_line(pslFileName):
    inFile = open(pslFileName,'r')
    for i in range(0,5): # skip over header
        line = inFile.readline() 
    line = inFile.readline()
    inFile.close()
    line = line.rstrip()
    line = line.split()
    blatLine = parse_blat_psl_line(line)
    return blatLine
#####################################################################
# assumes line is already list
def parse_blat_psl_line(line):
    blatLine = {}
    blatLine['match'] = int(line[0])
    blatLine['mismatch'] = int(line[1])
    blatLine['repmatch'] = int(line[2])
    blatLine['ncount'] = int(line[3])
    blatLine['qGapCount'] = int(line[4])
    blatLine['qGapBases'] = int(line[5])
    blatLine['tGapCount'] = int(line[6])
    blatLine['tGapBases'] = int(line[7])
    blatLine['strand'] = line[8]
    blatLine['qName'] = line[9]
    blatLine['qSize'] = int(line[10])
    blatLine['qStart'] = int(line[11])
    blatLine['qEnd'] = int(line[12])

    blatLine['tName'] = line[13]
    blatLine['tSize'] = int(line[14])
    blatLine['tStart'] = int(line[15])
    blatLine['tEnd'] = int(line[16])
    
    blatLine['blockCount'] = int(line[17])
    blatLine['blockSizes'] = line[18]
    blatLine['qStarts'] = line[19]
    blatLine['tStarts'] = line[20]
    return blatLine
#####################################################################
##############################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################




