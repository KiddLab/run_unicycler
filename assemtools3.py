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
def get_4l_record(myFile):
    #fastq style file...
    # just return sequence len
    # -1 if last record
    myLine1 = myFile.readline()
    if myLine1 == '':
        return ''
    myLine2 = myFile.readline()
    myLine3 = myFile.readline()
    myLine4 = myFile.readline()
    return [myLine1,myLine2,myLine3,myLine4]
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
    print('Checking samtools...')
    if shutil.which('samtools') is None:
        print('samtools not found in path! please fix (module load?)')
        sys.exit()

    print('Checking blat...')
    if shutil.which('blat') is None:
        print('blat not found in path! please fix (module load?)')
        sys.exit()

    print('Checking cutadapt...')
    if shutil.which('cutadapt') is None:
        print('cutadapt not found in path! please fix (module load?)')
        sys.exit()

    print('Checking picard environmental variable...')
    if 'PICARD_JARS' not in os.environ:
        print('PICARD_JARS not set! please fix (module load?)')
        sys.exit()

    print('Checking minimap2...')
    if shutil.which('minimap2') is None:
        print('minimap2 not found in path! please fix (module load?)')
        sys.exit()

    # requirements for running unicycler
    for n in ['bowtie2','racon','spades','tblastx','unicycler']:
        print('Checking %s...' % n)
        if shutil.which(n) is None:
            print('%s not found in path! please fix (module load?)' % n)
            sys.exit()

    print('Checking Pilon environmental variable...')
    if 'PILON_JARS' not in os.environ:
        print('PILON_JARS not set! please fix (module load?)')
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
def parse_paf_line(line):
    pafLine = {}
    pafLine['qName'] = line[0]
    pafLine['qLen'] = int(line[1])
    pafLine['qStart'] = int(line[2])
    pafLine['qEnd'] = int(line[3])
    pafLine['strand'] = line[4]
    pafLine['tName'] = line[5]
    pafLine['tLen'] = int(line[6])
    pafLine['tStart'] = int(line[7])
    pafLine['tEnd'] = int(line[8])
    pafLine['numMatch'] = int(line[9])
    pafLine['alignBlockLen'] = int(line[10])
    pafLine['mapQ'] = int(line[11])
    return pafLine

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
def run_cutadapt(myData):
    myData['cutadpt.fq1_1'] = myData['outDirBase'] + 'lib1.R1.cutadapt.fq.gz'
    myData['cutadpt.fq1_2'] = myData['outDirBase'] + 'lib1.R2.cutadapt.fq.gz'


    cmd = 'cutadapt -j %i --match-read-wildcards --discard-trimmed -g GATCGGAAGAGC -G GATCGGAAGAGC ' % myData['numThreads']
    cmd += ' -o %s -p %s %s %s' % (myData['cutadpt.fq1_1'], myData['cutadpt.fq1_2'],myData['fqR1'],myData['fqR2'])     
    if os.path.isfile(myData['cutadpt.fq1_1']) is False:
        print('running cutadpt for Illumina reads')
        print(cmd)
        runCMD(cmd)
    else:
        print('looks like cutadpt already ran for Illumina reads')
#####################################################################
def filter_contam_illumina(myData):
    tmpBam = myData['outDirBase'] + 'tmp.ilm.bam'
    tmpSam = myData['outDirBase'] + 'tmp.ilm.sam'
    myData['ilmContamBam'] = myData['outDirBase'] + 'ilmContam.sort.bam'
    myData['ilmContamBam_met'] = myData['outDirBase'] + 'ilmContam.sort.size_metrics.txt'    
    myData['ilmContamBam_hist'] = myData['outDirBase'] + 'ilmContam.sort.size_metrics.pdf'        
    myData['ilmContamBam_mappednames'] = myData['outDirBase'] + 'ilmContam.sort.mapped_names.txt'            
    
    myData['ilmFilt_1'] = myData['outDirBase'] + 'ilmFilter.R1.fq.gz'             
    myData['ilmFilt_2'] = myData['outDirBase'] + 'ilmFilter.R2.fq.gz'            
    
    if os.path.isfile(myData['ilmFilt_1']) and os.path.getsize(myData['ilmFilt_1']) > 0:
        print('Already ran Illumina filter contamination!')
        return


    bowTie2Index = myData['contam'].replace('.fa','')
    

    cmd = 'bowtie2 -x ' + bowTie2Index + ' -p %i ' % myData['numThreads']
    cmd += ' -1 %s -2 %s -S %s' % (myData['cutadpt.fq1_1'],myData['cutadpt.fq1_2'],tmpSam)
    
    print('\n%s\n' % cmd)
    runCMD(cmd)
    
    cmd = 'samtools view -b %s > %s' % (tmpSam,tmpBam)
    print('\n%s\n' % cmd)
    runCMD(cmd)
    
    cmd = 'rm %s' % tmpSam
    print('\n%s\n' % cmd)
    runCMD(cmd)
        
    cmd = 'java -Xmx4g -jar $PICARD_JARS/picard.jar SortSam I=%s O=%s SORT_ORDER=coordinate' % (tmpBam,myData['ilmContamBam'])
    print(cmd)
    runCMD(cmd)

    cmd = 'samtools index %s' % (myData['ilmContamBam'])
    print(cmd)
    runCMD(cmd)

    cmd = 'rm %s' % tmpBam
    print(cmd)
    runCMD(cmd)
       
    cmd = 'samtools view -F 4 %s | cut -f 1 | sort | uniq > %s' % (myData['ilmContamBam'],myData['ilmContamBam_mappednames'])
    print(cmd)
    runCMD(cmd)

    toDrop = {}
    inFile = open(myData['ilmContamBam_mappednames'],'r')
    for line in inFile:
        line = line.rstrip()
        toDrop[line] = 1
    inFile.close()
    print('Read in %i illumina read names to drop' % len(toDrop))

    inFq1 = gzip.open(myData['cutadpt.fq1_1'],'rt')
    inFq2 = gzip.open(myData['cutadpt.fq1_2'],'rt')
    outFq1 = gzip.open(myData['ilmFilt_1'],'wt')
    outFq2 = gzip.open(myData['ilmFilt_2'],'wt')    
    
    numDrop = 0
    numWrite = 0
    
    while True:
        r1 = get_4l_record(inFq1)
        r2 = get_4l_record(inFq2)
        if r1 == '':
            break
        
        n1 = r1[0].rstrip()
        n1 = n1[1:].split()[0]
        n2 = r2[0].rstrip()
        n2 = n2[1:].split()[0]
        
        if n1 != n2:
            print('Names do not match!',n1,n2)
            sys.exit()
        
        if n1 in toDrop:
            numDrop +=1
        else:
            numWrite +=1
            outFq1.write('%s%s%s%s' % (r1[0],r1[1],r1[2],r1[3]))
            outFq2.write('%s%s%s%s' % (r2[0],r2[1],r2[2],r2[3]))
            
    inFq1.close()
    inFq2.close()
    outFq1.close()
    outFq2.close()
    print('\nSearched for contamination in short reads')
    print('Total read pairs: %i' % (numDrop+numWrite) )
    print('Removed: %i  %f' % (numDrop, numDrop/(numDrop+numWrite)))
    print('Kept: %i  %f' % (numWrite, numWrite/(numDrop+numWrite)))
        
#####################################################################
def filter_contam_longread(myData):
    myData['longReadContam'] = myData['outDirBase'] + 'longread.contam.paf'
    myData['longReadFilt'] = myData['outDirBase'] + 'longread.filt.fq.gz'   
    myData['longReadFiltFail'] = myData['outDirBase'] + 'longread.fail.fq.gz'   
    
    # check to see if already ran
    if os.path.isfile(myData['longReadFilt']) and os.path.getsize(myData['longReadFilt']) > 0:
        print('Already ran long read filter contamination!')
        return
    
    if myData['longreadtype'] == 'ont':
       mapX = 'map-ont'
    else:
       print('uknown long read type, option not clear')
       print(myData['longreadtype'])
       sys.exit()
    
    cmd = 'minimap2 -t 1 -x %s %s %s > %s ' % (mapX,myData['contam'],myData['longread'],myData['longReadContam'] )
    print(cmd)
    runCMD(cmd)
    
    toDrop = {}
    inFile = open(myData['longReadContam'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pafLine = parse_paf_line(line)
        
        if((pafLine['numMatch']/pafLine['qLen']) > 0.25):
            toDrop[pafLine['qName']] = 0
    inFile.close()
    print('Read in %i names of long reads' % len(toDrop))
    
    inFile = gzip.open(myData['longread'],'rt')
    outPass = gzip.open(myData['longReadFilt'],'wt')
    outFail = gzip.open(myData['longReadFiltFail'],'wt')
    
    numDrop = 0
    numWrite = 0
    
    while True:
        r1 = get_4l_record(inFile)
        if r1 == '':
            break
        
        n1 = r1[0].rstrip()
        n1 = n1[1:].split()[0]
        
        if n1 in toDrop:
            numDrop +=1
            outFail.write('%s%s%s%s' % (r1[0],r1[1],r1[2],r1[3]))

        else:
            numWrite +=1
            outPass.write('%s%s%s%s' % (r1[0],r1[1],r1[2],r1[3]))
    inFile.close()
    outPass.close()
    outFail.close()
    print('\nSearched for contamination in long reads!')
    print('Total long reads: %i' % (numDrop+numWrite) )
    print('Removed: %i  %f' % (numDrop, numDrop/(numDrop+numWrite)))
    print('Kept: %i  %f' % (numWrite, numWrite/(numDrop+numWrite)))
#####################################################################
def run_unicycler_assem(myData):
    myData['originalAssem'] = myData['outDirBase'] + 'assembly.fasta'
    myData['assemFa'] = myData['originalAssem']    
    if os.path.isfile(myData['originalAssem']) is True:
        print('assembly exists, skipping run_unicycler_assem step!')
        return        
    print('\nstarting to run unicycler!')
    cmd = 'unicycler -t %i -1 %s -2 %s -l %s -o %s' % (myData['numThreads'],myData['ilmFilt_1'],myData['ilmFilt_2'],myData['longReadFilt'],myData['outDirBase'])
    print(cmd)
    runCMD(cmd)
#####################################################################
def do_rotate_circle(myData):
    myData['rotatedFa'] = myData['originalAssem'] + '.rotate.fa'
    myData['rotatedCleanFa'] = myData['originalAssem'] + '.rotate.clean.fa'

    print('initial assembly:',myData['assemFa'])
    print('portion to set at start:',myData['targetFa'])
    print('make version without the target:',myData['doClean'])

    # run blat of target vs assembly, to determine what the correct orientation is
    myData['blatOutFile'] = myData['assemFa'] + '.blatTMP'
    cmd = 'blat %s %s %s' % (myData['assemFa'],myData['targetFa'],myData['blatOutFile'])
    print(cmd)
    runCMD(cmd)
    blatLine = assemtools3.read_top_blat_line(myData['blatOutFile'])
    if blatLine['strand'] == '-':
        print('Need to make reverse complement of sequence!')
        myData['assemFa'] = myData['originalAssem'] + '.rc.fa'
        fastaSeqs = assemtools3.read_fasta_file_to_dict(myData['originalAssem'])
    #    print(fastaSeqs)
        name = list(fastaSeqs.keys())[0]
        seq = fastaSeqs[name]['seq']
        seq = revcomp(seq)
        seq = add_breaks_to_line(seq,n=100)
        print('seq name is',name)
        outFile = open(myData['assemFa'],'w')
        outFile.write('>%s\n' % (name))
        outFile.write(seq)
        outFile.close()
    
        # then, need to redo the blat
        myData['blatOutFile'] = myData['assemFa'] + '.blatTMP'
        cmd = 'blat %s %s %s' % (myData['assemFa'],myData['targetFa'],myData['blatOutFile'])
        print(cmd)
        runCMD(cmd)
        blatLine = read_top_blat_line(myData['blatOutFile'])

    # print out stats of where the query is 
    print('query size',blatLine['qSize'])
    print('target size',blatLine['tSize'])
    print('q hit:',blatLine['qStart'],'-',blatLine['qEnd'])
    print('t hit:',blatLine['tStart'],'-',blatLine['tEnd'])

    # check the size hit
    if(blatLine['qEnd'] - blatLine['qStart']) != blatLine['qSize']:
        print('Did not account for whole length in single hit.  What to do???')
        sys.exit()

    # for now, will assume that the vector sequence is not across the circle junction -- that would be a complication
    # will do the coordinates in 1 base system, not bed/psl
    tStart = blatLine['tStart'] + 1
    fastaSeqs = assemtools3.read_fasta_file_to_dict(myData['assemFa'])
    name = list(fastaSeqs.keys())[0]
    newName = name.split()[0] # simplify the name

    originalSeq = fastaSeqs[name]['seq']
    print(tStart)
    part1 = originalSeq[tStart-1:]  # starts at the original seq
    part2 = originalSeq[0:tStart-1]
    print('p1 len',len(part1))
    print('p2 len',len(part2))
    print('total',len(part1) + len(part2))

    outFile = open(myData['rotatedFa'],'w')
    outFile.write('>%s\n' % (newName))
    seq = part1 + part2
    seq = assemtools3.add_breaks_to_line(seq,n=100)
    outFile.write(seq)
    outFile.close()

    myData['rotatedFaBlatOutFile'] = myData['rotatedFa']+ '.blat'
    cmd = 'blat %s %s %s' % (myData['rotatedFa'],myData['targetFa'],myData['rotatedFaBlatOutFile'])
    print(cmd)
    assemtools3.runCMD(cmd)
    blatLine = assemtools3.read_top_blat_line(myData['rotatedFaBlatOutFile'])

    if myData['doClean'] is True:
        print('making version without the vector!')
        print('query size',blatLine['qSize'])
        print('target size',blatLine['tSize'])
        print('q hit:',blatLine['qStart'],'-',blatLine['qEnd'])
        print('t hit:',blatLine['tStart'],'-',blatLine['tEnd'])
    
        fastaSeqs = assemtools3.read_fasta_file_to_dict(myData['rotatedFa'])
        name = list(fastaSeqs.keys())[0]
        seq = fastaSeqs[name]['seq']
    
        #tEnd is in bedFormat -- so can use directly, will get the next base...
        newSeq = seq[blatLine['tEnd']:]
        newSeq = assemtools3.add_breaks_to_line(newSeq,n=100)
        outFile = open(myData['rotatedCleanFa'],'w')
        outFile.write('>%s\n' % (name))
        outFile.write(newSeq)
        outFile.close()
#####################################################################


