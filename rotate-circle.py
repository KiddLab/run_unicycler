import sys
import assemtools3
import argparse



parser = argparse.ArgumentParser(description='Post process circular alignment')


parser.add_argument('--assem', type=str,help='assembly fasta file',required=True)
parser.add_argument('--target',type=str,help='fasta file of sequnece to use for begining of rotated sequence',required=True)
parser.add_argument('--clean', help='whether or not to make assembly file with the target sequence removed, default is TRUE',
                               action='store_true')



args = parser.parse_args()

# setup dictionary for holding information
myData = {}
myData['assemFa']  = args.assem # assembly output
myData['originalAssem'] = args.assem 
myData['targetFa']  = args.target 
myData['doClean'] = args.clean

myData['rotatedFa'] = myData['originalAssem'] + '.rotate.fa'
myData['rotatedCleanFa'] = myData['originalAssem'] + '.rotate.clean.fa'



print('initial assembly:',myData['assemFa'])
print('portion to set at start:',myData['targetFa'])
print('make version without the target:',myData['doClean'])


assemtools3.check_prog_paths(myData)
###############################################################################

# run blat of target vs assembly, to determine what the correct orientation is
myData['blatOutFile'] = myData['assemFa'] + '.blatTMP'
cmd = 'blat %s %s %s' % (myData['assemFa'],myData['targetFa'],myData['blatOutFile'])
print(cmd)
assemtools3.runCMD(cmd)
blatLine = assemtools3.read_top_blat_line(myData['blatOutFile'])
if blatLine['strand'] == '-':
    print('Need to make reverse complement of sequence!')
    myData['assemFa'] = myData['originalAssem'] + '.rc.fa'
    fastaSeqs = assemtools3.read_fasta_file_to_dict(myData['originalAssem'])
#    print(fastaSeqs)
    name = list(fastaSeqs.keys())[0]
    seq = fastaSeqs[name]['seq']
    seq = assemtools3.revcomp(seq)
    seq = assemtools3.add_breaks_to_line(seq,n=100)
    print('seq name is',name)
    outFile = open(myData['assemFa'],'w')
    outFile.write('>%s\n' % (name))
    outFile.write(seq)
    outFile.close()
    
    # then, need to redo the blat
    myData['blatOutFile'] = myData['assemFa'] + '.blatTMP'
    cmd = 'blat %s %s %s' % (myData['assemFa'],myData['targetFa'],myData['blatOutFile'])
    print(cmd)
    assemtools3.runCMD(cmd)
    blatLine = assemtools3.read_top_blat_line(myData['blatOutFile'])

###############################################################################

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
    
    
    
    
    
    

    
    
    
    
    














