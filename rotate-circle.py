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
print('make version without the target:',myData['doClean'], flush=True)


#assemtools3.check_prog_paths(myData)
###############################################################################

assemtools3.do_rotate_circle(myData)    
    
    
    
    
    

    
    
    
    
    














