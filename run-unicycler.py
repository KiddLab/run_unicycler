import sys
import assemtools3
import argparse
import os


parser = argparse.ArgumentParser(description='run Unicycler program with Illumina and Nanopore data')


parser.add_argument('--outDirBase', type=str,help='Path to output directory',required=True)
parser.add_argument('--name', type=str,help='name of project',required=True)
parser.add_argument('--fqR1', type=str,help='Illumina Read 1 File',required=True)
parser.add_argument('--fqR2', type=str,help='Illumina Read 2 File',required=True)
parser.add_argument('--contam', type=str,help='fasta of contamination (E. coli), with bwa mem index',required=True)



args = parser.parse_args()

# setup dictionary for holding information and passing to functions
myData = {}
myData['outDirBase'] = args.outDirBase 
myData['name'] = args.name 
myData['fqR1'] = args.fqR1 
myData['fqR2'] = args.fqR2 
myData['contam'] = args.contam 


# setup needed files
if myData['outDirBase'][-1] != '/':
    myData['outDirBase'] += '/'

if os.path.isdir(myData['outDirBase']) is False:
    print('Output dir doest not exist, making it!')
    cmd = 'mkdir ' + myData['outDirBase']
    print(cmd)
    assemtools3.runCMD(cmd)

assemtools3.check_prog_paths(myData)
###############################################################################

# step 1 run cutadpt to remove adapter from Illumina reads
assemtools3.run_cutadapt(myData)

# step 2, remove reads that map to E. coli
assemtools3.filter_contam_illumina(myData)




