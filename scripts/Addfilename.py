import os
from Bio import SeqIO
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-input",
                    help="input path", type=str, default='.',metavar='.')
parser.add_argument("-orf",
                    help="format of ORFs", type=str, default='.faa',metavar='.faa')
################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.input)
list_fasta = glob.glob(os.path.join(input_path,'*'+args.orf))

################################################### Function ########################################################
def addname(file_name):
    Fasta_name = open(file_name, 'r')
    f = open('Temp/'+file_name + '.add', 'w')
    in_dir, input_file = os.path.split(file_name)
    for record in SeqIO.parse(Fasta_name, 'fasta'):
        f.write('>'+str(input_file)+'_'+str(record.id) + '\t' + str(record.description).replace('\t', ' ') + '\n' + str(
            str(record.seq)) + '\n')
    f.close()

################################################### Programme #######################################################
for file_name in list_fasta:
    filedir,file_name=os.path.split(file_name)
    print file_name
    addname(file_name)