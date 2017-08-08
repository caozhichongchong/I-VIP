import os
from Bio import SeqIO
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-input",
                    help="The input folder", type=str, default='extract',metavar='extract')
parser.add_argument("-target",
                    help="The taget files of attC", type=str, default='attc.hits.txt',metavar='attc.hits.txt')
################################################## Definition ########################################################
args = parser.parse_args()
in_dir= os.path.abspath(args.input)
attc_filelist=glob.glob(os.path.join(in_dir,'*'+str(args.target)))
################################################### Programme #######################################################
for attcfile in attc_filelist:
    f1=open(attcfile+'2.txt','ab')
    for line in open(attcfile,'rb'):
        if str(line)[0]=='#':
            pass
        else:
            line=str(line).replace(' # ','#')
            while line != str(line).replace('  ',' '):
                line=str(line).replace('  ',' ')
            line = str(line).replace(' ', '\t')
            line = str(line).replace('#',' # ')
            f1.write(line)
    f1.close()
