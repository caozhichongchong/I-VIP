import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument("-f",
                    help="file type or filename extension of your sequences",
                    type=str, default='.fa',metavar='.fa or .fasta or .fna')
parser.add_argument("-input",
                    help="The input folder of attC search results", 
                    type=str, default='extract',metavar='extract')
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result',metavar='I-VIP_result')
parser.add_argument('--q',
                    help="Optional: attC search setting \
                    (from least strict to most strict: \'max\', \'nohmm\', \'mid\', \'\', \'rfam\'),\
                     (default \'mid\')",
                    metavar="mid",
                    action='store', default='mid', type=str)
parser.add_argument('--c',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='Optional: set the e-value cutoff for attC search (default is 1.0)')
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running I-VIP (default 15)",
                    metavar="15", action='store', default=15, type=int)
parser.add_argument('--m',
                    help="Optional: set the search strategies for Module A \
                    (1: global search by Module A1; 2: local search by Module A2), \
                    (default \'1\' for global search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=2, type=int)
parser.add_argument('--fn',
                    help="Optional: set the format strategies for attC search results \
                    (1: for Module A1; 2: for Module A2), \
                    (default \'1\' for Module A1)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=2, type=int)
parser.add_argument('--cmsearch',
                    help="Optional: complete path to cmsearch if not in PATH,",
                    metavar="/usr/local/bin/cmsearch",
                    action='store', default='cmsearch', type=str)


################################################### Function ########################################################


def moduleA1(list_file):
    # Module A1 global search for attC
    for file in list_file:
        Newlist = []
        # calculate the total data size (bp) for query data
        Total_length = 0
        for record in SeqIO.parse(file, 'fasta'):
            Total_length += len(record.seq)
        # large data is firstly split into 2Mb files for attC search
        # for parallel comparison of cmserach results (e-value) 
        if Total_length >= 2000000:
            Length = 0
            for record in SeqIO.parse(file, 'fasta'):
                Length += len(record.seq)
                f1 = open(str(file) + '_' + str(int(Length / 2000000)), 'a')
                f1.write('>' + str(record.id) + '\t' + str(record.description) + '\n')
                f1.write(str(record.seq) + '\n')
                f1.close()
                if str(file) + '_' + str(int(Length / 2000000)) not in Newlist:
                    Newlist.append(str(file) + '_' + str(int(Length / 2000000)))
        else:
            Newlist.append(file)
        # attC search by
        Result = []
        for files in Newlist:
            filesize = max(int(os.path.getsize(files)) / 1000000,1)
            in_dir, files = os.path.split(files)
            cmd = args.cmsearch + ' --tblout ' + os.path.join(str(search_path), str(
                files) + '.Z.max.attc.hits.txt') + ' -Z %s --' % (filesize) \
                  + str(args.q) + ' --cpu ' + str(args.t) + ' -E ' \
                  + str(args.c) + ' database/attc_4.cm.hmm ' + os.path.join(in_dir, str(files))
            print(filesize, cmd)
            os.system(cmd)
            Result.append(os.path.join(str(search_path), str(files) + '.Z.max.attc.hits.txt'))
        # merge search results of split files
        if file not in Newlist:
            Newfilelist = ''
            for newfiles in Newlist:
                Newfilelist += ' ' + newfiles
            Resultlist = ''
            for files in Result:
                Resultlist += ' ' + files
            in_dir, file = os.path.split(file)
            os.system(
                'cat ' + str(Resultlist) + ' > ' + os.path.join(str(search_path),
                                                                str(file) + '.Z.max.attc.hits.txt \n'))
            os.system('rm -rf ' + str(Resultlist) + ' \n')
            os.system('rm -rf ' + str(Newfilelist) + ' \n')


################################################## Definition ########################################################
args = parser.parse_args()
in_dir= os.path.abspath(args.input)
# record Module A1 process
flog = open('ModuleA1.log', 'w')
search_path = args.r + '/output'


################################################### Programme #######################################################
# moduleA1
if args.m==1:
    Data_dir = os.path.abspath(args.i)
    list_data = glob.glob(os.path.join(Data_dir, '*' + args.f))
    moduleA1(list_data)
# load attC search result files
attc_filelist=glob.glob(os.path.join(in_dir,'*attc.hits.txt'))
# format attC results
for attcfile in attc_filelist:
    f1=open(attcfile+'2.txt','a')
    for line in open(attcfile,'r'):
        if str(line)[0]=='#':
            pass
        else:
            line=str(line).replace(' # ','#')
            while line != str(line).replace('  ',' '):
                line=str(line).replace('  ',' ')
            line = str(line).replace(' ', '\t')
            line = str(line).replace('#',' # ')
            filedir, filename = os.path.split(attcfile)
            filename = filename.split('.Z.max.attc.hits.txt')[0]
            if args.fn == 1:
                f1.write(filename+'_'+line)
            else:
                f1.write(line)
    f1.close()
print 'Finished ModuleA1!'
flog.write('Finished ModuleA1!\n')
flog.close()