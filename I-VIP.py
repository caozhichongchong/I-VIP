import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import glob
import copy


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input',metavar='input')
parser.add_argument("-f",
                    help="file type or filename extension of your sequences.\n \
                         To input genbank file, set \"-f .gbff\" or \"-f .gbff.gz\"",
                    type=str, default='.fa',metavar='.fa, .fasta, .fna, .gbff, .gbff.gz')
# optional files input and output
parser.add_argument("--o",
                    help="Optional: to provide your own CDS or ORFs files, \
                    please input the file type or filename extension of your CDS or ORFs",
                    type=str, default='.faa',metavar='.faa')
parser.add_argument('--ot',
                    help="Optional: to provide your own CDS or ORFs files, \
                     please input the orf type or the method you used to extract orfs,\
                     eg: 1 for genbank parsing; 2 for prodigal prediction",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=2, type=int)
parser.add_argument('--g',
                    default="None", action='store', type=str, metavar='gbff or None',
                    help="Optional: to provide your own genbank files for gene cassettes annotation,\
                     please set the directory of Genbank files \
                    (default: \'None\' for No gbff input)")
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result',metavar='I-VIP_result')
# optional search parameters
parser.add_argument('--d',
                    default=4000, action='store', type=int, metavar='4000',
                    help='Optional: set the distance cutoff for two cassettes to be clustered together\
                     (default is 4000)')
parser.add_argument('--c',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='Optional: set the e-value cutoff for attC search (default is 1.0)')
parser.add_argument('--m',
                    help="Optional: set the search strategies for Module A \
                    (1: global search by Module A1; 2: local search by Module A2), \
                    (default \'2\' for local search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=2, type=int)
parser.add_argument('--q',
                    help="Optional: attC search setting \
                    (from least strict to most strict: \'max\', \'nohmm\', \'mid\', \'\', \'rfam\'),\
                     (default \'mid\')",
                    metavar="mid",
                    action='store', default='mid', type=str)
parser.add_argument('--a',
                    help="Optional: annotate gene cassettes against antibiotic and metal resistance databases",
                    metavar="Y or N", action='store', default='Y', type=str, choices=["Y", "N"],)
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running I-VIP (default 1)",
                    metavar="1 or more", action='store', default=1, type=int)
parser.add_argument("--tx",
                    help="Optional: a file of taxonomy metadata (under your input folder)",
                    type=str, default='None',
                    metavar='taxon.txt or None')
parser.add_argument("--tc",
                    help="Optional: column number corresponding to the taxonomy, i.e., from phylum to strain",
                    type=str, default='None',
                    metavar='start column number, end column number, i.e., \'2,8\' or \'None\'')
# requirement for software calling
parser.add_argument('--u',
                    help="Optional: use two-step method for integrase and sul1 search,"+
                         " \'None\' for using one step, \'usearch\' or \'diamond\' for using two-step \
                         (complete path to usearch or diamond if not in PATH, \
                         please make sure the search tools can be directly called), (default: \'None\')",
                    metavar="None or usearch",
                    action='store', default='None', type=str)
parser.add_argument('--cmsearch',
                    help="Optional: complete path to cmsearch if not in PATH,",
                    metavar="/usr/local/bin/cmsearch",
                    action='store', default='cmsearch', type=str)
parser.add_argument('--prodigal',
                    help="Optional: complete path to prodigal if not in PATH,",
                    metavar="/usr/local/bin/prodigal",
                    action='store', default='prodigal', type=str)
parser.add_argument('--blastp',
                    help="Optional: complete path to blastp if not in PATH,",
                    metavar="/usr/local/bin/blastp",
                    action='store', default='blastp', type=str)


################################################## Definition ########################################################
args = parser.parse_args()
in_dir= os.path.abspath(args.i)
if args.tx == 'none':
    txset = 'None'
elif args.i == 'example':
    txset = os.path.join(in_dir + '/example',args.tx)
else:
    txset = args.tx
if args.tc == 'none':
    tcset = 'None'
else:
    tcset = args.tc

workingdir = os.path.abspath(os.path.dirname(__file__))

# set up input files
#print [args.f, any (format in args.f for format in ['gbff','gbff.gz'])]
if any (format in args.f for format in ['gbff','gbff.gz']):
    os.system('python '+workingdir+'/scripts/GbffParser.py -i ' + str(in_dir) + ' -f ' + str(args.f) + ' \n')
    print 'python '+workingdir+'/scripts/GbffParser.py -i ' + str(in_dir) + ' -f ' + str(args.f) + ' \n'
    input_extension='.fa'
    orf_extension='.faa'
    genbank = in_dir
    orf_format = 1
else:
    input_extension = args.f
    orf_extension = args.o
    genbank = args.g
    orf_format = args.ot
# fetch all input sequences
list_data = glob.glob(os.path.join(in_dir,'*'+input_extension))
search_path = args.r + '/output'
try:
    os.mkdir(args.r)
except OSError:
    pass
try:
    os.mkdir(search_path)
except OSError:
    pass
try:
    os.mkdir(args.r + '/Temp')
except OSError:
    pass
# load the integrase and sul1 database
Length = dict()
for line in open(os.path.join(workingdir+'/database/' + 'Integrase.fasta.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
for line in open(os.path.join(workingdir+'/database/' + 'sul1_database.txt.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
# record I-VIP process
flog = open('I-VIP.log', 'wb')


################################################### Programme #######################################################
print "\
------------------------------------------------------------------------\n\
I-VIP (Integron Identification and Visualization)\n\
Requirement: cmsearch, blast\n\
Optional: prodigal, usearch, diamond, and cytoscape\n\
Please download the whole package I-VIP.py, scripts and database\n\
Please remove  \":\", \"____\", and \"#\" from your filename\n\
It\'s highly recommended to keep less than 10,000 files in your input folder\n\
Copyright:An-Ni Zhang, Prof. Tong Zhang, University of Hong Kong\n\
Citation:\n\
2. Cury J, Jove T, Touchon M, Neron B, Rocha EP: Identification and analysis of integrons and \
cassette arrays in bacterial genomes. Nucleic acids research 2016, 44:4539-4550. (attC database)\n\
3. Yang Y, Jiang X, Chai B, Ma L, Li B, Zhang A, Cole JR, Tiedje JM, Zhang T: ARGs-OAP: online analysis\
pipeline for antibiotic resistance genes detection from metagenomic data using an integrated \
structured ARG-database. Bioinformatics 2016. (optional: antibiotic resistance database)\n\
4. Li L-G, Xia Y, Zhang T: Co-occurrence of antibiotic and metal resistance genes revealed in \
complete genome collection. The ISME Journal 2016.(optional: metal resistance database)\n\
Contact anniz44@mit.edu and/or caozhichongchong@gmail.com\n\
------------------------------------------------------------------------\n\
\nStart I-VIP log (in I-VIP.log)\n\
"


if list_data == []:
    print 'Find no ' + str(input_extension) + ' files in ' + str(in_dir) + '!\n'
    print 'Stop running I-VIP!\n'
    flog.write('Find no ' + str(input_extension) + ' files in ' + str(in_dir) + '!\nStop running I-VIP!\n')
else:
    print 'I-VIP: Searching integrons in ' + str(len(list_data)) + ' files of ' + str(in_dir) + '\n'
    if len(list_data) >= 10000:
        print 'It\'s highly recommended to keep less than 10,000 files in your input folder!\n'
    # Step1 format orfs and search for integrase and sul1
    cmd1 = 'python '+workingdir+'/scripts/Preparation.py -i ' + str(in_dir) + ' -f ' + str(input_extension) + ' --o '+str(orf_extension) + \
            ' --r ' + str(args.r) + ' --t ' + str(args.t) + ' --u ' + str(args.u) + ' --ot ' + str(orf_format) + \
           ' --prodigal '+ str(args.prodigal) + ' --blastp ' + str(args.blastp) +' \n'
    if args.u != 'None':
        if 'usearch' in args.u:
            print 'Start search integrase and sul1 genes by usearch!\n'
            flog.write('Start search integrase and sul1 genes by usearch!\n')
        elif "diamond" in args.u:
            print 'Start search integrase and sul1 genes by diamond!\n'
            flog.write('Start search integrase and sul1 genes by diamond!\n')
        else:
            print "Wrong input for -u (usearch, diamond or None), proceed one-step search using blastp\n"
            cmd1 = 'python '+workingdir+'/scripts/Preparation.py -i ' + str(in_dir) + ' -f ' + str(input_extension) + ' --o ' + str(
                input_extension) + ' --t ' + str(args.t) + str(args.u) + ' --ot ' + str(orf_format) + ' --u None --prodigal ' + str(
                args.prodigal) + ' --r ' + str(args.r) +' --blastp ' + str(args.blastp)
            flog.write("Wrong input for -u (usearch, diamond or None), proceed one-step search using blastp\n")
    else:
        print 'Start search integrase and sul1 genes by blastp!\n'
        flog.write('Start search integrase and sul1 genes by blastp!\n')
    print cmd1
    flog.write(cmd1)
    os.system(cmd1)
    print 'Step1 Finished: format orfs + search for integrase and sul1\n'
    flog.write('Step1 Finished: format orfs + search for integrase and sul1\n')
    # Step2 attC search
    if args.m == 1:
        cmd2 = 'python '+workingdir+'/scripts/ModuleA1.py -i ' + str(in_dir) + ' -input ' + str(search_path) + ' -f ' + str(input_extension) + \
               ' --m ' + str(args.m) + ' --q ' + str(args.q) + ' --c ' + str(args.c) + ' --t ' + str(args.t) \
               + ' --r ' + str(args.r) + ' --cmsearch ' + str(args.cmsearch) + '\n'
    else:
        cmd2 = 'python '+workingdir+'/scripts/ModuleA2.py -input ' + str(in_dir) + ' --f ' + str(input_extension) + ' --o ' + str(orf_extension) + \
               ' --ot ' + str(args.r) + '/Temp/all.orf.length --t ' + str(args.t) + ' --c ' + str(args.c) + \
               ' --r ' + str(args.r) + ' --q ' + str(args.q) + ' --cmsearch ' + str(args.cmsearch) + '\n'
    print cmd2
    flog.write(cmd2)
    os.system(cmd2)
    print 'Step2 Finished: search for attC\n'
    flog.write('Step2 Finished: search for attC\n')
    # Step3 Integron combination
    cmd3='python '+workingdir+'/scripts/Integron_identification.py -i ' + str(in_dir) + ' -f '+str(input_extension)+' -o '+str(orf_extension)+\
         ' --d '+str(args.d)+ ' --r ' + str(args.r) +\
         ' --m ' + str(args.m) + ' --ot ' + str(args.r) + '/Temp/all.orf.length --c '+str(args.c) + '\n'
    print cmd3
    flog.write(cmd3)
    os.system(cmd3)
    print 'Step3 Finished: integron identification and classification\n'
    flog.write('Step3 Finished: integron identification and classification\n')
    # Step4 Integron extraction
    cmd4='python '+workingdir+'/scripts/Integron_extraction.py -i ' + str(in_dir) + ' --g '+\
        str(genbank) + ' --r ' + str(args.r) + ' --f ' + str(input_extension) + ' --o ' + str(orf_extension) + '\n'
    print cmd4
    flog.write(cmd4)
    os.system(cmd4)
    print 'Step4 Finished: integron extraction\n'
    flog.write('Step4 Finished: integron extraction\n')
    # Step5 Gene cassettes annotation
    if txset != 'None':
        # normalize the taxonomy information
        cmd5 = 'python '+workingdir+'/scripts/Taxon_normalization.py -i ' + str(in_dir) + ' --tc ' + \
               str(tcset) + ' --tx ' + str(txset) + '\n'
    else:
        cmd5 = ''
    cmd5 += 'python '+workingdir+'/scripts/Integron_annotation.py -i ' + str(in_dir) + ' -f ' + \
               str(args.f) + ' --r ' + str(args.r) + ' --a ' + str(args.a) \
            + ' --t ' + str(args.t) + ' --tx ' + str(txset) + '.norm --blastp ' + str(args.blastp) +'\n'
    print cmd5
    flog.write(cmd5)
    os.system(cmd5)
    print 'Step5 Finished: gene cassettes annotation\n'
    flog.write('Step5 Finished: gene cassettes annotation\n')
    os.system('rm -rf Temp \n')
    # Step6 Integron visualization
    cmd6 = 'python '+workingdir+'/scripts/Visualization.py -i ' + str(in_dir) + ' --tc ' + \
               str(tcset) + ' --tx ' + str(txset) + '.norm --r ' + \
           str(args.r) + ' --ot ' + str(args.r) + '/Temp/all.orf.length \n'
    print cmd6
    flog.write(cmd6)
    os.system(cmd6)
    print 'Step6 Finished: integron visualization\n'
    flog.write('Step6 Finished: integron visualization\n')
    #os.system('rm -rf ' + str(args.r) +  '/Temp \n')
flog.close()
