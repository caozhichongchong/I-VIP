import os
from Bio import SeqIO
import argparse
import glob
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-fasta",
                    help="format of fasta sequences", type=str, default='.fa',metavar='.fa')
parser.add_argument("-orf",
                    help="format of ORFs", type=str, default='.faa',metavar='.faa')
parser.add_argument('--gbff',
                    default="gbff", action='store', type=str, metavar='gbff',
                    help="Set the directory of Genbank files (default: gbff)\nNo gbff input: None")

parser.add_argument('--orftype',
                    help="The orf type, eg: 0 for genbank parsing; 1 for prodigal prediction",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=0, type=int)
parser.add_argument('--distance',
                    default=4000, action='store', type=int, metavar='4000',
                    help='Distance cutoff for two cassettes to be clustered together (default is 4000)')
parser.add_argument('--cutoff',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='E-value cutoff for attC site identification (default is 1.0)')
parser.add_argument('--module',
                    help="Alternative attC search strategies for Module A \
                    (0: global search by Module A1; 1: local search by Module A2), (default 0)",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=0, type=int)
parser.add_argument('--quick',
                    help="attC search and filtering setting \
                    (from least strict to most strict: \'max\', \'nohmm\', \'mid\', \'\', \'rfam\'), (default \'mid\')",
                    metavar="mid",
                    action='store', default='mid', type=str)
parser.add_argument('--usearch',
                    help="Use usearch ahead of blastp for integrase and sulI alignment,"+
                         " \'None\' for not using usearch, \'usearch\' or other version of usearch"+
                         "for using usearch (default: \'None\')",
                    metavar="version of usearch",
                    action='store', default='usearchv8', type=str)
parser.add_argument('--thread',
                    help="The thread number assigned for running I-VIP (default 15)",
                    metavar="15", action='store', default=15, type=int)

################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.fasta)
in_dir, input_file = os.path.split(input_path)
in_dir = os.path.abspath(in_dir)
list_data = glob.glob(os.path.join(in_dir,'*'+args.fasta))
search_path = 'output'
try:
    os.mkdir(search_path)
except OSError:
    pass
try:
    os.mkdir('Temp')
except OSError:
    pass
#IntI and SulI length
Length = dict()
for line in open(os.path.join('database/' + 'Integrase.fasta.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
for line in open(os.path.join('database/' + 'sul1_sequences_in_SARG.txt.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
################################################### Function ########################################################
def Calculate_length(list_of_files):
    AA_seq=dict()
    for file_name in list_of_files:
        try:
            Fasta_name = open(file_name, 'r')
            filedir, file_name = os.path.split(file_name)
            f = open(file_name + '.length', 'w')
            for record in SeqIO.parse(Fasta_name, 'fasta'):
                f.write(str(record.id) + '\t' + str(record.description).replace('\t', ' ') + '\t' + str(
                    len(record.seq)) + '\n')
                AA_seq.setdefault(str(record.id),str(record.seq))
            f.close()
        except (IOError):
            print 'Files were missing: '+'file_name'
        finally:
            print 'Cleaning Up...'
            del file_name
    return AA_seq


def Cmsearch(list_file):
    print 'Cmsearch'
    for file in list_file:
        Newlist = []
        Total_length=0
        for record in SeqIO.parse(file, 'fasta'):
            Total_length+=len(record.seq)
        #for large fasta, first split into 2Mb files for cmsearch and then merge search results
        if Total_length>=2097152:
            Length=0
            for record in SeqIO.parse(file, 'fasta'):
                Length+=len(record.seq)
                f1=open(str(file)+'_'+str(int(Length/2097152)),'ab')
                f1.write('>'+str(record.id)+'\t'+str(record.description)+'\n')
                f1.write(str(record.seq)+'\n')
                f1.close()
                if str(file)+'_'+str(int(Length/2097152)) not in Newlist:
                    Newlist.append(str(file)+'_'+str(int(Length/2097152)))
        else:
            Newlist.append(file)
        Result=[]
        for files in Newlist:
            in_dir,files=os.path.split(files)
            cmd= 'cmsearch --tblout '+os.path.join(str(in_dir)+'/'+str(search_path),str(files) + '.Z.max.attc.hits.txt')+' -Z `du -m '+str(files) + \
           ' | cut -f1` --'+str(args.quick)+' --cpu '+str(args.thread)+' -E '+str(args.cutoff)+' database/attc_4.cm.hmm '+str(files)
            print cmd
            os.system(cmd)
            Result.append(os.path.join(str(search_path),str(files) + '.Z.max.attc.hits.txt'))
        if file not in Newlist:
            Newfilelist=''
            for newfiles in Newlist:
                Newfilelist+=' '+newfiles
            Resultlist=''
            for files in Result:
                Resultlist+=' '+files
            in_dir, file = os.path.split(file)
            os.system('cat '+str(Resultlist)+' > '+os.path.join(str(search_path),str(file) + '.Z.max.attc.hits.txt\n'))
            os.system('rm -rf ' + str(Resultlist) + ' \n')
            os.system('rm -rf ' + str(Newfilelist) + ' \n')

def filter_blast_list(file, Cutoff_identity, Cutoff_hitlength):
    try:
        f = open(file + '.filter.txt', 'w')
        for line in open(file, 'rb'):
            if float(str(line).split('\t')[2]) >= Cutoff_identity:
                if float(str(line).split('\t')[3]) >= Cutoff_hitlength * float(
                        Length.get(str(line).split('\t')[1])) / 100:
                    f.write(str(line))
        f.close()
        return file + '.filter.txt'
    except IOError:
        print str(file) + ' missing!'


def Extractaa(AA_seq,root,searchfile,gene):
    f1=open(os.path.join(root,'all.'+str(gene)+'.aa'),'wb')
    for line in open(searchfile,'rb'):
        try:
            f1.write('>'+str(line).split('\t')[0].split(' ')[0]+'\n'+str(AA_seq[str(line).split('\t')[0].split(' ')[0]])+'\n')
        except KeyError:
            print 'AA not found for '+str(line).split('\t')[0].split(' ')[0]
    f1.close()


################################################### Programme #######################################################
print 'Requirement: database folder containing all attC, integrase and sulI databases'
print 'Requirement: scripts folder containing all scripts'
print 'Requirement: both sequences and orfs of input data'
print 'Attention: It\'s recommended to remove "_" and ":" from input data name'
#Step1 format orfs
cmd1='python scripts/Addfilename.py -orf '+str(args.orf)+'\ncat Temp/*.add > Temp/all.orfs.aa\n'
print cmd1
os.system(cmd1)
print 'Step1 Finished: format orfs'
#Step2 search for integrase and sulI
if args.usearch !='None':
    print 'Start search Integrase and sulI genes by usearch!'
    cmd2p = str(args.usearch)+" -ublast Temp/all.orfs.aa " + \
            "-db database/Integrase.fasta -evalue 1e-0 -accel 0.5 -blast6out " \
            + "Temp/all.orfs.aa.Int.usearch.txt -threads "+str(args.thread)+"\n"
    print cmd2p
    os.system(cmd2p)
    ORF_seq = Calculate_length(['Temp/all.orfs.aa'])
    Extractaa(ORF_seq,'Temp','Temp/all.orfs.aa.Int.usearch.txt','Int')
    cmd2p = str(args.usearch)+" -ublast Temp/all.orfs.aa " + \
            "-db database/sul1_sequences_in_SARG.txt -evalue 1e-0 -accel 0.5 -blast6out " \
            + "Temp/all.orfs.aa.sul1.usearch.txt -threads "+str(args.thread)+"\n"
    print cmd2p
    os.system(cmd2p)
    Extractaa(ORF_seq, 'Temp', 'Temp/all.orfs.aa.sul1.usearch.txt', 'sul1')
    print 'Start search Integrase and sulI genes by blastp!'
    Int_file='Temp/all.Int.aa'
    SulI_file='Temp/all.sul1.aa'
else:
    SulI_file=Int_file = 'Temp/all.orfs.aa'
cmd2p = "blastp -query "+str(Int_file)+" -db database/Integrase.fasta -out " +str(search_path)\
            + "/all.Int.blast.txt  -outfmt 6 -max_target_seqs 1 -evalue 1e-2 -num_threads "+\
        str(args.thread)+" \n"
print cmd2p
os.system(cmd2p)
cmd2p = "blastp -query "+str(SulI_file)+" -db database/sul1_sequences_in_SARG.txt -out " +str(search_path)\
            + "/all.sul1.blast.txt  -outfmt 6 -max_target_seqs 1 -evalue 1e-2 -num_threads "+\
        str(args.thread)+" \n"
print cmd2p
os.system(cmd2p)
fsul1=open(os.path.join(search_path, 'all.sul1.blast.txt'))
for line in fsul1:
    Filename =str(line).split('\t')[0].split(args.orf)[0].split(str(args.fasta))[0]
    fnew=open(os.path.join(search_path, Filename+str(args.fasta)+'.sul1.txt'),'ab')
    fnew.write(line)
filter_blast_list(str(search_path)+'/all.Int.blast.txt', 80, 50)
fint1=open(os.path.join(search_path, 'all.Int.blast.txt.filter.txt'))
for line in fint1:
    Filename = str(line).split('\t')[0].split(args.orf)[0].split(str(args.fasta))[0]
    fnew=open(os.path.join(search_path, Filename+str(args.fasta)+'.Int.txt'),'ab')
    fnew.write(line)
print 'Step2 Finished: search for integrase and sulI'
#Step3 attC search
if args.module==0:
    Cmsearch(list_data)
    cmd4 = 'python scripts/attCformat.py -input '+str(search_path)+'\n'
    print cmd4
    os.system(cmd4)
    pass
else:
    cmd3='python scripts/Integrase_extraction_seq.py --fasta '+str(args.fasta)+' --orf '+str(args.orf)+\
         ' --AA_type '+str(args.orftype)+' --thread '+str(args.thread)+' --cutoff '+str(args.cutoff)+\
        ' --quick '+str(args.quick)
    print cmd3
    os.system(cmd3)
print 'Step3 Finished: search for attC'
#Step4 Integron combination
cmd4='python scripts/Result_Combine_seq.py --fasta '+str(args.fasta)+' --distance '+str(args.distance)+\
     ' --AA_type '+str(args.orftype)+' --attc_evalue '+str(args.cutoff) +' --orf '+ str(args.orf) +'\n'
print cmd4
os.system(cmd4)
print 'Step4 Finished: integron identification and classification'
#Step5 Integron extraction
cmd5='python scripts/Integron_extraction_seq.py  --AA_format '+str(args.orftype)+' --gbffdir '+\
    str(args.gbff) +' --fasta '+ str(args.fasta) +' --orf '+ str(args.orf) +'\n'
os.system(cmd5)
print cmd5
print 'Step5 Finished: integron extraction'




