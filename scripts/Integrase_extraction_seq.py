import os
from Bio import SeqIO
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i',
                    help="The irratation time",
                    metavar="1", action='store', default=1, type=int)
parser.add_argument("-input",
                    help="The input integron locus file", type=str, default='output/all.Int.blast.txt',metavar='output/all.Int.blast.txt')
parser.add_argument("-target",
                    help="The taget sequences file", type=str, default='None',metavar='None')
parser.add_argument("--fasta",
                    help="format of fasta sequences", type=str, default='.fa',metavar='.fa')
parser.add_argument("--orf",
                    help="format of ORFs", type=str, default='.faa',metavar='.faa')
parser.add_argument("--func_type",
                    help="Function, eg: 0 for sequence extraction; 1 for cmsearch result analysis; 2 for attC result formating",
                    metavar="0 or 1 or 2",
                    choices=[0, 1, 2],
                    action='store', default=0, type=int)
parser.add_argument('--AA_type',
                    help="The AA type, eg: 0 for genbank parsing; 1 for prodigal prediction",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=0, type=int)
parser.add_argument('--thread',
                    help="The thread number assigned for running I-VIP",
                    metavar="15", action='store', default=15, type=int)
parser.add_argument('--cutoff',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='E-value cutoff for attC site identification (default is 1.0)')
parser.add_argument('--outputdir',
                    default="output", action='store', type=str, metavar='output',
                    help="Set the output directory for integron identification (default: output)")
parser.add_argument('--head',
                    help="Extract the extremity of sequence",
                    type=str, default='Warning.txt',metavar='Warning.txt')
parser.add_argument('--quick',
                    help="attC search and filtering setting \
                    (from least strict to most strict: \'max\', \'nohmm\', \'mid\', \'\', \'rfam\'), (default \'mid\')",
                    metavar="mid",
                    action='store', default='mid', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.fasta)
in_dir, input_file = os.path.split(input_path)
Int_file = (os.path.join(in_dir, args.input))
if args.func_type == 0:
    if args.i==1:
        Distance=20000
        try:
            os.mkdir('extract')
        except OSError:
            pass
        f1 = open(os.path.join('extract','Finished_Integrase_seq.fa'), 'ab')
    else:
        Distance=4000
        f1 = open(os.path.join('extract',str(args.i)+'-1-Integrase_seq.fa'), 'ab')
        f2 = open(os.path.join('extract',str(args.i) + '-2-Integrase_seq.fa'), 'ab')
    if args.target=='None':
        cmd1='ls *'+str(args.fasta)+' > Target.txt'
        os.system(cmd1)
        Target_file = 'Target.txt'
    else:
        Target_file=args.target
elif args.func_type == 1:
    if args.i==1:
        Distance=20000
        File_attC ='Finished_Integrase_seq.fa.Z.max.attc.hits.txt2.txt'
    else:
        Distance=4000
        File_attC1 = str(args.i)+'-1-Integrase_seq.fa.Z.max.attc.hits.txt2.txt'
        File_attC2 = str(args.i) + '-2-Integrase_seq.fa.Z.max.attc.hits.txt2.txt'
else:
    Attc_file_name=glob.glob(os.path.join('extract','*.new.txt2.txt'))
    try:
        os.mkdir(args.outputdir)
    except OSError:
        pass

################################################### Function ########################################################
def Long_integron_1(line, Sequence):
    if int(str(line).split('\t')[7])<=4000 or int(str(line).split('\t')[8])>=40000-4000:
        if str(line).split('\t')[0] not in Sequence:
            Sequence.append(str(line).split('\t')[0])


def Long_integron_2(line, Sequence):
    if str(line).split('\t')[0] not in Sequence:
        Sequence.append(str(line).split('\t')[0])


def Adjust_locus(line,Locus):
    if Locus>= 20000 + (int(args.i) - 1) * Distance: #locus of integrase
        return [int(str(line).split('\t')[7])-(20000+(int(args.i)-1)*Distance+1)+Locus,
                int(str(line).split('\t')[8]) - (20000 + (int(args.i) - 1) * Distance + 1) + Locus]
    else:
        return [int(str(line).split('\t')[7]),int(str(line).split('\t')[8])]


def filter_blast_list(file,Cutoff_identity,Cutoff_hitlength):
    try:
        f = open(file + '.filter.txt', 'w')
        for line in open(file, 'rb'):
            if AA_Length.get(str(line).split('\t')[0],'None')!='None':
                if float(str(line).split('\t')[2])>=Cutoff_identity:
                    if float(str(line).split('\t')[3])>=Cutoff_hitlength*float(AA_Length.get(str(line).split('\t')[0]))/100:
                        f.write(str(line))
        f.close()
        return file + '.filter.txt'
    except IOError:
        print str(file)+' missing!'

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
            cmd = 'cmsearch --tblout ' + str(files) + '.Z.max.attc.hits.txt -Z `du -m ' + str(files) + \
                  ' | cut -f1` --' + str(args.quick) + ' --cpu ' + str(args.thread) + ' -E ' + str(
                args.cutoff) + ' database/attc_4.cm.hmm ' + \
                  str(files)
            print cmd
            os.system(cmd)
            Result.append(str(files) + '.Z.max.attc.hits.txt')
        if file not in Newlist:
            Newfilelist=''
            for newfiles in Newlist:
                Newfilelist+=' '+newfiles
            Resultlist=''
            for files in Result:
                Resultlist+=' '+files
            os.system('cat '+str(Resultlist)+' > '+str(file) + '.Z.max.attc.hits.txt\n')
            os.system('rm -rf ' + str(Resultlist) + ' \n')
            os.system('rm -rf ' + str(Newfilelist) + ' \n')


################################################### Programme #######################################################
#for sequence-based blast result
ORFs_name = dict()
ORFs_anno = dict()
AA_Length = dict()
for line in open(os.path.join(in_dir, 'all.orfs.aa.length'), 'rb'):
    AA_Length.setdefault(str(line).split('\t')[0],float(str(line).split('\t')[2]))
    if args.AA_type == 0:
        ORFs_name.setdefault(str(line).split('\t')[0],
                             [float(str(line).split('\t')[1].split(' ')[-2]),
                              float(str(line).split('\t')[1].split(' ')[-1])])
    elif args.AA_type == 1:
        ORFs_name.setdefault(str(line).split('\t')[0],
                             [float(str(line).split('\t')[1].split(' # ')[1]),
                              float(str(line).split('\t')[1].split(' # ')[2])])
    ORFs_anno.setdefault(str(line).split('\t')[0],
                         str(line).split('\t')[1])
if args.i == 1:
    Int_file = filter_blast_list(Int_file, 80, 50)
else:
    Int_file=str(Int_file) + '.filter.txt'


#input integron-like sequences extracted from the head of input data
Warninglist=[]
try:
    for line in open(os.path.join(str(in_dir) + '/extract', args.head), 'r'):
        Warninglist.append(str(line).split('\n')[0])
except IOError:
    pass

#only chromosome integron sequences are extracted by locus comparison, also can be checked by the chr.fa name in the sequence.
#remove duplicated ones by locus comparison
#Locus2+Distance >Len(record.seq) is okay, only extract sequence to the end
if args.func_type==0:
    File_list = []
    Sequence_list=[]
    Extract_list=[]
    for line in open(Target_file,'r'):
        File_list.append(str(line).split('\n')[0].split(str(args.fasta))[0]+str(args.fasta))
        if args.i>1:
            Sequence_list.append(str(line).split('\n')[0].split(str(args.fasta)+'__')[1])
    if File_list!=[]:
        for line in open(Int_file,'rb'):
            ORF=str(line).split('\t')[0]
            for key in File_list:
                #the relation between file name and ORF name
                if key.split(args.fasta)[0] in ORF:
                    if str(key) + str(args.fasta) + '__' + str(ORF) not in Extract_list:
                        print str(key) + str(args.fasta) + '__' + str(ORF)
                        Extract_list.append(str(key) + str(args.fasta) + '__' + str(ORF))
                        if args.i>1:
                            if all(ORF not in Sequence for Sequence in Sequence_list):
                                break
                        Locus1 = int(ORFs_name.get(ORF, 'None')[0]) - 1
                        Locus2 = int(ORFs_name.get(ORF, 'None')[1]) - 1
                        if Locus1 > Locus2:
                            Locustemp = Locus1
                            Locus1 = Locus2
                            Locus2 = Locustemp
                        if args.i == 1:
                            if Locus1 >= Distance:
                                Locus1 = Locus1 - Distance
                            else:
                                Locus1 = 0
                                fwarning = open(
                                    os.path.join(str(in_dir) + '/extract', str(args.i) + '.warning.txt'),
                                    'ab')
                                fwarning.write(str(key + str(args.fasta)) + '__' + str(ORF) + '\n')
                                fwarning.close()
                            Locus2 = Locus2 + Distance
                            for record in SeqIO.parse(str(key), "fasta"):
                                if record.id in str(ORF):
                                    f1.write(">" + str(key) + str(args.fasta) + '__' + str(ORF) + \
                                             '\t' + str(ORFs_anno.get(ORF, 'None')) + '\n')
                                    f1.write(str(record.seq)[int(Locus1):int(Locus2)] + '\n')
                        else:
                            if Locus1 >= 20000 + (int(args.i) - 1) * Distance:
                                Locus1 = Locus1 - 20000 - (int(args.i) - 1) * Distance
                            else:
                                Locus1 = 0
                                fwarning = open(
                                    os.path.join(str(in_dir) + '/extract', str(args.i) + '.warning.txt'), 'ab')
                                fwarning.write(str(key + str(args.fasta)) + '__' + str(ORF) + '\n')
                                fwarning.close()
                            Locus2 = Locus2 + 20000 + (int(args.i) - 2) * Distance - 200
                            for record in SeqIO.parse(str(key), "fasta"):
                                if record.id in str(ORF):
                                    if str(key + str(args.fasta)) + '__' + str(ORF) not in Warninglist:
                                        f1.write(">" + str(key + str(args.fasta)) + '__' + str(ORF) + \
                                                 '\t' + str(ORFs_anno.get(ORF, 'None')) + '\n')
                                        f1.write(str(record.seq)[int(Locus1):int(Locus1 + Distance + 200)] + '\n')
                                    else:
                                        f1.write(">" + str(key + str(args.fasta)) + '__' + str(ORF) + \
                                                 '\t' + str(ORFs_anno.get(ORF, 'None')) + '\n')
                                        f1.write('\n')
                                    f2.write(">" + str(key + str(args.fasta)) + '__' + str(ORF) + \
                                             '\t' + str(ORFs_anno.get(ORF, 'None')) + '\n')
                                    f2.write(str(record.seq)[int(Locus2):int(Locus2 + Distance + 200)] + '\n')

        f1.close()
        if args.i!=1:
            f2.close()
            Cmlist=['extract/'+str(args.i)+'-1-Integrase_seq.fa','extract/'+str(args.i)+'-2-Integrase_seq.fa']
            Cmsearch(Cmlist)
        else:
            Cmsearch(['extract/Finished_Integrase_seq.fa'])
        cmd4='python scripts/attCformat.py -input extract\n'
        cmd5='python scripts/Integrase_extraction_seq.py -i '+str(args.i)+' --func_type 1'+ \
                 ' --fasta '+args.fasta+' --orf '+args.orf+ ' --AA_type '+str(args.AA_type)+' --outputdir '+\
             str(args.outputdir)+' --head '+str(args.i)+'.warning.txt\n'
        print cmd4
        print cmd5
        os.system(cmd4)
        os.system(cmd5)
    else:
        cmd7='python scripts/Integrase_extraction_seq.py --fasta '+str(args.fasta)+' --func_type 2'+ \
             ' --fasta '+args.fasta+' --orf '+args.orf+ ' --AA_type '+str(args.AA_type)+' --outputdir '\
             +str(args.outputdir)+' --head '+str(args.i)+'.warning.txt\n'
        print cmd7
        os.system(cmd7)
elif args.func_type == 1:
    Sequence_list=[]
    f3=open(os.path.join('extract', str(args.i)+'-Target.txt'), 'wb')
    if args.i == 1:
        f1 = open(os.path.join('extract', str(File_attC).replace('txt2.txt', 'new.txt2.txt')), 'wb')
        for line in open(os.path.join('extract', str(File_attC)), 'rb'):
            Long_integron_1(str(line), Sequence_list)
            for line1 in open(Int_file, 'rb'):
                ORF = str(line1).split('\t')[0]
                if ORF in str(line).split('\t')[0]:
                    Locus = int(ORFs_name.get(ORF, 'None')[0])
                    Locusnew = Adjust_locus(str(line), Locus)
                    f1.write(str(line).split('\t')[0] + '\t' + '\t'.join(
                        str(line).split('\t')[1:7]) + \
                             '\t' + str(Locusnew[0]) + '\t' + str(Locusnew[1]) + '\t' + \
                             str(line).split('\t')[9] + '\t' + '\t'.join(
                        str(line).split('\t')[10:]))
    else:
        f1 = open(os.path.join('extract', str(File_attC1).replace('txt2.txt', 'new.txt2.txt')), 'wb')
        f2 = open(os.path.join('extract', str(File_attC2).replace('txt2.txt', 'new.txt2.txt')), 'wb')
        for line in open(os.path.join('extract', str(File_attC1)), 'rb'):
            Long_integron_2(str(line), Sequence_list)
            for line1 in open(Int_file, 'rb'):
                ORF = str(line1).split('\t')[0]
                if ORF in str(line).split('\t')[0]:
                    Locus = int(ORFs_name.get(ORF, 'None')[0])
                    Locusnew = Adjust_locus(str(line), Locus)
                    if Locusnew[0] > 0:
                        f1.write(str(line).split('\t')[0] + '\t' + '\t'.join(
                            str(line).split('\t')[1:7]) + \
                                 '\t' + str(Locusnew[0]) + '\t' + str(Locusnew[1]) + '\t' + \
                                 str(line).split('\t')[9] + '\t' + '\t'.join(
                            str(line).split('\t')[10:]))
                    else:
                        try:
                            Sequence_list.remove(str(line).split('\t')[0])
                        except ValueError:
                            pass

        for line in open(os.path.join('extract', str(File_attC2)), 'rb'):
            Long_integron_2(str(line), Sequence_list)
            for line1 in open(Int_file, 'rb'):
                ORF = str(line1).split('\t')[0]
                if ORF in str(line).split('\t')[0]:
                    Locus = int(ORFs_name.get(ORF, 'None')[0])
                    Locusnew = Adjust_locus(str(line), Locus)
                    if Locusnew[0] > 0:
                        f2.write(str(line).split('\t')[0] + '\t' + '\t'.join(
                            str(line).split('\t')[1:7]) + \
                                 '\t' + str(Locusnew[0]) + '\t' + str(Locusnew[1]) + '\t' + \
                                 str(line).split('\t')[9] + '\t' + '\t'.join(
                            str(line).split('\t')[10:]))
                    else:
                        try:
                            Sequence_list.remove(str(line).split('\t')[0])
                        except ValueError:
                            pass
        f2.close()
    f1.close()
    for key in Sequence_list:
        if key not in Warninglist:
            f3.write(str(key)+'\n')
            Sequence_list.remove(key)
    f3.close()
    if Sequence_list!=[]:
        cmd6='python scripts/Integrase_extraction_seq.py -i '+str(int(args.i)+1)+' -input '+\
             str(args.input)+' -target extract/'+str(args.i)+'-Target.txt --fasta '+str(args.fasta)+' --func_type 0'+ \
             ' --fasta '+args.fasta+' --orf '+args.orf+ ' --AA_type '+str(args.AA_type)+' --outputdir '\
             +str(args.outputdir)+' --head '+str(args.i)+'.warning.txt\n'
        print cmd6
        os.system(cmd6)
    else:
        cmd7='python scripts/Integrase_extraction_seq.py --fasta '+str(args.fasta)+' --func_type 2'+ \
             ' --fasta '+args.fasta+' --orf '+args.orf+ ' --AA_type '+str(args.AA_type)+' --outputdir '\
             +str(args.outputdir)+' --head '+str(args.i)+'.warning.txt\n'
        print cmd7
        os.system(cmd7)
else:
    for Attc_file in Attc_file_name:
        for line in open(Attc_file,'r'):
            print os.path.join(args.outputdir,str(line).split(args.fasta)[0]+str(args.fasta)+'.Z.max.attc.hits.txt2.txt')
            f1=open(os.path.join(args.outputdir,str(line).split(args.fasta)[0]+str(args.fasta)+'.Z.max.attc.hits.txt2.txt'),'ab')
            f1.write(str(line))
            f1.close()







