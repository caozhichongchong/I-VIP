import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-input",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument("-target",
                    help="The taget sequence file", type=str, default='None',metavar='None')
parser.add_argument("--f",
                    help="file type or filename extension of your sequences",
                    type=str, default='.fa', metavar='.fa or .fasta or .fna')
parser.add_argument("--o",
                    help="Optional: to provide your own CDS or ORFs files, \
                                        please input the file type or filename extension of your CDS or ORFs",
                    type=str, default='.faa', metavar='.faa')
parser.add_argument('--ot',
                    help="Optional: to provide your own CDS or ORFs files, \
                                         please input the orf type or the method you used to extract orfs,\
                                         eg: 1 for genbank parsing; 2 for prodigal prediction",
                    metavar="a file containing target files and their specific ORF format",
                    action='store', default='Temp/all.orf.length', type=str)
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result',metavar='I-VIP_result')
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running I-VIP (default 15)",
                    metavar="15", action='store', default=15, type=int)
parser.add_argument('--q',
                    help="Optional: attC search setting \
                    (from least strict to most strict: \'max\', \'nohmm\', \'mid\', \'\', \'rfam\'),\
                     (default \'mid\')",
                    metavar="mid",
                    action='store', default='mid', type=str)
parser.add_argument('--c',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='Optional: set the e-value cutoff for attC search (default is 1.0)')
parser.add_argument('--cmsearch',
                    help="Optional: complete path to cmsearch if not in PATH,",
                    metavar="/usr/local/bin/cmsearch",
                    action='store', default='cmsearch', type=str)
#for ModuleA2 self-calling, please do not tune these parameters
parser.add_argument('-i',
                    help="The irratation time for ModuleA2.py",
                    metavar="1", action='store', default=1, type=int)
parser.add_argument("--func_type",
                    help="Function step of ModuleA2.py, \
                    eg: 0 for sequence extraction and extension; "
                         "1 for cmsearch result analysis; "
                         "2 for cmsearch result summarizing",
                    metavar="0 or 1 or 2",
                    choices=[0, 1, 2],
                    action='store', default=0, type=int)
parser.add_argument('--head',
                    help="The integron-like sequences that locate on the extremity of input sequence",
                    type=str, default='Warning.txt',metavar='Warning.txt')

################################################## Definition ########################################################
args = parser.parse_args()
input_path = args.input
outputdir = args.r + '/output'
# load the intI1 and sul1 search results
Int_file = (args.r + '/output/all.int.blast.txt.filter.txt')
Sul1_file = (args.r + '/output/all.sul1.blast.txt.filter.txt')
os.system('cat '+Int_file+' '+ Sul1_file+' > '+ args.r + '/output/all.int.sul1.blast.txt.filter.txt')
Int_file = (args.r + '/output/all.int.sul1.blast.txt.filter.txt')
# load input files for ModuleA2.py
if args.func_type == 0:
    # extract integron-like sequences
    if args.i==1:
        # first time extraction
        Distance=20000
        try:
            os.mkdir(str(input_path)+'/extract')
        except OSError:
            pass
        f1 = open(os.path.join(str(input_path)+'/extract','Finished_Integrase_seq.fa'), 'ab')
    else:
        # when the attC site was found at the extremity of an integron-like sequence, it will be extended at that extremity
        # extend an integron-like sequence for i times
        Distance=4000
        f1 = open(os.path.join(str(input_path)+'/extract',str(args.i)+'-1-Integrase_seq.fa'), 'ab')
        f2 = open(os.path.join(str(input_path)+'/extract',str(args.i) + '-2-Integrase_seq.fa'), 'ab')
    if args.target=='None':
        # set the target sequence files for first time extraction
        cmd1='ls '+str(input_path)+'/*'+str(args.f)+' > '+str(input_path)+'/extract/1-Target.txt'
        os.system(cmd1)
        Target_file = str(input_path)+'/extract/1-Target.txt'
    else:
        # set the target sequence files for sequence extension
        Target_file = args.target
elif args.func_type == 1:
    # analyze attC results
    if args.i==1:
        # set distance and attC result files of the last extension
        Distance=20000
        # File_attCstores the attC search results of integron-like sequences
        File_attC ='Finished_Integrase_seq.fa.Z.max.attc.hits.txt2.txt'
    else:
        # set distance and attC result files of the last extension
        Distance=4000
        # File_attC1, and File_attC2 store the attC search results of integron-like sequences (for each side)
        File_attC1 = str(args.i)+'-1-Integrase_seq.fa.Z.max.attc.hits.txt2.txt'
        File_attC2 = str(args.i) + '-2-Integrase_seq.fa.Z.max.attc.hits.txt2.txt'
else:
    # Module2.py will be finished
    # summarizing and formatting all the attC results
    # load all attC search results
    Attc_file_name=glob.glob(os.path.join(str(input_path)+'/extract','*.new.txt2.txt'))
    try:
        os.mkdir(outputdir)
    except OSError:
        pass
# record Module A2 process
flog = open('ModuleA2.log', 'wb')

################################################### Function ########################################################


def Long_integron_1(line, Sequence):
    # decide whether to extend an integron-like sequence during first time extraction
    if int(str(line).split('\t')[7])<=4000 or int(str(line).split('\t')[8])>=40000-4000:
        if str(line).split('\t')[0] not in Sequence:
            Sequence.append(str(line).split('\t')[0])


def Long_integron_2(line, Sequence):
    # decide whether to extend an integron-like sequence after first time extraction
    if str(line).split('\t')[0] not in Sequence:
        Sequence.append(str(line).split('\t')[0])


def Adjust_locus(line,Locus):
    # correct the loci of integron-like sequences in the original target sequences
    if Locus>= 20000 + (int(args.i) - 1) * Distance: #locus of integrase
        return [int(str(line).split('\t')[7])-(20000+(int(args.i)-1)*Distance+1)+Locus,
                int(str(line).split('\t')[8]) - (20000 + (int(args.i) - 1) * Distance + 1) + Locus]
    else:
        return [int(str(line).split('\t')[7]),int(str(line).split('\t')[8])]


def Cmsearch(list_file):
    #search for attC
    for file in list_file:
        Newlist = []
        #calculate the total data size (bp) for query data
        Total_length = 0
        for record in SeqIO.parse(file, 'fasta'):
            Total_length += len(record.seq)
        # large data is firstly split into 2Mb files for attC search
        # for parallel comparison of cmserach results (e-value)
        if Total_length >= 2000000:
            Length = 0
            for record in SeqIO.parse(file, 'fasta'):
                Length += len(record.seq)
                f1 = open(str(file) + '_' + str(int(Length / 2000000)), 'ab')
                f1.write('>' + str(record.id) + '\t' + str(record.description) + '\n')
                f1.write(str(record.seq) + '\n')
                f1.close()
                if str(file) + '_' + str(int(Length / 2000000)) not in Newlist:
                    Newlist.append(str(file) + '_' + str(int(Length / 2000000)))
        else:
            Newlist.append(file)
        #attC search by
        Result = []
        for files in Newlist:
            cmd = args.cmsearch + ' --tblout ' + str(files) + '.Z.max.attc.hits.txt' + ' -Z `du -m ' \
                  + str(files)+ ' | cut -f1` --' + str(args.q) + ' --cpu ' + str(args.t) + ' -E ' \
                  + str(args.c) + ' database/attc_4.cm.hmm ' + str(files)
            os.system(cmd)
            Result.append(str(files) + '.Z.max.attc.hits.txt')
        # merge search results of split files
        if file not in Newlist:
            Newfilelist = ''
            for newfiles in Newlist:
                Newfilelist += ' ' + newfiles
            Resultlist = ''
            for files in Result:
                Resultlist += ' ' + files
            os.system(
                'cat ' + str(Resultlist) + ' > ' + str(file) + '.Z.max.attc.hits.txt \n')
            os.system('rm -rf ' + str(Resultlist) + ' \n')
            os.system('rm -rf ' + str(Newfilelist) + ' \n')


def ORFinfo(file, ORFs_anno):
    #set up ORF information
    for line in open(file, 'rb'):
        ORFs_anno.setdefault(str(line).split('\t')[0],
                             [str(line).split('\t')[1], float(str(line).split('\t')[-3]),
                                  float(str(line).split('\t')[-2])])


################################################### Programme #######################################################
# input integron-like sequences extracted from the head of input data
Warninglist=[]
try:
    for line in open(os.path.join(str(input_path) + '/extract', args.head), 'r'):
        Warninglist.append(str(line).split('\n')[0])
except IOError:
    pass

# start main programme
if args.func_type==0:
    # load ORF format of files
    ORFs_anno = dict()
    ORFinfo(args.ot, ORFs_anno)
    # extract integron-like sequences
    File_list = [] # all target sequences for first extraction
    Sequence_list = [] # all target sequences for extension
    Extract_list = [] # sequences to be extracted
    # set up target sequences
    for line in open(Target_file, 'rb'):
        in_dir, input_file = os.path.split(line.split('\n')[0])
        File_list.append(str(input_file).split(str(args.f))[0])
        if args.i >1 :
            Sequence_list.append(str(line).split('\n')[0])
    # locate the sequences for extraction
    for line in open(Int_file, 'rb'):
        ORF=str(line).split('\t')[0]
        # an integrase or sul1 ORF
        for key in File_list:
            if key in ORF:
                # identify sequences containing integrases and sul1 genes
                Sequenceid = str(key) + str(args.f) + '____' + str(ORF)
                if Sequenceid not in Extract_list:
                    # input the sequence in to-be-extracted list
                    # name the sequence as filename + filename extension + '____' + ORF name
                    Extract_list.append(Sequenceid)
                    if args.i > 1:
                        # for extension, to avoid repeatedly extracting ORFs-containing sequences
                        if not all(ORF not in Sequence for Sequence in Sequence_list):
                            break
                    # set loci of target sequences for extraction
                    Locus1 = int(ORFs_anno.get(ORF, 'None')[1]) - 1
                    Locus2 = int(ORFs_anno.get(ORF, 'None')[2]) - 1
                    if Locus1 > Locus2:
                        Locustemp = Locus1
                        Locus1 = Locus2
                        Locus2 = Locustemp
                    if args.i == 1:
                        # first extraction
                        # set extracting locus
                        if Locus1 >= Distance:
                            # if Locus1 has not reached the extremity
                            Locus1 = Locus1 - Distance
                        else:
                            # if Locus1 has reached the extremity
                            Locus1 = 0
                            # to record that for this sequence, one extremity has been reached
                            fwarning = open(os.path.join(str(input_path) + '/extract', str(args.i) + '.warning.txt'),
                                            'ab')
                            fwarning.write(Sequenceid + '\n')
                            fwarning.close()
                        Locus2 = Locus2 + Distance + 1
                        # accept Locus2 + Distance >Len(record.seq), but only extract sequence to the end
                        # extract the integron-like sequences
                        for record in SeqIO.parse(os.path.join(str(input_path), key + str(args.f)), "fasta"):
                            if record.id in Sequenceid:
                                f1.write((">" + Sequenceid + \
                                         '\t' + str(ORFs_anno.get(ORF, 'None')[0])).replace('\r','').replace('\n','') + '\n')
                                f1.write(str(record.seq)[int(Locus1):int(Locus2)].replace('\r','').replace('\n','') + '\n')
                    else:
                        # sequence extension
                        # set loci of target sequences for extension
                        if Locus1 >= 20000 + (int(args.i) - 1) * Distance:
                            # if Locus1 has not reached the extremity
                            Locus1 = Locus1 - 20000 - (int(args.i) - 1) * Distance
                        else:
                            # if Locus1 has reached the extremity
                            Locus1 = 0
                            # to record that for this sequence, one extremity has been reached
                            fwarning = open(os.path.join(str(input_path) + '/extract', str(args.i) + '.warning.txt'),
                                            'ab')
                            fwarning.write(Sequenceid + '\n')
                            fwarning.close()
                        Locus2 = Locus2 + 20000 + (int(args.i) - 2) * Distance - 200
                        # extend the integron-like sequences
                        for record in SeqIO.parse(os.path.join(str(input_path), key + str(args.f)), "fasta"):
                            if record.id in Sequenceid:
                                if Sequenceid not in Warninglist:
                                    # not reaching the  extremity, output sequence
                                    f1.write((">" + Sequenceid + \
                                             '\t' + str(ORFs_anno.get(ORF, 'None')[0])).replace('\r','').replace('\n','') + '\n')
                                    f1.write(str(record.seq)[int(Locus1):int(Locus1 + Distance + 201)].replace('\r','').replace('\n','') + '\n')
                                else:
                                    # reaching the  extremity, forbid output sequence
                                    f1.write((">" + Sequenceid + \
                                             '\t' + str(ORFs_anno.get(ORF, 'None')[0])).replace('\r','').replace('\n','') + '\n')
                                    f1.write('A\n')
                                f2.write((">" + Sequenceid + \
                                         '\t' + str(ORFs_anno.get(ORF, 'None')[0])).replace('\r','').replace('\n','') + '\n')
                                if str(record.seq)[int(Locus2):int(Locus2 + Distance + 201)]=='':
                                    # reaching the  extremity, no output sequence
                                    f2.write('A\n')
                                else:
                                    f2.write(str(record.seq)[int(Locus2):int(Locus2 + Distance + 201)].replace('\r','').replace('\n','') + '\n')
    f1.close()
    # search attC in integron-like sequences (Cmlist containing the extracted fasta files)
    if args.i!=1:
        f2.close()
        Cmlist=[str(input_path)+'/extract/'+str(args.i)+'-1-Integrase_seq.fa',str(input_path)+'/extract/'+str(args.i)+'-2-Integrase_seq.fa']
        Cmsearch(Cmlist)
    else:
        Cmsearch([str(input_path)+'/extract/Finished_Integrase_seq.fa'])
    # format attC search results
    cmd4 = 'python scripts/ModuleA1.py -input ' + str(input_path) + '/extract '+ ' -f ' + str(args.f) + \
           ' --r ' + str(args.r) +' --m 2 --fn 2 \n'
    # call ModuleA2.py for attC results analysis
    cmd5 = 'python scripts/ModuleA2.py -input ' + str(input_path) + ' -i ' + str(
        args.i) + ' --func_type 1' + ' --r ' + str(args.r) + \
           ' --f ' + args.f + ' --o ' + args.o + ' --ot ' + str(args.ot) + \
             ' --t ' + str(args.t) + ' --c ' + str(args.c) + \
           ' --q ' + str(args.q) + ' --cmsearch ' + str(args.cmsearch) + \
           ' --head ' + str(args.i) + '.warning.txt \n'
    print cmd4
    print cmd5
    del ORFs_anno
    ORFs_anno = dict()
    os.system(cmd4)
    os.system(cmd5)
    flog.write(cmd4 + '\n' + cmd5 + '\n')
elif args.func_type == 1:
    # load ORF format of files
    ORFs_anno = dict()
    ORFinfo(args.ot, ORFs_anno)
    # analyze attC results to decide whether to further extend the sequences
    # Sequence_list to store all target sequences for extension and output to f3
    Sequence_list = []
    f3 = open(os.path.join(str(input_path) + '/extract', str(args.i+1) + '-Target.txt'), 'ab')
    Target_file=str(input_path) + '/extract/'+str(args.i)+'-Target.txt'
    for line in open(Target_file, 'r'):
        in_dir, input_file = os.path.split(line.split('\n')[0])
    # start analyze attC reasults
    if args.i == 1:
        # after first extraction
        # File_attC stores the attC search results of integron-like sequences
        f1 = open(os.path.join(str(input_path) + '/extract', str(File_attC).replace('txt2.txt', 'new.txt2.txt')), 'ab')
        # correct the loci of attC search results into the loci on the original sequences
        for line in open(os.path.join(str(input_path) + '/extract', str(File_attC)), 'rb'):
            Long_integron_1(str(line), Sequence_list)
            # decide whether an integron-like sequence needs to be extended
            for line1 in open(Int_file, 'rb'):
                ORF = str(line1).split('\t')[0]
                if ORF in str(line).split('\t')[0]:
                    Locus = int(ORFs_anno.get(ORF, 'None')[1])
                    Locusnew = Adjust_locus(str(line), Locus)
                    f1.write(str(line).split('\t')[0] + '\t' + '\t'.join(
                        str(line).split('\t')[1:7]) + \
                             '\t' + str(Locusnew[0]) + '\t' + str(Locusnew[1]) + '\t' + \
                             str(line).split('\t')[9] + '\t' + '\t'.join(
                        str(line).split('\t')[10:]))
    else:
        # after extension
        # File_attC1, and File_attC2 store the attC search results of integron-like sequences (for each side)
        f1 = open(os.path.join(str(input_path) + '/extract', str(File_attC1).replace('txt2.txt', 'new.txt2.txt')), 'ab')
        f2 = open(os.path.join(str(input_path) + '/extract', str(File_attC2).replace('txt2.txt', 'new.txt2.txt')), 'ab')
        # File_attC1
        for line in open(os.path.join(str(input_path) + '/extract', str(File_attC1)), 'rb'):
            Long_integron_1(str(line), Sequence_list)
            # decide whether an integron-like sequence needs to be extended
            # correct the loci of attC search results into the loci on the original sequences
            for line1 in open(Int_file, 'rb'):
                ORF = str(line1).split('\t')[0]
                if ORF in str(line).split('\t')[0]:
                    Locus = int(ORFs_anno.get(ORF, 'None')[1])
                    Locusnew = Adjust_locus(str(line), Locus)
                    if Locusnew[0] > 0:
                        # correct loci
                        f1.write(str(line).split('\t')[0] + '\t' + '\t'.join(
                            str(line).split('\t')[1:7]) + \
                                 '\t' + str(Locusnew[0]) + '\t' + str(Locusnew[1]) + '\t' + \
                                 str(line).split('\t')[9] + '\t' + '\t'.join(
                            str(line).split('\t')[10:]))
                    else:
                        # integron-like sequences should not be extended
                        try:
                            Sequence_list.remove(str(line).split('\t')[0])
                        except ValueError:
                            pass
        # File_attC2
        for line in open(os.path.join(str(input_path) + '/extract', str(File_attC2)), 'rb'):
            Long_integron_2(str(line), Sequence_list)
            # decide whether an integron-like sequence needs to be extended
            # correct the loci of attC search results into the loci on the original sequences
            for line1 in open(Int_file, 'rb'):
                ORF = str(line1).split('\t')[0]
                if ORF in str(line).split('\t')[0]:
                    Locus = int(ORFs_anno.get(ORF, 'None')[1])
                    Locusnew = Adjust_locus(str(line), Locus)
                    if Locusnew[0] > 0:
                        # correct loci
                        f2.write(str(line).split('\t')[0] + '\t' + '\t'.join(
                            str(line).split('\t')[1:7]) + \
                                 '\t' + str(Locusnew[0]) + '\t' + str(Locusnew[1]) + '\t' + \
                                 str(line).split('\t')[9] + '\t' + '\t'.join(
                            str(line).split('\t')[10:]))
                    else:
                        # integron-like sequences should not be extended
                        try:
                            Sequence_list.remove(str(line).split('\t')[0])
                        except ValueError:
                            pass
        f2.close()
    f1.close()
    # output sequences for the next extension
    Output = 0
    for key in Sequence_list:
        if key not in Warninglist:
            f3.write(str(key) + '\n')
            Output = 1
    f3.close()
    if Sequence_list != [] and Output == 1:
        # if next extension is necessary, call ModuleA2
        cmd6 = 'python scripts/ModuleA2.py -input ' + str(input_path) + ' -i ' + str(
            int(args.i) + 1) + ' --r ' + str(args.r) + \
               ' -target ' + str(input_path) + '/extract/' + str(args.i+1) + '-Target.txt --f ' + str(
            args.f) + ' --func_type 0' + \
            ' --o ' + args.o + ' --ot ' + str(args.ot) + \
               ' --t ' + str(args.t) + ' --c ' + str(args.c) + \
               ' --q ' + str(args.q) + ' --cmsearch ' + str(args.cmsearch) + \
               ' --head ' + str(args.i) + '.warning.txt \n'
        print cmd6
        flog.write(cmd6 + '\n')
        del ORFs_anno
        ORFs_anno = dict()
        os.system(cmd6)
    else:
        # if all extension is finished, call ModuleA2 for sammarizing
        cmd7 = 'python scripts/ModuleA2.py -input ' + str(input_path) + ' --f ' + str(
            args.f) + ' --func_type 2' + ' --r ' + str(args.r) +\
             ' --o ' + args.o + ' --ot ' + str(args.ot) + \
               ' --t ' + str(args.t) + ' --c ' + str(args.c) + \
               ' --q ' + str(args.q) + ' --cmsearch ' + str(args.cmsearch) + \
               ' --head ' + str(args.i) + '.warning.txt \n'
        print cmd7
        flog.write(cmd7 + '\n')
        del ORFs_anno
        ORFs_anno = dict()
        os.system(cmd7)
else:
    # Module2.py will be finished
    # summarizing and formatting all the attC results
    for Attc_file in Attc_file_name:
        for line in open(Attc_file, 'r'):
            f1 = open(os.path.join(outputdir,
                                   str(line).split(args.f)[0] + str(args.f) + '.Z.max.attc.hits.txt2.txt'),
                      'ab')
            f1.write(str(line))
            f1.close()
    # os.system('-rm -rf ' +str(input_path)+'/extract \n')
    print 'ModuleA2.py finished!'
    flog.write('Finished preparation!!\n')
    flog.close()






