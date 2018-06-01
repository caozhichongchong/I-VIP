import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input',metavar='input')
parser.add_argument("-f",
                    help="file type or filename extension of your sequences",
                    type=str, default='.fa',metavar='.fa or .fasta or .fna')
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result',metavar='I-VIP_result')
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
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running I-VIP (default 15)",
                    metavar="15", action='store', default=15, type=int)
parser.add_argument('--u',
                    help="Optional: use two-step method for integrase and sul1 search," +
                         " \'None\' for using one step, \'usearch\' or \'diamond\' for using two-step \
                         (complete path to usearch or diamond if not in PATH, \
                         please make sure the search tools can be directly called), (default: \'None\')",
                    metavar="None or usearch",
                    action='store', default='usearchv8', type=str)
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
input_path = os.path.abspath(args.i)
list_data = glob.glob(os.path.join(input_path,'*'+args.f))
search_path = args.r + '/output'
# record preparation process
flog = open('Preparation.log', 'wb')
# load integrases and sul1 reference length
Length = dict()
for line in open(os.path.join('database/' + 'Integrase.fasta.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
for line in open(os.path.join('database/' + 'sul1_database.txt.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
################################################### Function ########################################################



def addname(filedir, file_name):
    Fasta_name = open(os.path.join(filedir,file_name), 'r')
    f = open(args.r + '/Temp/'+file_name + '.add', 'w')
    in_dir, input_file = os.path.split(file_name)
    for record in SeqIO.parse(Fasta_name, 'fasta'):
        if len(str(record.seq).replace(' ',''))>0:
            # remove empty ORF sequences, otherwise there could be a problem for blast, usearch and diamond
            f.write('>'+str(input_file)+'_'+str(record.id) + '\t' + str(record.description).replace('\t', ' ') + '\n' + str(
                str(record.seq)) + '\n')
    f.close()
    return args.r + '/Temp/'+file_name + '.add'


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
        flog.write(str(file) + ' missing!' + '\n')


def Extractaa(root, searchfile, orffile, gene):
    # extract the query aa sequences according to a usearch or diamond alignment output
    # generate a smaller data of potential intI1 or sul1 for blastp search
    # input the query ORF sequences
    AA_seq = dict()
    try:
        for record in SeqIO.parse(open(os.path.join(root, orffile), 'r'), 'fasta'):
            AA_seq.setdefault(str(record.id), str(record.seq))
        f1 = open(os.path.join(root, orffile + '.' + str(gene) + '.aa'), 'wb')
        try:
            for line in open(os.path.join(root, searchfile), 'rb'):
                try:
                    AA = str(line).split('\t')[0].split(' ')[0]
                    if AA_seq[AA] != '':
                        # avoid duplicate ORF
                        f1.write('>' + AA + '\n' +
                                 str(AA_seq[AA]) + '\n')
                        AA_seq[AA]=''
                except KeyError:
                    print 'AA not found for ' + AA
                    flog.write('AA not found for ' + AA + '\n')
        except IOError:
            pass
        f1.close()
        return os.path.join(root, orffile + '.' + str(gene) + '.aa')
    except (IOError):
        print 'Files were missing: ' + orffile
        flog.write('Files were missing: ' + orffile + '\n')


def searchintsul(filename):
    # search integrase and sul1 for each file
    if args.u != 'None':
        # two-step search
        if 'usearch' in args.u:
            # Start search integrase and sul1 genes by usearch
            cmdint1 = str(args.u) + " -ublast " + args.r + "/Temp/"+filename+ \
                    " -db database/Integrase.fasta -evalue 1e-0 -accel 0.5 -blast6out " \
                    + args.r + "/Temp/"+filename+".int.usearch.txt -threads " + str(args.t) + "\n"
            cmdsulI1 = str(args.u) + " -ublast " + args.r + "/Temp/"+filename+ \
                    " -db database/sul1_database.txt -evalue 1e-0 -accel 0.5 -blast6out " \
                    + args.r + "/Temp/"+filename+".sul1.usearch.txt -threads " + str(args.t) + "\n"
        elif "diamond" in args.u:
            # Start search integrase and sul1 genes by diamond!
            cmdint1 = str(args.u) +" blastp --query " + args.r + "/Temp/"+filename + \
                      " --db database/Integrase.fasta --out " + args.r + "/Temp/"+filename+\
                      ".int.usearch.txt --outfmt 6 --max-target-seqs 1 --evalue 1e-2 --threads " + str(args.t) + " \n"
            cmdsulI1 = str(args.u) +" blastp --query " + args.r + "/Temp/"+filename + \
                      " --db database/sul1_database.txt --out " + args.r + "/Temp/"+filename+\
                      ".sul1.usearch.txt --outfmt 6 --max-target-seqs 1 --evalue 1e-2 --threads " + str(args.t) + " \n"
        os.system(cmdint1)
        os.system(cmdsulI1)
        Int_file = Extractaa( args.r + '/Temp', filename+".int.usearch.txt", filename,'int')
        SulI_file = Extractaa(args.r + '/Temp', filename + ".sul1.usearch.txt", filename, 'sul1')
    else:
        # one-step search
        SulI_file = Int_file = args.r + '/Temp/' + filename
    # blast search
    cmdint2 = str(args.blastp) +" -query " + str(Int_file) + " -db database/Integrase.fasta -out " + str(search_path) \
            + "/"+filename+".int.blast.txt  -outfmt 6 -max_target_seqs 1 -evalue 1e-2 -num_threads " + \
            str(args.t) + " \n"
    os.system(cmdint2)
    cmdsulI2 = str(args.blastp) +" -query " + str(SulI_file) + " -db database/sul1_database.txt -out " + str(search_path) \
            + "/"+filename+".sul1.blast.txt  -outfmt 6 -max_target_seqs 1 -evalue 1e-2 -num_threads " + \
            str(args.t) + " \n"
    os.system(cmdsulI2)


def Calculate_length(file_name, ORF_format):
    # calculate the sequence length of ORFs
    if ORF_format == 'None':
        print 'No ORF format information for ' + file_name + '!'
        flog.write('No ORF format information for ' + file_name + '!' + '\n')
    else:
        f1=open(args.r + '/Temp/all.orf.length','ab')
        try:
            Fasta_name = open(file_name, 'r')
            for record in SeqIO.parse(Fasta_name, 'fasta'):
                if ORF_format == 1:
                    l1 = float(str(record.description).split('\t')[1].split(' ')[-2])
                    l2 = float(str(record.description).split('\t')[1].split(' ')[-1])
                elif ORF_format == 2:
                    l1 = float(str(record.description).split(' # ')[1])
                    l2 = float(str(record.description).split(' # ')[2])
                # output in the format of ORFname, annotation, loci1, loci2, length
                f1.write(str(record.id) + '\t \t' +
                        str(l1) + '\t' + str(l2) + '\t'+str(
                    len(record.seq)) + '\n')
        except (IOError):
            print 'Files were missing: ' + file_name
            flog.write('Files were missing: ' + file_name + '\n')
        f1.close()



################################################### Programme #######################################################
# check the format of ORF, genbank parsing or prodigal prediction
fot=open(args.r + '/Temp/ORF_format.log','wb')
for file_name in list_data:
    # check orf file
    filedir, file_name = os.path.split(file_name)
    try:
        orf_name=file_name.replace(args.f,args.o)
        f1 = open(os.path.join(filedir,orf_name),'rb')
        for line in f1:
            if args.ot == 1:
                try:
                    l1 = float(str(line).split(' ')[-2])
                    l2 = float(str(line).split(' ')[-1])
                    # output orf_name contains the dir information
                    orf_name = addname(filedir, orf_name)
                    Calculate_length(orf_name, 1)
                    fot.write(file_name + '\t1\n')
                except ValueError:
                    try:
                        l1 = float(str(line).split(' # ')[1])
                        l2 = float(str(line).split(' # ')[2])
                        # output orf_name contains the dir information
                        orf_name = addname(filedir, orf_name)
                        Calculate_length(orf_name, 2)
                        fot.write(file_name + '\t2\n')
                    except ValueError:
                        print "There\'s some problem with the ORF files!\n" \
                              "Please check whether they are in the form of genome parsing (--ot 1)" \
                              " or prodigal prediction (--ot 2)\n"
            elif args.ot == 2:
                try:
                    l1 = float(str(line).split(' # ')[1])
                    l2 = float(str(line).split(' # ')[2])
                    # output orf_name contains the dir information
                    orf_name = addname(filedir, orf_name)
                    Calculate_length(orf_name, 2)
                    fot.write(file_name + '\t2\n')
                except ValueError:
                    print "There\'s some problem with the ORF files!\n" \
                          "Please check whether they are in the form of genome parsing (--ot 1)" \
                          " or prodigal prediction (--ot 2)\n"
            break
        f1.close()
    except IOError or ValueError:
        # if no orf file or not valid orf file
        # prodigal predict orf
        os.system(args.prodigal+" -q  -a "+os.path.join(filedir,orf_name)
                  +" -i "+os.path.join(filedir,file_name)+"\n")
        # output orf_name contains the dir information
        orf_name = addname(filedir, orf_name)
        Calculate_length(orf_name, 2)
        fot.write(file_name + '\t2\n')
    orf_dir, orf_name=os.path.split(orf_name)
    searchintsul(orf_name)
fot.close()
# organizing integrase and sul1 blast results
os.system('cat '+search_path+'/*.sul1.blast.txt > '+search_path+'/all.sul1.blast.txt \n')
os.system('cat '+search_path+'/*.int.blast.txt > '+search_path+'/all.int.blast.txt \n')
# filter blast results
filter_blast_list(str(search_path) + '/all.sul1.blast.txt', 50, 50)
filter_blast_list(str(search_path) + '/all.int.blast.txt', 80, 50)
# separate filtered blast results into individual files, corresponding to the target sequences
fsul1 = open(os.path.join(search_path, 'all.sul1.blast.txt.filter.txt'))
for line in fsul1:
    Filename = str(line).split('\t')[0].split(args.o)[0].split(str(args.f))[0]
    fnew = open(os.path.join(search_path, Filename + str(args.f) + '.sul1.txt'), 'ab')
    fnew.write(line)
    fnew.close()
fint1 = open(os.path.join(search_path, 'all.int.blast.txt.filter.txt'))
for line in fint1:
    Filename = str(line).split('\t')[0].split(args.o)[0].split(str(args.f))[0]
    fnew = open(os.path.join(search_path, Filename + str(args.f) + '.int.txt'), 'ab')
    fnew.write(line)
    fnew.close()
print 'Finished preparation!'
flog.write('Finished preparation!!\n')
flog.close()