import os
import glob
import argparse
from Bio import SeqIO
import copy


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument('--g',
                    default="None", action='store', type=str, metavar='gbff or None',
                    help="Optional: to provide your own genbank files for gene cassettes annotation,\
                                         please set the directory of Genbank files \
                                        (default: \'None\' for No gbff input)")
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result',metavar='I-VIP_result')
parser.add_argument("--f",
                    help="file type or filename extension of your sequences",
                    type=str, default='.fa', metavar='.fa or .fasta or .fna')
parser.add_argument("--o",
                    help="Optional: to provide your own CDS or ORFs files, \
                                        please input the file type or filename extension of your CDS or ORFs",
                    type=str, default='.faa', metavar='.faa')


################################################## Definition ########################################################
args = parser.parse_args()
file_path = os.path.abspath(args.i)
AA_path = args.r + '/Temp'
gbff_path = args.g
input_path = os.path.abspath(args.r + '/result')
input_format = '.Integron.txt'
list_fasta = glob.glob(os.path.join(input_path,'*'+input_format))
resultdir = args.r + "/Integron"
# load ORF format of files
ORFformat = dict()
for lines in open(args.r + '/Temp/ORF_format.log', 'r'):
    ORFformat.setdefault(lines.split('\t')[0], int(lines.split('\t')[-1].split('\n')[0]))
# generate dirs for results, including ORFs, Integron_seqs (integron),
# Integrase_seqs (integrase sequences), Genbank_annotation (genbank annotation), Integron_structure (integron structure)
try:
    os.mkdir(resultdir)
except OSError:
    pass
try:
    # dir for class 1 integrons (Types A and B)
    os.mkdir(str(resultdir) + '/ClassI')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/ClassI/ORFs')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/ClassI/Integron_seqs')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/ClassI/Integrase_seqs')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/ClassI/Genbank_annotation')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/ClassI/Integron_structure')
except OSError:
    pass
try:
    # dir for integrons of other classes (Types C, D, E)
    os.mkdir(str(resultdir) + '/Other')
except OSError:
    pass
try:
    os.mkdir(str(resultdir)+'/Other/ORFs')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/Other/Integron_seqs')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/Other/Integrase_seqs')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/Other/Genbank_annotation')
except OSError:
    pass
try:
    os.mkdir(str(resultdir) + '/Other/Integron_structure')
except OSError:
    pass


################################################### Integron class #########################################################
__metaclass__=type
class Integron:
    # create a class to store integrons
    'a class to store integrons, with the information of \
    name, type, int, sulI, l1 (locus1), l2 (locus2), structure'
    def setName(self,name):
        self.name=name

    def setType(self,type):
        self.type=type

    def setInt(self):
        self.int=[]

    def addInt(self,Int):
        self.int.append(float(Int))

    def setL1(self,Locus1):
        # locus1 of integron
        self.l1=float(Locus1)
        self.newl1 = self.l1

    def changeL1(self,Locus1):
        # additional locus of potential intergrase
        if self.l1>float(Locus1):
            self.newl1=float(Locus1)
            self.l1 = self.newl1

    def setL2(self,Locus2):
        # locus2 of integron
        self.l2=float(Locus2)
        self.newl2 = self.l2

    def changeL2(self,Locus2):
        # additional locus of potential intergrase
        if self.l2<float(Locus2):
            self.newl2 = float(Locus2)
            self.l2 = self.newl2

    def setIntlocus(self):
        self.intlocus=[]

    def addIntlocus(self,Int):
        # no duplicated locus
        if float(Int) not in self.intlocus:
            self.intlocus.append(float(Int))

    def setattClocus(self):
        self.attclocus=[]

    def addattClocus(self,attc):
        # no duplicated locus
        if float(attc) not in self.attclocus:
            self.attclocus.append(float(attc))

    def setIntanno(self):
        self.intanno1 = []
        self.intanno2=dict()

    def addIntanno(self,locus,annotation):
        self.intanno1.append(locus)
        if self.intanno2.get(locus, 'None') == 'None':
            self.intanno2.setdefault(locus, [str(annotation).replace('[', '').replace(']', '')])
        else:
            self.intanno2[locus].append(str(annotation).replace('[', '').replace(']', ''))

    def setInttag(self):
        self.inttag1 = []
        self.inttag2=dict()

    def addInttag(self,locus,tag):
        self.inttag1.append(float(locus))
        if self.inttag2.get(float(locus), 'None') == 'None':
            self.inttag2.setdefault(float(locus), [str(tag).replace('[', '').replace(']', '')])
        else:
            self.inttag2[float(locus)].append(str(tag).replace('[', '').replace(']', ''))

    def setSul1(self):
        self.sul1=[]

    def addSul1(self,sul1):
        self.sul1.append(float(sul1))

    def setSul1locus(self):
        self.sul1locus=[]

    def addSul1locus(self,Sul1l):
        # no duplicated locus
        if float(Sul1l) not in self.sul1locus:
            self.sul1locus.append(float(Sul1l))

    def setStructure(self,length):
        # set up integron structure for integron classification
        # decide whether to include one extra ORF on the extremity
        Before=-1
        After=1
        if self.type == 'A':
            if max(self.int)==(length-1):
                # integrase as the end
                After=0
            if min(self.int)==0:
                # integrase as the start
                Before=0
            if max(self.sul1)==(length-1):
                # sul1 as the end
                After=0
            if min(self.sul1)==0:
                # sul1 as the start
                Before=0
        elif self.type=='C' or self.type=='B':
            if max(self.int)==(length-1):
                # integrase as the end
                After=0
            if min(self.int)==0:
                # integrase as the start
                Before=0
        elif self.type=='D':
            if max(self.sul1)==(length-1):
                # sul1 as the end
                After=0
            if min(self.sul1)==0:
                # sul1 as the start
                Before=0
        self.structure = [Before, After]

    def setOutput(self):
        # this integron has not been output
        self.output = 'FALSE'

    def doOutput(self):
        # this integron has been output
        self.output = 'TRUE'

    def deleteKey(self):
        # delete Type E integrons with only one attC
        if self.type=='E' and len(self.attclocus)==1:
            self.delete='TRUE'
        else:
            self.delete='FALSE'


################################################### Function #########################################################


def int_integron(file, AA_format):
    # main function of integron extraction
    for line in open(os.path.join(in_dir, file), 'r'):
        if ':' in str(line).split('\t')[1]:
            # set up all integrons in the format of a new class
            Keyintegron = Integron()
            Keyintegron.setName(str(line).split('\t')[1])
            Keyintegron.setType(str(line).split('\t')[0])
            Keyintegron.setattClocus()
            Keyintegron.setInt()
            Keyintegron.setIntlocus()
            Keyintegron.setIntanno()
            Keyintegron.setInttag()
            Keyintegron.setSul1()
            Keyintegron.setSul1locus()
            # set up locus for an integron, label the actucal position in the site list
            # position start from 0, thus max position = total position -1
            # locusmin and locusmax for integron range, locusele for the loci of integron elements
            locusmin = locusmax = locusele = 0
            for key in str(line).split('\t')[4].split(';'):
                if 'attC' in key or 'SulI' in key or 'Integrase' in key or 'IntI1' in key:
                    # set integron max locus position
                    locusmax = locusele
                    if 'attC' in key:
                        # set attC loci
                        Keyintegron.addattClocus(str(line).split('\t')[3].split(';')[locusele])
                    if 'Integrase' in key:
                        # set integrase loci
                        Keyintegron.addIntlocus(str(line).split('\t')[3].split(';')[locusele])
                        Keyintegron.addInt(locusele)
                        Keyintegron.addInttag(str(line).split('\t')[3].split(';')[locusele],
                                              'Integrase')
                    if 'IntI1' in key:
                        # set intI1 loci
                        Keyintegron.addIntlocus(str(line).split('\t')[3].split(';')[locusele])
                        Keyintegron.addInt(locusele)
                        Keyintegron.addInttag(str(line).split('\t')[3].split(';')[locusele],
                                              'IntI1')
                    if 'SulI' in key:
                        # set sulI loci
                        Keyintegron.addSul1(locusele)
                        Keyintegron.addSul1locus(str(line).split('\t')[3].split(';')[locusele])
                locusele += 1
            # Locus1 for the integron min locus, Locus2 for the integron max locus
            Locus1 = float(str(line).split('\t')[3].split(';')[locusmin])
            Locus2 = float(str(line).split('\t')[3].split(';')[locusmax])
            Keyintegron.setL1(min(Locus1,Locus2))
            Keyintegron.setL2(max(Locus1,Locus2))
            # set integron structure by integron type and number of elements
            Keyintegron.setStructure(len(str(line).split('\t')[3].split(';')))
            # this integron has not been output
            Keyintegron.setOutput()
            # delete Type E integrons with only one attC
            Keyintegron.deleteKey()
            # delete Type D with only one attC locus (two attC sites may have the same loci)
            if Keyintegron.delete!='TRUE':
                if gbff_path != 'None':
                    # optional: extraction of integrase annotation and sequences from genbank
                    # integrase annotation will be added into the class Integron
                    print "Extract integrase annotation from genbank files!"
                    extract_integrase_gbff(list_gbff, input_file, Keyintegron)
                # extract ORFs
                print "Extract integron elements and ORFs!"
                extract_ORFs(input_file, Keyintegron, AA_format)
                # extract integron and integrase sequences
                print "Extract integron sequences!"
                extract_Seqs(input_file, Keyintegron)
                if gbff_path != 'None':
                    # optional: extraction of integron annotation and sequences from genbank
                    print "Extract ORF annotation from genbank files!"
                    extract_gbff(list_gbff, input_file, Keyintegron)


def write_integrase(key, LocusORF1, LocusORF2, i, f, record, end, file):
    # output integrase
    if end == 'None':
        for locusint in key.intlocus:
            # find the integrase annotation for the ORF
            if locusint-10 <= (LocusORF1+LocusORF2)/2 <= locusint+10:
                # find the right integrase
                # integrase name: integron_id (key.name) + ':' + number_of_ORF_on_integron (i)
                # output integrase annotation 1. key.intanno2.get(locusint) is previously imported from genbank
                # output integrase annotation 2. record.description
                if 'IntI1' in str(key.inttag2.get(locusint)).replace('[', '').replace(']', ''):
                    # annotated as intI1
                    f.write(
                        '>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(locusint)) + '\t'
                        + str(record.description) + '\n')
                    f.write(str(record.seq) + '\n')
                    return 'IntI1'
                else:
                    # annotated as integrases
                    fnotA = open(os.path.join(str(resultdir) + '/Other/Integrase_seqs',
                                              str(file).replace(str(input_format), '.txt')),
                                 'a')
                    fnotA.write(
                        '>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(locusint)) + '\t'
                        + str(record.description) + '\n')
                    fnotA.write(str(record.seq) + '\n')
                    fnotA.close()
                    return 'Integrase'
        return 'ORF'
    else:
        try:
            f.write('>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(key.intanno1[int(end)])) +
                    '\t' + str(record.description) + '\n')
        except IndexError:
            f.write('>' + str(key.name) + ':0' + '\tNone\t' + str(record.description) + '\n')
        f.write(str(record.seq) + '\n')
        return 'None'


def write_AA(key, LocusORF1, LocusORF2, i, f, record, Tag):
    # output ORF
    # name of ORF: integron_id (key.name) + ':' + number_of_ORF_on_integron (i)
    if len(record.seq)==0:
        flog.write('empty ORF: '+str(key.name)+ ':' + str(i) +'\t'+str(record.id)+'\n')
    else:
        Sequence = record.seq
        f.write('>' + str(key.name) + ':' + str(i) + '\t'+str(Tag)+'\t' + str(LocusORF1) + '\t' +
             str(LocusORF2) + '\t' + str(record.description) + '\n')
        f.write(str(Sequence) + '\n')


def write_structure(key,Locus,i,f,Tag,temp):
    # output integron structure
    if temp!= 'None':
        # output f4temp
        f.write(temp)
    elif Tag == 'attC':
        # output attC, name of attC: integron_id (key.name) + ':' + 'attC'
        f.write(str(key.type) + '\t' + str(key.name) + ':' +str(Tag)+ '\t'+str(Tag)+'\t' + str(Locus) + '\n')
    elif Tag != 'None':
        # output ORF, integrase, intI1, and sul1
        f.write(str(key.type) + '\t' + str(key.name) + ':' + str(i) + '\t'+str(Tag)+'\t' + str(Locus) + '\n')
        return ''


def compare_structure(key, Locusmean, i, f, f4temp, Tag):
    # output the structure of ORF and empty the temp list of structure
    if f4temp != '':
        write_structure(key, Locusmean, i, f, Tag, f4temp)
    # format output for ORF
    f4temp = str(key.type) + '\t' + str(key.name) + ':' + str(i) + '\t' + str(Tag) + '\t' + str(
        Locusmean) + '\n'
    # adjust the structure of ORFs with attC
    try:
        # a temp list to remove attC that has been output
        templist = copy.deepcopy(key.attclocus)
        locusattc = key.attclocus[0]
        for locusattc in key.attclocus:
            if Locusmean > locusattc:
                # write attC sites ahead of ORF
                write_structure(key, locusattc, i, f, 'attC', 'None')
                templist.remove(locusattc)
            elif Locusmean <= locusattc:
                # when no attC is found to be in front of the ORF, write ORF
                write_structure(key, Locusmean, i, f, Tag, f4temp)
                f4temp = ''
                break
        key.attclocus = templist
    except IndexError:
        # no attC left on the integron
        # output the structure of ORF and empty the temp list of structure
        pass
    return f4temp


def write_forward(key,Beforetemp1,Beforelength1,f1,f3,f4,file):
    # write one extra forward ORF that is not in the range of integron
    if key.structure[0] == -1 and str(Beforetemp1) != '':
        # for integrons with no integrase and sul1 at the left extremity
        # output integrase, not like the extra forward ORF is an integrase
        Tag = write_integrase(key, Beforelength1[0], Beforelength1[1], '-1', f3, Beforetemp1,
                              'None', file)
        # output ORF
        write_AA(key, Beforelength1[0], Beforelength1[1], '-1', f1, Beforetemp1, Tag)
        # output integron structure
        write_structure(key, Beforelength1[2], '-1', f4, Tag, 'None')
        # add the forward ORF loci to the integorn range
        key.changeL1(Beforelength1[0])
    return ''


def extract_ORFs(file,key,AA_format):
    # extract ORFs on the integron
    # set output folder
    # key stored one integron
    if key.type not in ['A','B']:
        f1 = open(os.path.join(str(resultdir)+'/Other/ORFs',str(file).replace(str(input_format),'.AA')),'a')
        f3 = open(os.path.join(str(resultdir) + '/Other/Integrase_seqs', str(file).replace(str(input_format), '.AA')), 'a')
        f4 = open(os.path.join(str(resultdir) + '/Other/Integron_structure', str(file).replace(str(input_format), '.txt')), 'a')
    else:
        f1 = open(os.path.join(str(resultdir) + '/ClassI/ORFs', str(file).replace(str(input_format), '.AA')), 'a')
        f3 = open(os.path.join(str(resultdir) + '/ClassI/Integrase_seqs', str(file).replace(str(input_format), '.AA')),
                  'a')
        f4 = open(os.path.join(str(resultdir) + '/ClassI/Integron_structure', str(file).replace(str(input_format), '.txt')),
                  'a')
    # set input ORF file
    f2 = SeqIO.parse(os.path.join(AA_path, str(file).replace(str(input_format), args.o)) + '.add', 'fasta')
    # start extracting ORFs
    if key.output != 'TRUE':
        # integrons that have not been output
        # store forward and afterward extra ORFs
        Beforetemp1 = '' # store forward ORF
        Beforelength1 =[] # store forward ORF loci
        Aftertemp = -1 # whether to output extra afterward ORF
        i = 1 # the number of ORFs
        Locusmean=0 # ORF average locus
        f4temp='' # temp file for integron structure output
        for record in f2:
            contigname = key.name.replace(args.f, args.o).split(':')[0]
            if contigname in record.id:
                # for each ORF, set the locusORF
                if AA_format == 2:
                    # prodigal output
                    LocusORF1 = float(str(record.description).split(' # ')[1])
                    LocusORF2 = float(str(record.description).split(' # ')[2])
                elif AA_format == 1:
                    # gbff parse output
                    LocusORF1 = float(str(record.description).split(' ')[-2])
                    LocusORF2 = float(str(record.description).split(' ')[-1])
                # ORF average locus
                Locusmeanold = Locusmean
                Locusmean = (LocusORF1 + LocusORF2) / 2
                if Locusmean != Locusmeanold:
                    if key.l1 > LocusORF2: # the ORF is still ahead of integron range
                        # update forward ORFs
                        Beforetemp1 = record
                        Beforelength1 = [min(LocusORF1, LocusORF2), max(LocusORF1, LocusORF2),
                                         (LocusORF1 + LocusORF2) / 2]
                        # set whether to output extra afterward ORF
                        Aftertemp = key.structure[-1]
                    elif key.l1 <= Locusmean <= key.l2:
                        # the ORF is in the range of integron loci
                        # Three conditions are used to completely extract all possible ORFs
                        key.doOutput()
                        # integron has been output
                        Tag = 'ORF'
                        # add ORF loci to the integorn range
                        key.changeL1(min(LocusORF1, LocusORF2))
                        key.changeL2(max(LocusORF1, LocusORF2))
                        if i == 1:
                            # first ORF on the integron
                            # writing the forward ORFs
                            Beforetemp1 = write_forward(key, Beforetemp1,  Beforelength1,  f1, f3,
                                                        f4, file)
                        # start writing ORFs in the middle
                        if any(locussulI - 10 <= Locusmean <= locussulI + 10 for locussulI in key.sul1locus):
                            # if the ORF is sul1
                            Tag = 'SulI'
                        else:
                            # if the ORFs is integrase, intI1 or ORF
                            Tag = write_integrase(key, LocusORF1, LocusORF2, i, f3, record, 'None',file)
                        # output ORF
                        write_AA(key, LocusORF1, LocusORF2, i, f1, record, Tag)
                        # adjust the structure of ORFs with attC
                        f4temp = compare_structure(key, Locusmean, i, f4, f4temp, Tag)
                        i += 1
                        Aftertemp = key.structure[-1]
                    if Locusmean > key.l2 and Aftertemp >= 0 and str(Beforetemp1) == '':
                        # write all afterward ORFs and attC
                        if f4temp != '':
                            f4temp = write_structure(key, Locusmean, '+1', f4, 'None', f4temp)
                            # write all remaining ORFs after attC
                        try:
                            # write all remaining attC sites
                            locusattc = key.attclocus[0]
                            for locusattc in key.attclocus:
                                write_structure(key, locusattc, '+1', f4, 'attC', 'None')
                            key.attclocus = []
                            break
                        except IndexError:
                            pass
                        if Aftertemp == 0:
                            # sulI or Integrase at the right extremity
                            break
                        if Aftertemp == 1:
                            # write one extra ORF
                            write_structure(key, Locusmean, '+1', f4, 'None', f4temp)
                            write_AA(key, LocusORF1, LocusORF2, '+1', f1, record, 'ORF')
                            write_structure(key, Locusmean, '+1', f4, 'ORF', 'None')
                            break
                        #extend integron range
                        key.changeL2(max(LocusORF1, LocusORF2))
            elif Aftertemp >= 0 and str(Beforetemp1) == '':
                # integron structure at the end of the last sequence
                # write all afterward ORFs and attC
                if f4temp != '':
                    f4temp = write_structure(key, Locusmean, '+1', f4, 'None', f4temp)
                    # write all remaining ORFs after attC
                try:
                    # write all remaining attC sites
                    locusattc = key.attclocus[0]
                    for locusattc in key.attclocus:
                        write_structure(key, locusattc, '+1', f4, 'attC', 'None')
                    key.attclocus = []
                    break
                except IndexError:
                    pass
                Aftertemp = -1
                break
    f1.close()
    f3.close()
    f4.close()


def extract_Seqs(file,key):
    # extract integron sequences
    if key.type in ['A','B']:
        f1 = open(os.path.join(str(resultdir) + '/ClassI/Integron_seqs', str(file).replace(str(input_format), '.fasta')),
                  'a')
    else:
        f1 = open(os.path.join(str(resultdir) + '/Other/Integron_seqs', str(file).replace(str(input_format), '.fasta')), 'a')
    f3 = open(os.path.join(str(resultdir), 'Contigs_with_integrons.fasta'),
                  'a')
    # input target sequence file
    f2 = SeqIO.parse(os.path.join(file_path, str(file).replace(str(input_format), args.f)), 'fasta')
    for record in f2:
        if key.output == 'FALSE':
            flog.write('Failed output integron: '+str(key.name)+'\t'+str(key.l1)+'\t'+str(key.l2)+'\n')
        if record.id in key.name:
            # extract [integron range - 5, integron range + 5]
            f1.write('>' + str(key.name) + '\t' +str(key.type)+'\t'+ str(record.description) + '\n')
            f1.write(str(record.seq)[max(int(key.newl1 - 5),0):min(int(key.newl2 +5),len(record.seq)-1)] + '\n')
            f3.write('>' + str(key.name) + '\t' + str(key.type) + '\t' + str(record.description) + '\n')
            f3.write(str(record.seq) + '\n')
    f1.close()
    f3.close()


def extract_gbff(gbff_list, file, key):
    # optional: extract genbank annotation
    for gbfffile in gbff_list:
            if key.type in ['A', 'B']:
                f1 = open(
                    os.path.join(str(resultdir) + '/ClassI/Genbank_annotation', str(file).replace(str(input_format), '.annotation')),
                    'a')
            else:
                f1 = open(
                    os.path.join(str(resultdir) + '/Other/Genbank_annotation',
                                 str(file).replace(str(input_format), '.annotation')),
                    'a')
            for record in SeqIO.parse(gbfffile, "genbank"):
                if record.id in str(key.name):
                    # find the target sequence
                    for feature in record.features:
                        annotation = []  # annotation of one gene
                        start = float(feature.location.start.position)
                        end = float(feature.location.end.position)
                        meanlocus = (start + end) / 2
                        if key.newl1 <= meanlocus <= key.newl2:
                            annotation.append(str(feature.type) + ';' +
                                              str(feature.qualifiers.get('product', 'None')[0]) + ';' +
                                              str(feature.qualifiers.get("gene", 'None')[0]) + ';' +
                                              str(feature.qualifiers.get("note", 'None')[0]))
                            locus_tag = feature.qualifiers.get('locus_tag', 'None')[0]
                            f1.write(str(key.name) + '\t' + str(start) + "\t" + str(end) + "\t" + str(
                                locus_tag) + "\t" + str(annotation) + "\n")
            f1.close()


def extract_integrase_gbff(gbff_list, file, key):
    # optional: extraction of integrase annotation and sequences from genbank
    if key.type in ['A', 'C', 'B']:
        # for Types A, B, C integrons with integrases
        for gbfffile in gbff_list:
                for record in SeqIO.parse(gbfffile, "genbank"):
                    if record.id in str(key.name):
                        # find the target sequence that harboured the integrases
                        for feature in record.features:
                            # annotation of ORF
                            start = float(feature.location.start.position)
                            end = float(feature.location.end.position)
                            for locusint in key.intlocus:
                                if min(start, end) <= locusint <= max(start, end):
                                    # find the right ORF of integrases
                                    annotation = [str(feature.type) + ';' +
                                                  str(feature.qualifiers.get('product', 'None')[0]) + ';' +
                                                  str(feature.qualifiers.get("gene", 'None')[0]) + ';' +
                                                  str(feature.qualifiers.get("note", 'None')[0])]
                                    key.addIntanno(locusint, annotation)
                                    # add integrase annotation


################################################### Programme #######################################################
flog = open('Integron_Extraction.log', 'a')
# main programme body
for file_name in list_fasta:
    try:
        in_dir, input_file = os.path.split(file_name)
        # genbank annotation extraction, optional
        if gbff_path != 'None':
            # genbank accession number, e.g., GCA_900321765.1
            accession = os.path.splitext(input_file.split(input_format)[0])[0]
            list_gbff = glob.glob(os.path.join(gbff_path, '*' + accession + '*'))
        # load ORFs
        # AA_format, 1 for genbank parsing; 2 for prodigal predicting
        AA_format = ORFformat.get(input_file.replace(input_format,args.f), 'None')
        if AA_format == 'None':
            print 'No ORF format information for ' +\
                  os.path.join(file_path, input_file.replace(input_format,args.f)) + '!'
            flog.write('No ORF format information for ' +
                  os.path.join(file_path, input_file.replace(input_format,args.f)) + '!\n')
        else:
            # main function of integron extraction
            int_integron(input_file,AA_format)
    except IOError:
        flog.write('Files were missing?\n')
        flog.write(str(input_file)+'\n')
# summarize output results into single files
cmd='cat '+str(resultdir) + '/Other/Genbank_annotation/* > '+str(resultdir) + '/all.Genbank_annotation.txt\n'+\
'cat '+str(resultdir) + '/Other/Integrase_seqs/* > '+str(resultdir) + '/all.Integrase_seqs.fasta\n'+\
'cat '+str(resultdir) + '/Other/Integron_structure/* > '+str(resultdir) + '/all.Integron_structure.txt\n'+\
'cat '+str(resultdir) + '/Other/ORFs/* > '+str(resultdir) + '/all.ORFs.fasta\n'+\
'cat '+str(resultdir) + '/Other/Integron_seqs/* > '+str(resultdir) + '/all.Integron_seqs.fasta\n'+\
'cat '+str(resultdir) + '/ClassI/Integrase_seqs/* > '+str(resultdir) + '/all.Integrase_seqs.ClassI.fasta\n'+\
'cat '+str(resultdir) + '/ClassI/Integron_structure/* > '+str(resultdir) + '/all.Integron_structure.ClassI.txt\n'+\
'cat '+str(resultdir) + '/ClassI/ORFs/* > '+str(resultdir) + '/all.ORFs.ClassI.fasta\n'+\
'cat '+str(resultdir) + '/ClassI/Integron_seqs/* > '+str(resultdir) + '/all.Integron_seqs.ClassI.fasta\n'+\
'cat '+str(resultdir) + '/ClassI/Genbank_annotation/* > '+str(resultdir) + '/all.Genbank_annotation.ClassI.txt\n'
os.system(cmd)
print 'Finished integron extraction!'
flog.write('Finished integron extraction!\n')
flog.close()