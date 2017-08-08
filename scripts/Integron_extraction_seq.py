import os
import glob
import argparse
from Bio import SeqIO

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-input",
                    help="A single file name or dir name (default: result)", type=str, default='result',metavar='your file or dir')
parser.add_argument("--input_format",
                    help="The input format for files in a dir", type=str, default='.Integron.txt',metavar='.Integron.txt')
parser.add_argument('--resultdir',
                    default="Integron", action='store', type=str, metavar='Integron',
                    help="Set the result directory for integron sequences and annotation (default: Integron)")
parser.add_argument('--AAdir',
                    default=".", action='store', type=str, metavar='.',
                    help="Set the directory of ORFs files (default: .)")
parser.add_argument('--AA_format',
                    help="The input AA type, eg: 0 for a Genbank parsed file; 1 for a prodigal predicted file",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=0, type=int)
parser.add_argument('--gbffdir',
                    default="gbff_copy", action='store', type=str, metavar='gbff_copy',
                    help="Set the directory of Genbank files (default: gbff_copy)\nNo gbff input: None")
parser.add_argument("--fasta",
                    help="format of fasta sequences", type=str, default='.fa', metavar='.fa')
parser.add_argument("--orf",
                    help="format of ORFs", type=str, default='.faa', metavar='.faa')
################################################## Definition ########################################################
args = parser.parse_args()
AA_path = os.path.abspath(args.AAdir)
if args.gbffdir!='None':
    gbff_path = os.path.abspath(args.gbffdir)
else:
    gbff_path=args.gbffdir
AA_format=args.AA_format
input_path = os.path.abspath(args.input)
in_dir = os.path.abspath(input_path)
list_fasta = glob.glob(os.path.join(in_dir,'*'+args.input_format))
try:
    os.mkdir(args.resultdir)
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/Other')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir)+'/Other/ORFs')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/Other/Seqs')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/Other/Intseqs')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/Other/Anno')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/Other/Intstruc')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/ClassI')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/ClassI/ORFs')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/ClassI/Seqs')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/ClassI/Intseqs')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/ClassI/Anno')
except OSError:
    pass
try:
    os.mkdir(str(args.resultdir) + '/ClassI/Intstruc')
except OSError:
    pass

################################################### Function #########################################################
__metaclass__=type
class Integron:
    'name,type,int,sulI,l1,l2,structure'
    def setName(self,name):
        self.name=name
    def setType(self,type):
        self.type=type
    def setInt(self,Int):
        self.int=[]
    def addInt(self,Int):
        self.int.append(float(Int))
    def setL1(self,Locus1): #locus of integron
        self.l1=float(Locus1)
        self.newl1 = self.l1
    def changeL1(self,Locus1): #additional locus of potential intergrase
        if self.l1>float(Locus1):
            self.newl1=float(Locus1)
    def setL2(self,Locus2): #locus of integron
        self.l2=float(Locus2)
        self.newl2 = self.l2
    def changeL2(self,Locus2): #additional locus of potential intergrase
        if self.l2<float(Locus2):
            self.newl2 = float(Locus2)
    def setInt(self):
        self.int=[]
    def addInt(self,Int):
        self.int.append(float(Int))
    def setIntlocus(self):
        self.intlocus=[]
    def addIntlocus(self,Int): #no duplicated locus
        if float(Int) not in self.intlocus:
            self.intlocus.append(float(Int))
    def setattClocus(self):
        self.attclocus=[]
    def addattClocus(self,attc): #no duplicated locus
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
    def addSul1locus(self,Sul1l): #no duplicated locus
        if float(Sul1l) not in self.sul1locus:
            self.sul1locus.append(float(Sul1l))
    def setStructure(self,length): #integron classification
        Before=-1
        After=1
        if self.type == 'A':
            if max(self.int)==(length-1):
                After=0
            if min(self.int)==0:
                Before=0
            if max(self.sul1)==(length-1):
                After=0
            if min(self.sul1)==0:
                Before=0
        elif self.type=='C' or self.type=='B':
            if max(self.int)==(length-1):
                After=0
            if min(self.int)==0:
                Before=0
        elif self.type=='D':
            if max(self.sul1)==(length-1):
                After=0
            if min(self.sul1)==0:
                Before=0
        self.structure = [Before, After]
    def setOutput(self):
        self.output = 'FALSE'
    def doOutput(self):
        self.output = 'TRUE'
    def deleteKey(self):
        if self.type=='E' and len(self.attclocus)==1:
            self.delete='TRUE'
        else:
            self.delete='FALSE'



def int_integron(file,AA_format):
    for line in open(os.path.join(in_dir,file),'rb'):
        if ':' in str(line).split('\t')[1]: #set up Integron
            Keyintegron=Integron()
            Keyintegron.setName(str(line).split('\t')[1])
            print Keyintegron.name
            Keyintegron.setType(str(line).split('\t')[0])
            Keyintegron.setattClocus()
            locusmin = locusmax = locusattc = locusstart =0
            for key in str(line).split('\t')[4].split(';'):
                if 'attC' in key or 'SulI' in key or 'Integrase' in key or 'IntI1' in key: #locus setup
                    locusmax = locusattc
                    if locusstart==0:
                        locusmin=locusattc
                        locusstart=1
                    if 'attC' in key: #set attClocus
                        Keyintegron.addattClocus(str(line).split('\t')[3].split(';')[locusattc])
                locusattc += 1
            Locus1=float(str(line).split('\t')[3].split(';')[locusmin])
            Locus2=float(str(line).split('\t')[3].split(';')[locusmax])
            Keyintegron.setL1(min(Locus1,Locus2))
            Keyintegron.setL2(max(Locus1,Locus2))
            Keyintegron.setInt()
            Keyintegron.setIntlocus()
            Keyintegron.setIntanno()
            Keyintegron.setInttag()
            Keyintegron.setSul1()
            Keyintegron.setSul1locus()
            if Keyintegron.type in ['A','B','C','D']: #setup integrase and sulI location
                Keyintegron.setIntlocus()
                locusint=0
                for key in str(line).split('\t')[4].split(';'):
                    if 'Integrase' in key or 'IntI1' in key:
                        Keyintegron.addIntlocus(str(line).split('\t')[3].split(';')[locusint])
                        Keyintegron.addInt(locusint)
                        Keyintegron.addInttag(str(line).split('\t')[3].split(';')[locusint], \
                                              str(line).split('\t')[4].split(';')[locusint].split(' ')[1])
                    if 'SulI' in key:
                        Keyintegron.addSul1(locusint)
                        Keyintegron.addSul1locus(str(line).split('\t')[3].split(';')[locusint])
                    locusint += 1
            else:
                Keyintegron.addInt(-1)
                Keyintegron.addSul1(-1)
            Keyintegron.setStructure(len(str(line).split('\t')[3].split(';')))
            Keyintegron.setOutput()
            Keyintegron.deleteKey()  # delete Type D with only one attC locus (two attC sites may have the same locus)
            if Keyintegron.delete!='TRUE':
                if gbff_path != 'None': #extraction of integron annotation and sequences
                    extract_integrase_gbff(list_gbff, input_file, Keyintegron)
                extract_ORFs(input_file, Keyintegron,AA_format)
                extract_Seqs(input_file, Keyintegron)
                if gbff_path != 'None': #extraction of integron annotation and sequences
                    extract_gbff(list_gbff, input_file, Keyintegron)


def write_integrase(key,LocusORF1,LocusORF2,i,f,record,end,file):
    if end =='None':
        for locusint in key.intlocus:  # write integrase AA
            if  locusint-5<= (LocusORF1+LocusORF2)/2 <= locusint+5:
                if 'IntI1' in str(key.inttag2.get(locusint)).replace('[', '').replace(']', ''):
                    f.write(
                        '>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(locusint)) + '\t' \
                        + str(record.description) + '\n')
                    f.write(str(record.seq) + '\n')
                else:
                    fnotA = open(os.path.join(str(args.resultdir) + '/Other/Intseqs',
                                           str(file).replace(str(args.input_format), '.txt')),
                              'ab')
                    fnotA.write(
                        '>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(locusint)) + '\t' \
                        + str(record.description) + '\n')
                    fnotA.write(str(record.seq) + '\n')
                    fnotA.close()
                flog.write('>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(locusint)) + '\t' \
                    + str(record.description) + '\n')
                return 'Integrase'
        return 'ORF'
    else:
        try:
            f.write('>' + str(key.name) + ':' + str(i) + '\t' + str(key.intanno2.get(key.intanno1[int(end)])) + \
                     '\t' + str(record.description) + '\n')
        except IndexError:
            f.write('>' + str(key.name) + ':0' + '\tNone\t' + str(record.description) + '\n')
        f.write(str(record.seq) + '\n')
        return 'None'


def write_AA(key,LocusORF1,LocusORF2,i,f,record,Tag,file):
    if len(record.seq)==0:#for empty ORFs, no output
        #f5 = SeqIO.parse(os.path.join(AA_path, str(file).replace(str(args.input_format), '.fa')), 'fasta')
        flog.write('empty ORF: '+str(key.name)+ ':' + str(i) +'\t'+str(record.id)+'\n')
        #for record in f5:
        #    Sequence=str(record.seq)[int(LocusORF1 - 1):int(LocusORF2 - 1)]
    else:
        Sequence=record.seq
        f.write('>' + str(key.name) + ':' + str(i) + '\t'+str(Tag)+'\t' + str(LocusORF1) + '\t' +
             str(LocusORF2) + '\t' + str(record.description) + '\n')
        f.write(str(Sequence) + '\n')


def write_structure(key,Locus,i,f,Tag,temp):
    if Tag=='attC':
        f.write(str(key.type) + '\t' + str(key.name) + ':' +str(Tag)+ '\t'+str(Tag)+'\t' + str(Locus) + '\n')
    elif Tag !='None':
        f.write(str(key.type) + '\t' + str(key.name) + ':' + str(i) + '\t'+str(Tag)+'\t' + str(Locus) + '\n')
    elif temp!='None':
        f.write(temp)
        return ''


def compare_structure(key,Locusmean,i,f,f4temp,Tag):
    try:
        locusattcold = key.attclocus[0]  # attC and ORF structure and locus
        templist = key.attclocus
        if len(key.attclocus)>=2:
            for locusattc in key.attclocus[1:]:  # at least two attC sites
                if Locusmean > locusattc:  # write before attC sites
                    write_structure(key, locusattcold, i, f, 'attC','None')
                    templist.remove(locusattcold)
                    if str(key.name)  + ':' + str(i) + '\t'+str(Tag) not in f4temp:
                        f4temp = f4temp + str(key.type) + '\t' + str(key.name) + ':' + str(i) + '\t'+str(Tag)+'\t' + str(
                            Locusmean) + '\n'
                elif locusattcold <= Locusmean <= locusattc:  # write before attC sites
                    write_structure(key,locusattcold, i, f, 'attC','None')
                    if str(key.name)  + ':' + str(i) + '\t'+str(Tag) not in f4temp:
                        f4temp = f4temp + str(key.type) + '\t' + str(key.name) + ':' + str(i) + '\t'+str(Tag)+'\t' + str(
                            Locusmean) + '\n'
                    f.write(f4temp)
                    f4temp = ''
                    templist.remove(locusattcold)
                    break
                elif Locusmean < locusattcold:  # write ORF
                    f.write(f4temp)
                    f4temp=''
                    write_structure(key, Locusmean, i, f, Tag,'None')
                    break
                locusattcold = locusattc
            key.attclocus = templist
        elif len(key.attclocus)==1:
            locusattc = key.attclocus[0]
            if locusattc <= Locusmean:  # order of attC and ORF
                f.write(f4temp)
                write_structure(key, locusattcold, i, f, 'attC', 'None')
                write_structure(key, Locusmean, i, f, Tag, 'None')
                key.attclocus.remove(locusattcold)
            else:
                f.write(f4temp)
                write_structure(key, Locusmean, i, f, Tag, 'None')
                #write_structure(key, locusattcold, i, f, 'attC', 'None')
            f4temp = ''
    except IndexError:
        f.write(f4temp)
        write_structure(key, Locusmean, i, f, Tag, 'None')
        f4temp = ''
    return f4temp


def write_forward(key,Beforetemp1,Beforelength1,f1,f3,f4,file):
    if str(Beforetemp1) != '':  # write forward ORFs
        if key.structure[0] == -1 and str(Beforetemp1) != '':
            write_AA(key, Beforelength1[0], Beforelength1[1], '-1', f1, Beforetemp1, 'ORF',file)
            write_integrase(key, Beforelength1[0], Beforelength1[1], '-1', f3, Beforetemp1,
                            'None',file)
            write_structure(key, Beforelength1[2], '-1', f4, 'ORF','None')
            key.changeL1(Beforelength1[0])
    return 0


def extract_ORFs(file,key,AA_format):
    if key.type!='A':
        f1 = open(os.path.join(str(args.resultdir)+'/Other/ORFs',str(file).replace(str(args.input_format),'.AA')),'ab')
        f3 = open(os.path.join(str(args.resultdir) + '/Other/Intseqs', str(file).replace(str(args.input_format), '.AA')), 'ab')
        f4 = open(os.path.join(str(args.resultdir) + '/Other/Intstruc', str(file).replace(str(args.input_format), '.txt')), 'ab')
    else:
        f1 = open(os.path.join(str(args.resultdir) + '/ClassI/ORFs', str(file).replace(str(args.input_format), '.AA')), 'ab')
        f3 = open(os.path.join(str(args.resultdir) + '/ClassI/Intseqs', str(file).replace(str(args.input_format), '.AA')),
                  'ab')
        f4 = open(os.path.join(str(args.resultdir) + '/ClassI/Intstruc', str(file).replace(str(args.input_format), '.txt')),
                  'ab')
    f2 = SeqIO.parse(os.path.join(AA_path, str(file).replace(str(args.input_format), args.orf)), 'fasta')
    if key.output != 'TRUE': #not extracted integrons
        Beforetemp1 = '' #forward and afterward ORFs temp storing
        Beforelength1 =[]
        Aftertemp = -1
        i = 1 #ORFs number/locus
        Locusmean=0
        f4temp=''
        for record in f2:# for each ORF
            temp='_'.join((str(record.id).split('_')[0:-1]))
            if str(key.name).split(':')[0] == temp:
                print record.id
                if AA_format==1: #prodigal output
                    LocusORF1 = float(str(record.description).split(' # ')[1])
                    LocusORF2 = float(str(record.description).split(' # ')[2])
                else: #gbff parse output
                    LocusORF1 = float(str(record.description).split(' ')[-2])
                    LocusORF2 = float(str(record.description).split(' ')[-1])
                Locusmeanold = Locusmean
                Locusmean = (LocusORF1 + LocusORF2) / 2
                if Locusmean != Locusmeanold:
                    if key.l1 > LocusORF2:  # set forward ORFs
                        Beforetemp1 = record
                        Beforelength1 = [min(LocusORF1, LocusORF2), max(LocusORF1, LocusORF2),
                                         (LocusORF1 + LocusORF2) / 2]
                        Aftertemp = key.structure[-1]
                    elif key.l1 == key.l2 and str(Beforetemp1)!='':  # B with only one attC
                        key.doOutput()
                        Beforetemp1 =write_forward(key, Beforetemp1,  Beforelength1,  f1, f3, f4,file)
                        write_structure(key, key.attclocus[0], '-1', f4, 'attC', 'None')
                        key.attclocus.remove(key.attclocus[0])
                        Aftertemp = key.structure[-1]
                    elif key.l1 <= LocusORF2 <= key.l2 or key.l1 <= LocusORF1 <= key.l2 or key.l1 <= Locusmean <= key.l2:
                        # ORF between two attCs or one attC + SulI
                        # Locusmean, LocusORF1 and LocusORF2 are all used here to completely extract all possible ORFs
                        key.doOutput()
                        Tag = 'ORF'
                        key.changeL1(min(LocusORF1, LocusORF2))
                        key.changeL2(max(LocusORF1, LocusORF2))
                        if i == 1:  # start writing forward ORFs
                            Beforetemp1 = write_forward(key, Beforetemp1,  Beforelength1,  f1, f3,
                                                        f4,file)
                        # start writing ORFs in the middle
                        if any(locussulI - 5 <= (LocusORF1 + LocusORF2) / 2 <= locussulI + 5 for locussulI in key.sul1locus):
                            Tag = 'SulI'
                        else:
                            Tag = write_integrase(key, LocusORF1, LocusORF2, i, f3, record, 'None',file)
                        write_AA(key, LocusORF1, LocusORF2, i, f1, record, Tag,file)
                        f4temp = compare_structure(key, Locusmean, i, f4, f4temp,Tag)
                        i += 1
                        Aftertemp = key.structure[-1]
                    elif Aftertemp >= 0 and str(Beforetemp1)=='':  # write all afterward ORFs
                        if f4temp!='':
                            f4temp=write_structure(key, Locusmean, '+1', f4, 'None', f4temp) #write all remaining ORFs after all attC sites
                        if Aftertemp==0: #sulI or Integrase
                            try:#write all remaining attC sites
                                for locusattc in key.attclocus:
                                    write_structure(key, locusattc, '+1', f4, 'attC','None')
                                key.attclocus = []
                                break
                            except KeyError:
                                break
                        if Aftertemp == 1:#attC ends
                            write_structure(key, Locusmean, '+1', f4, 'None', f4temp)
                            write_AA(key, LocusORF1, LocusORF2, '+1', f1, record, 'ORF',file)
                            try:#write all remaining attC sites
                                for locusattc in key.attclocus:
                                    write_structure(key, locusattc, '+1', f4, 'attC','None')
                                key.attclocus=[]
                            except KeyError:
                                pass
                            write_structure(key, Locusmean, '+1', f4, 'ORF', 'None')
                        Aftertemp -= 1
                        key.changeL2(max(LocusORF1, LocusORF2))
    f1.close()
    f3.close()
    f4.close()



def extract_Seqs(file,key):
    if key.type=='A':
        f1 = open(os.path.join(str(args.resultdir) + '/ClassI/Seqs', str(file).replace(str(args.input_format), '.fasta')),
                  'ab')
    else:
        f1 = open(os.path.join(str(args.resultdir) + '/Other/Seqs', str(file).replace(str(args.input_format), '.fasta')), 'ab')
    f2 = SeqIO.parse(os.path.join(AA_path, str(file).replace(str(args.input_format), args.fasta)), 'fasta')
    for record in f2:
        if str(key.name).split(':')[0] == record.id:
            print str(key.name).split(':')[0]
            print record.id
            if key.output == 'FALSE':
                flog.write('Failed output integron: '+str(key.name)+'\t'+str(key.l1)+'\t'+str(key.l2)+'\n')
            f1.write('>' + str(key.name) + '\t' +str(key.type)+'\t'+ str(record.description) + '\n')
            locus1=int(key.newl1 - 1)
            locus2=int(key.newl2 - 1)
            if int(key.newl1 - 1)<0:
                locus1=0
            if int(key.newl2 - 1)>(len(str(record.seq))-1):
                locus2=len(str(record.seq))-1
            f1.write(str(record.seq)[locus1:locus2] + '\n')
    f1.close()



def extract_gbff(gbff_list, file, key):
    for gbfffile in gbff_list:
        if str(file).split('_')[1] in gbfffile:
            if key.type=='A':
                f1 = open(
                    os.path.join(str(args.resultdir) + '/ClassI/Anno', str(file).replace(str(args.input_format), '.annotation')),
                    'ab')
            else:
                f1 = open(
                    os.path.join(str(args.resultdir) + '/Other/Anno',
                                 str(file).replace(str(args.input_format), '.annotation')),
                    'ab')
            i = 0
            for record in SeqIO.parse(gbfffile, "genbank"):
                if str(key.name).split(':')[0] == record.id:
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
    for gbfffile in gbff_list:
        if str(file).split('_')[1] in gbfffile:
            i=0
            for record in SeqIO.parse(gbfffile, "genbank"):
                if str(key.name).split(':')[0] == record.id:
                    for feature in record.features:  # annotation of one gene
                        start = float(feature.location.start.position)
                        end = float(feature.location.end.position)
                        if key.type in ['A', 'B-non-IntI1', 'B-IntI1']:
                            for locusint in key.intlocus:
                                if min(start, end) <= locusint <= max(start, end):
                                    annotation = [str(feature.type) + ';' +
                                                  str(feature.qualifiers.get('product', 'None')[0]) + ';' +
                                                  str(feature.qualifiers.get("gene", 'None')[0]) + ';' +
                                                  str(feature.qualifiers.get("note", 'None')[0])]
                                    key.addIntanno(locusint, annotation)  # integrase annotation


################################################### Programme #######################################################
#optimal input of args.AA_format
flog = open(os.path.join(str(args.AAdir), 'Integron_Extraction_Pla.log'), 'ab')
#main programme body
for file_name in list_fasta:
    try:
        in_dir, input_file = os.path.split(file_name)
        print file_name
        if gbff_path != 'None':
            list_gbff = glob.glob(os.path.join(gbff_path, '*' + str(input_file).split('_')[1] + '*'))
        int_integron(input_file,AA_format)
    except IOError:
        flog.write('Files were missing?\n')
        flog.write(str(input_file)+'\n')
    finally:
        flog.write('Cleaning Up...\n')
        del file_name
#annotation of ORFs
flog.close()
#merge files
cmd='cat '+str(args.resultdir) + '/Other/Anno/* > '+str(args.resultdir) + '/all.Anno.txt\n'+\
'cat '+str(args.resultdir) + '/Other/Intseqs/* > '+str(args.resultdir) + '/all.Intseqs.fasta\n'+\
'cat '+str(args.resultdir) + '/Other/Intstruc/* > '+str(args.resultdir) + '/all.Intstruc.txt\n'+\
'cat '+str(args.resultdir) + '/Other/ORFs/* > '+str(args.resultdir) + '/all.ORFs.fasta\n'+\
'cat '+str(args.resultdir) + '/Other/Seqs/* > '+str(args.resultdir) + '/all.Seqs.fasta\n'+\
'cat '+str(args.resultdir) + '/ClassI/Intseqs/* > '+str(args.resultdir) + '/all.Intseqs.ClassI.fasta\n'+\
'cat '+str(args.resultdir) + '/ClassI/Intstruc/* > '+str(args.resultdir) + '/all.Intstruc.ClassI.txt\n'+\
'cat '+str(args.resultdir) + '/ClassI/ORFs/* > '+str(args.resultdir) + '/all.ORFs.ClassI.fasta\n'+\
'cat '+str(args.resultdir) + '/ClassI/Seqs/* > '+str(args.resultdir) + '/all.Seqs.ClassI.fasta\n'+\
'cat '+str(args.resultdir) + '/ClassI/Anno/* > '+str(args.resultdir) + '/all.Anno.ClassI.txt\n'
os.system(cmd)