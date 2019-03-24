import os
import glob
import argparse
from Bio import SeqIO

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument('--ot',
                    help="Optional: to provide your own CDS or ORFs files, \
                                         please input the orf type or the method you used to extract orfs,\
                                         eg: 1 for genbank parsing; 2 for prodigal prediction",
                    metavar="a file containing target files and their specific ORF format",
                    action='store', default='Temp/all.orf.length', type=str)
parser.add_argument("--tx",
                    help="a file of taxonomy metadata (under your input folder)",
                    type=str, default='genbank_taxon.txt',
                    metavar='genbank_taxon.txt or None.norm')
parser.add_argument("--tc",
                    help="column number corresponding to the taxonomy, i.e., from phylum to strain",
                    type=str, default='3,9', metavar='start column number, end column number')
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result', metavar='I-VIP_result')


################################################## Definition ########################################################
args = parser.parse_args()
if args.tx != 'None.norm':
    tx_dir, tx_file = os.path.split(os.path.abspath(args.tx))
    # load pathogen list
    Pathogen=[]
    for lines in open('database/Pathogen.txt','r'):
        Pathogen.append(lines.split('\r')[0].split('\n')[0])
in_dir = os.path.abspath(args.i)
# load all annotation results
annotationdir = args.r + "/Integron"
list_file = [ "all.Integron_structure.ClassI.annotated.formated.txt",
            "all.Integron_structure.annotated.formated.txt"]
resultdir = args.r + "/Phylogram"
try:
    os.mkdir(resultdir)
except OSError:
    pass


################################################### Function #########################################################


def node_add(node,shape,label,labelcolor,color,length):
    f1=open(os.path.join(resultdir,str(anno_file).replace('.txt','.node.txt')),'a')
    f1.write(str(node)+'\t'+str(shape)+'\t'+str(label)+'\t'+str(labelcolor)+'\t'+str(color)+'\t'+str(length)+'\n')


def edge_add(node1,node2,edgetype,color,length):
    f1 = open(os.path.join(resultdir, str(anno_file).replace('.txt','.edge.txt')), 'a')
    f1.write(str(node1) + '\t' + str(node2) + '\t' + str(edgetype) + '\t' + str(color)+ '\t' + str(length) + '\n')


def pathogen(taxa):
    for taxaname in taxa:
        if any(str(key).lower() in taxaname.lower() for key in Pathogen):
            # for potential human pathogens, set the label color to red
            return 'red'
        else:
            pass
            # for non-pathogens, set the label color to black
    return 'black'


def network(line,nodetable,edgetable,Repli,Taxon):
    if Taxon != 'None.norm':
        try:
            # contruct the phylogenetic tree by adding taxonomy nodes and edges
            for i in range((c1 - 1), c2):
                if i == (c1 - 1):
                    # for phylum level taxonomy, connect it to the bacterial node
                    if 'Bacteria' + '-' + str(line).split('\t')[i] not in edgetable:
                        # edge color as the class taxonomy
                        edge_add('Bacteria', str(line).split('\t')[i], 'taxon', str(line).split('\t')[c1], 1)
                        edgetable.append('Bacteria' + '-' + str(line).split('\t')[i])
                # for higher level of taxonomy, connect it to its lower level
                node1 = str(line).split('\t')[i].split('\r')[0].split('\n')[0]
                if len(node1.split(' ')) >= 2:
                    # shorten the taxonomy label
                    labelshort = node1.split(' ')[-2] + ' ' + \
                                 node1.split(' ')[-1]
                else:
                    labelshort = node1
                if i == (c2 - 1):
                    # connect the last taxonomy to the integron
                    # set integron type
                    Integron_type = str(line).split('\t')[0]
                    Type_node = str(node1 + '_' + Integron_type)
                    # set the label color as red if the host is a pathogen, black for non-pathogen
                    labelcolor = pathogen(str(line).split('\t')[(c2 - 3):(c2)])
                    # add edge for the integron type
                    if str(node1) + '-' + str(Type_node) not in edgetable and node1 != '':
                        # edge color as the class taxonomy
                        edge_add(node1, Type_node, 'taxon', str(line).split('\t')[c1], 1)
                        edgetable.append(str(node1) + '-' + str(Type_node))
                    # add node for the integron type
                    if str(Type_node) not in nodetable:
                        node_add(Type_node, Integron_type,
                                 '', Integron_type
                                 , Integron_type, 50)
                        nodetable.append(Type_node)
                else:
                    # for the taxonomy from class to species level
                    labelcolor = 'black'
                if node1 not in nodetable:
                    # add the taxonomy on the phylogenetic tree
                    # node color as class
                    node_add(node1, 'taxon', labelshort, labelcolor, str(line).split('\t')[c1],
                             50)
                    nodetable.append(node1)
                if i <= (c2 - 2):
                    # for the taxonomy from class to species level, add the edges
                    # a higher taxonomy
                    node2 = str(line).split('\t')[i + 1].split('\r')[0].split('\n')[0]
                    # connect the node1 and node2 by adding the edge
                    if str(node1) + '-' + str(node2) not in edgetable and node1 != '':
                        edge_add(str(node1), str(node2), 'taxon', str(line).split('\t')[c1], 1)
                        edgetable.append(str(node1) + '-' + str(node2))
        except IndexError:
            print "Wrong input for --tc \nPlease input in the format of: start column number, end column number\n"
            print  "For example: 2,8\nProceed without taxonomy information\n"
            # no taxonomy information
            # connect the integron to the same node
            node1 = 'Phylogram'
            # set integron type
            Integron_type = str(line).split('\t')[0]
            Type_node = str(node1 + '_' + Integron_type)
            # add edge for the integron type
            if str(node1) + '-' + Type_node not in edgetable and node1 != '':
                # edge color as the class taxonomy
                edge_add(node1, Type_node, 'None', 'None', 1)
                edgetable.append(str(node1) + '-' + Type_node)
            # add node for the integron type
            if Type_node not in nodetable:
                node_add(Type_node, Integron_type,
                         '', Integron_type
                         , Integron_type, 50)
                nodetable.append(Type_node)
    else:
        # no taxonomy information
        # connect the integron to the same node
        node1 = 'Phylogram'
        # set integron type
        Integron_type = str(line).split('\t')[0]
        Type_node = str(node1 + '_' + Integron_type)
        # add edge for the integron type
        if str(node1) + '-' + Type_node not in edgetable and node1 != '':
            # edge color as the class taxonomy
            edge_add(node1, Type_node, 'None', 'None', 1)
            edgetable.append(str(node1) + '-' + Type_node)
        # add node for the integron type
        if Type_node not in nodetable:
            node_add(Type_node, Integron_type,
                     '', Integron_type
                     , Integron_type, 50)
            nodetable.append(Type_node)
    # add the integron to the phylogram
    # set integron direction
    Locus = 0
    # check the integrase and sul1 locus
    Intlocus = SulIlocus = 0
    for key in str(line).split('\t')[cint].split(';'):
        Locus += 1
        if key !='':
            try:
                Type=key.split('[')[1].split(']')[0].split(',')[1].replace('\'','')
                if 'Integrase' in Type or 'IntI1' in Type:
                    # set integrase locus
                    Intlocus = Locus
                elif 'SulI' in Type:
                    # set sul1 locus
                    SulIlocus = Locus
            except IndexError:
                pass
    # add the integron to the phylogram
    # Finalkey: record the last integron element
    Finalkey=''
    # merge the replicated integrons
    Repli_node = Type_node
    # for an integron in the forward direction
    if (SulIlocus != 0 and Intlocus != 0 and Intlocus <= SulIlocus) or (
                    Intlocus == 0 and SulIlocus >= Locus / 2) or \
            (SulIlocus == 0 and 0 < Intlocus <= Locus / 2) or (Intlocus == SulIlocus == 0):
        # load the gene cassettes one by one, j stores the locus of a gene cassette
        j = 0
        while j <= (len(str(line).split('\t')[cint].split(';')) - 1):
            # key stores the annotation of a gene cassette
            # key = [integron_element_name, type, locus, annotation]
            key = str(line).split('\t')[cint].split(';')[j]
            j += 1
            if key != '':
                try:
                    Type = key.split('[')[1].split(']')[0].split(',')[1].replace('\'', '')
                    Label = key.split('[')[1].split(']')[0].split(',')[-1].replace(' ', '').replace('\'', '')
                    # add annotation
                    annotation = key.split('[')[1].split(']')[0].split(',')[3:-1]
                    if annotation == []:
                        # for integrases, sul1 and attC, ORF with no annotation
                        Label = Type
                    if Label in ARG_shortname:
                        # shorten the ARG phenotype label
                        Label = ARG_shortname[Label]
                    # record the replicated integron structure
                    Repli_node += ':' + Label
                    # record the last integron element
                    Finalkey = key
                except IndexError:
                    pass
    else:
        # for an integron in the reverse direction
        # start load the last integron element
        j = len(str(line).split('\t')[cint].split(';')) - 1
        while j >= 0:
            key = str(line).split('\t')[cint].split(';')[j]
            j -= 1
            if key != '':
                try:
                    Type = key.split('[')[1].split(']')[0].split(',')[1].replace('\'', '')
                    Label = key.split('[')[1].split(']')[0].split(',')[-1].replace(' ', '').replace('\'', '')
                    # add annotation
                    annotation = key.split('[')[1].split(']')[0].split(',')[3:-1]
                    if annotation == []:
                        # for integrases, sul1 and attC, ORF with no annotation
                        Label = Type
                    if Label in ARG_shortname:
                        # shorten the ARG phenotype label
                        Label = ARG_shortname[Label]
                    # record the replicated integron structure
                    Repli_node += ':' + Label
                    # record the last integron element
                    Finalkey = key
                except IndexError:
                    pass
    # count replicated integrons number
    if Repli_node in Repli:
        # for replicated integrons, simply add up the number of occurrence
        Repli[Repli_node][0]+=1
    else:
        # connect the first integron element to the phylogram tree
        if (SulIlocus!=0 and Intlocus!=0 and Intlocus<=SulIlocus) or (Intlocus==0 and SulIlocus>=Locus/2) or\
            (SulIlocus==0 and 0<Intlocus<=Locus/2) or (Intlocus==SulIlocus==0):
            # for forward integrons
            edge_add(Type_node, str(line).split('\t')[cint].split(';')[0], 'structure', 'None', 10)
        else:
            # for reverse integrons
            edge_add(Type_node, str(line).split('\t')[cint].split(';')[-2], 'structure', 'None', 10)
        # contruct the integron structure
        # Lastkey: record a previous added integron element
        Lastkey = ''
        for key in str(line).split('\t')[cint].split(';'):
            if key != '':
                try:
                    integron_element_name = key.split('[')[1].split(']')[0].split(',')[0].replace('\'', '')
                    Type = key.split('[')[1].split(']')[0].split(',')[1].replace('\'', '')
                    Label = key.split('[')[1].split(']')[0].split(',')[-1].replace(' ', '').replace('\'', '')
                    # add annotation
                    annotation = key.split('[')[1].split(']')[0].split(',')[3:-1]
                    if annotation != []:
                        annotation.append(key.split('[')[1].split(']')[0].split(',')[-1])
                    if 'ORF' in Type:
                        if len(annotation) == 2:
                            # metal resistance
                            Type = 'MRG'
                        elif len(annotation) == 3:
                            # antibiotic resistance
                            Type = 'ARG'
                        else:
                            # for ORF of no annotation
                            Label = Type
                        # shorten the ARG phenotype label
                        if Label in ARG_shortname:
                            Label = ARG_shortname[Label]
                            Color = Label
                        else:
                            Label = ''
                            Color = Type
                    else:
                        # for integrases, sul1 and attC
                        Label = ''
                        Color=Type
                    # add the integron node
                    if Type.replace(' ','') == 'attC':
                        node_add(key, Type, Label, 'black', Color, 50)
                    else:
                        node_add(key, Type, Label, 'black', Color, 100)
                    if Lastkey != '':
                        # add edge for integron structure
                        if (SulIlocus != 0 and Intlocus != 0 and Intlocus <= SulIlocus) or (
                                        Intlocus == 0 and SulIlocus >= Locus / 2) or \
                                (SulIlocus == 0 and 0 < Intlocus <= Locus / 2) or (Intlocus == SulIlocus == 0):
                            # for forward integron
                            edge_add(Lastkey, key, 'structure', 'None', 10)
                        else:
                            # for reverse integron
                            edge_add(key, Lastkey, 'structure', 'None', 10)
                    # update the previous integron element
                    Lastkey = key
                except IndexError:
                    pass
        # connect the occurrence number to the last integron element of the replicated integron
        Repli.setdefault(Repli_node, [1, Finalkey])


def netword_add(Repli):
    # output the number of occurrence for replicated integrons
    for Repli_node in Repli:
        if Repli[Repli_node][0] > 1:
            node_add(Repli_node, 'Occurrence', '*'+str(Repli[Repli_node][0]), 'yellow', 'Occurrence', 150)
            edge_add(Repli[Repli_node][1], Repli_node, 'Occurrence', 'Occurrence', 0)


################################################### Programme #######################################################
# load ORF length information
#ORFs_name=dict()
#for line in open(args.ot,'r'):
#    ORFs_name.setdefault(str(line).split('\t')[0],
#                         abs(float(str(line).split('\t')[-3]) - float(str(line).split('\t')[-2])))
# shorten the ARG phenotype name for display
ARG_shortname=dict()
for line in open('database/ARG_shortname.txt','r'):
    ARG_shortname.setdefault(str(line).split('\t')[0],str(line).split('\t')[1].split('\r')[0].split('\n')[0])
# set column number for the integron gene cassettes and taxonomy
if args.tx != 'None.norm':
    # taxonomy column
    c1 = int(str(args.tc).split(',')[0].replace(' ', '')) + 5
    c2 = int(str(args.tc).split(',')[1].replace(' ', '')) + 5
# integron gene cassettes column
cint = 3
# draw integorn phylogram
for anno_file in list_file:
    # initialize node and edge tables
    nodetable = []
    edgetable = []
    # write labels in the edge and node tables
    node_add('Node', 'Shape', 'Label', 'Labelcolor', 'Color', 'Length')
    if args.tx != 'None.norm':
        node_add('Bacteria', 'taxon', 'Bacteria', 'black', 'None', 50)
    else:
        node_add('Phylogram', 'taxon', 'Phylogram', 'black', 'None', 50)
    edge_add('Node1', 'Node2', 'Type', 'Color', 'Length')
    # merge the replicated integrons
    Repli = dict()
    for line in open(os.path.join(annotationdir,anno_file),'r'):
        if 'Type' not in str(line).split('\t')[0]:
            # main function of contructing phylogram
            network(str(line),nodetable,edgetable,Repli,args.tx)
    # output the number of occurrence for replicated integrons
    netword_add(Repli)