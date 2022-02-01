import os
import glob
import argparse
from Bio import SeqIO

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument("--tx",
                    help="a file of taxonomy metadata (under your input folder)",
                    type=str, default='genbank_taxon.txt',
                    metavar='genbank_taxon.txt or None')
parser.add_argument("--tc",
                    help="column number corresponding to the taxonomy, i.e., from phylum to strain",
                    type=str, default='2,8', metavar='start column number, end column number')
################################################## Definition ########################################################


args = parser.parse_args()
anno_dir, anno_file = os.path.split(os.path.abspath(args.tx))
anno_dir = os.path.abspath(args.i)
try:
    c1=int(str(args.tc).split(',')[0].replace(' ',''))
    c2=int(str(args.tc).split(',')[1].replace(' ',''))
except ValueError:
    print('Wrong input for --tc \nPlease input in the format of: start column number, end column number\n"\
                      "Example: 2,8\nProceed with default value of 2,8\n')
################################################### Function #########################################################


def normalize(f1,f2):
    # correct "NA" taxonamy into "the taxonomy in a higher level" + "_" + "taxonomy_level"
    for line in f1:
        try:
            if line.split('\r')[0].split('\n')[0]!='':
                newline='\t'.join(str(line).split('\t')[0:(c1-1)])
                for i in range((c1-1),(c2)):
                    if str(line).split('\t')[i].split('\r')[0].split('\n')[0] not in ['NA','','None'] \
                           and 'environmental samples' not in str(line).split('\t')[i].split('\r')[0].split('\n')[0]:
                        newline+='\t'+str(line).split('\t')[i].split('\r')[0].split('\n')[0] + '_' + str(i)
                    else:
                        if c1==i: # replace 'NA', 'None' or empty annotation at the phylum (lowest taxonomy level)
                            if str(line).split('\t')[(c1-1)].split('\r')[0].split('\n')[0] in ['NA','','None'] \
                           or 'environmental samples' in str(line).split('\t')[(c1-1)].split('\r')[0].split('\n')[0]:
                                newline += '\tunclassified Bacteria_' + str(i)
                            else:
                                newline += '\t' + str(line).split('\t')[(c1 - 1)].split('\r')[0].split('\n')[0]\
                                           + '_' + str(i)
                        elif c1-1==i:
                            newline += '\tunclassified Bacteria_' + str(i)
                        else:
                            Out=False
                            for j in reversed(list(range((c1-1), i))):  # replace 'NA', 'None' or empty annotation
                                if str(line).split('\t')[j].split('\r')[0].split('\n')[0] not in ['NA','','None'] and \
                                        'environmental samples' not in str(line).split('\t')[j].split('\r')[0].split('\n')[0] :
                                    newline += '\t' + str(line).split('\t')[j].split('\r')[0].split('\n')[0] + '_' + str(i)
                                    Out=True
                                    break
                            if j==c1-1 and Out==False:
                                newline += '\tunclassified Bacteria_' + str(i)
                try:
                    newline += '\t' + '\t'.join(str(line).split('\t')[(c2):len(str(line).split('\t'))])
                    newline = newline.split('\r')[0].split('\n')[0] + '\n'
                except IndexError:
                    newline += '\n'
        except KeyError:
            print('The input c1 or c2 is out of the range of column number!\n' \
                  'Please input the right column numbers\n')
            break
        f2.write(newline)
    f1.close()
    f2.close()


def normalize2(f3):
    # correct taxonomy level for each taxonomy, by voting
    Taxon = dict()
    for lines in f3:
        for i in range((c1-1), (c2)):
            try:
                node2 = str(lines).split('\t')[i + 1]
                node1 = str(lines).split('\t')[i]
                if node2 not in Taxon:
                    Taxon.setdefault(node2, [[node1 + '--' + node2, 1]])
                else:
                    Set = 0
                    for key in Taxon[node2]:
                        if node1 + '--' + node2 == key[0]:
                            key[1] += 1
                            Set = 1
                            break
                    if Set == 0:
                        Taxon[node2].append([node1 + '--' + node2, 1])
            except IndexError:
                pass
    Taxon2 = dict()
    for node2 in Taxon:
        newkey = ''
        max = 0
        for key in Taxon[node2]:
            if key[1] > max:
                max = key[1]
                newkey = key
        Taxon2.setdefault(node2, newkey)
    f3.close()
    return Taxon2


def normalize3(f3,f4,Taxon2):
    # output normalized taxonomy information
    for line in f3:
        newline = '\t'.join(str(line).split('\t')[0:(c1-1)])
        node2 = str(line).split('\t')[(c2-1)].split('\r')[0].split('\n')[0]
        templine = [node2]
        for i in reversed(list(range((c1-1), (c2-1)))):
            node1 = Taxon2[node2][0].split('--')[0]
            templine.append(node1)
            node2 = node1
        for newnode in reversed(templine):
            newline+='\t'+newnode
        try:
            newline += '\t' + '\t'.join(str(line).split('\t')[(c2):len(str(line).split('\t'))])
            newline = newline.split('\r')[0].split('\n')[0] + '\n'
        except IndexError:
            newline += '\n'
        f4.write(newline)
    f3.close()
    f4.close()


################################################### Programme #######################################################


# main programme
# correct "NA" taxonamy into "the taxonomy in a higher level" + "_" + "taxonomy_level"
normalize(open(os.path.join(anno_dir, anno_file),'r'),
          open(os.path.join(anno_dir, anno_file+'.temp'),'w'))
# correct taxonomy level for each taxonomy, by voting
Taxon2=normalize2(open(os.path.join(anno_dir, anno_file + '.temp'), 'r'))
# output normalized taxonomy information
normalize3(open(os.path.join(anno_dir, anno_file + '.temp'), 'r'),
           open(os.path.join(anno_dir, anno_file)+ '.norm','w'),Taxon2)
