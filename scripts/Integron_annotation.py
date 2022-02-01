import os
import glob
import argparse


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument("-f",
                    help="file type or filename extension of your sequences.\n \
                         To input genbank file, set \"-f .gbff\" or \"-f .gbff.gz\"",
                    type=str, default='.fa',metavar='.fa or .fasta or .fna')
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result', metavar='I-VIP_result')
parser.add_argument('--a',
                    help="Optional: annotate gene cassettes against antibiotic and metal resistance databases",
                    metavar="Y or N", action='store', default='N', type=str)
parser.add_argument("--tx",
                    help="a file of taxonomy metadata (under your input folder)",
                    type=str, default='genbank_taxon.txt',
                    metavar='genbank_taxon.txt or None')
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running I-VIP (default 15)",
                    metavar="15", action='store', default=15, type=int)
parser.add_argument('--blastp',
                    help="Optional: complete path to blastp if not in PATH,",
                    metavar="/usr/local/bin/blastp",
                    action='store', default='blastp', type=str)
################################################## Definition ########################################################


args = parser.parse_args()
resultdir = args.r + "/Integron"
# optimal input of taxonomy information
if args.tx != 'None':
    tx_dir, tx_file = os.path.split(os.path.abspath(args.tx))
    tx_dir = os.path.abspath(args.i)
    taxon = dict()
    try:
        for line in open(os.path.join(tx_dir, tx_file),'r'):
            taxon.setdefault(str(line).split('\t')[0].split(args.f)[0],
                             line.split('\r')[0].split('\n')[0])
    except IOError:
        pass
tx_dir = os.path.abspath(args.i)


################################################### Function #########################################################


def filter_blast_list(file, Cutoff_identity, Cutoff_hitlength, GC_annotation):
    try:
        f = open(file + '.filter.txt', 'w')
        for line in open(file, 'r'):
            if float(str(line).split('\t')[2]) >= Cutoff_identity:
                if float(str(line).split('\t')[3]) >= Cutoff_hitlength * float(
                        Length.get(str(line).split('\t')[1])) / 100:
                    ORF = str(line).split('\t')[0]
                    Reference= str(line).split('\t')[1]
                    f.write(str(line))
                    if Reference in Structure:
                        Annotation = Reference + '\t' + Structure.get(Reference)
                    else:
                        Annotation = Reference
                    if ORF not in GC_annotation:
                        GC_annotation.setdefault(ORF, Annotation)
                    else:
                       pass
        f.close()
        return
    except IOError:
        print(str(file) + ' is missing!\n')
        flog.write(str(file) + ' is missing!\n')


def ORF_annotation(Intstruc, GC_MRG, GC_ARG, output):
    try:
        f = open(output, 'w')
        for line in open(Intstruc, 'r'):
            ORF = str(line).split('\t')[1]
            line = line.split('\r')[0].split('\n')[0]
            if ORF in GC_ARG:
                # add ARG annotation (reference, genotype, phenotype)
                line += '\t'+ GC_ARG[ORF]
            else:
                line += '\tNone\tNone\tNone'
            if ORF in GC_MRG:
                # add MRG annotation (reference, phenotype)
                line += '\t'+ GC_MRG[ORF] + '\t'+ GC_MRG[ORF].split(':')[0]
            else:
                line += '\tNone\tNone'
            f.write(line+'\n')
        f.close()
    except IOError:
        print(str(Intstruc) + ' is missing!\n')
        flog.write(str(Intstruc) + ' is missing!\n')


################################################### Programme #######################################################
# annotate gene cassettes again antibiotic and metal resistant databases
flog = open('Integron_annotation.log','w')
if args.a == 'Y':
    cmdARG = str(args.blastp) + " -query " + str(resultdir) + "/all.ORFs.fasta -db database/SARG.db.fasta -out " \
             + str(
        resultdir) + "/all.ORFs.fasta.ARG.blast.txt -outfmt 6 -evalue 1e-5 -num_threads " + \
             str(args.t) + " \n"
    cmdARG += str(args.blastp) + " -query " + str(resultdir) + "/all.ORFs.ClassI.fasta -db database/SARG.db.fasta -out " \
              + str(
        resultdir) + "/all.ORFs.ClassI.fasta.ARG.blast.txt -outfmt 6 -evalue 1e-5 -num_threads " + \
              str(args.t) + " \n"
    cmdMRG = str(args.blastp) + " -query " + str(resultdir) + "/all.ORFs.fasta -db database/MRG.db.fasta -out " \
             + str(
        resultdir) + "/all.ORFs.fasta.MRG.blast.txt -outfmt 6 -evalue 1e-5 -num_threads " + \
             str(args.t) + " \n"
    cmdMRG += str(args.blastp) + " -query " + str(resultdir) + "/all.ORFs.ClassI.fasta -db database/MRG.db.fasta -out " \
              + str(
        resultdir) + "/all.ORFs.ClassI.fasta.MRG.blast.txt -outfmt 6 -evalue 1e-5 -num_threads " + \
              str(args.t) + " \n"
    print(cmdARG)
    print(cmdMRG)
    os.system(cmdARG)
    os.system(cmdMRG)
    # load ARG and MRG reference length
    Length = dict()
    for line in open(os.path.join('database/' + 'MRG.db.fasta.length.txt'), 'r'):
        Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
    for line in open(os.path.join('database/' + 'SARG.db.fasta.length.txt'), 'r'):
        Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
    # load ARG structure (genotype and phenotype)
    Structure = dict()
    for line in open(os.path.join('database/' + 'SARG.structure.txt'), 'r'):
        Structure.setdefault(str(line).split('\t')[0], str(line).split('\t')[1] + '\t'
                             +str(line).split('\t')[2].split('\r')[0].split('\n')[0])
    # filter blast results and annotate ORF on integrons of other classes
    GC_ARG = dict()
    GC_MRG = dict()
    filter_blast_list(str(resultdir) + "/all.ORFs.fasta.MRG.blast.txt", 80, 90, GC_MRG)
    filter_blast_list(str(resultdir) + "/all.ORFs.fasta.ARG.blast.txt", 90, 80, GC_ARG)
    # output gene cassette annotation on integrons of other classes
    ORF_annotation(str(resultdir) + "/all.Integron_structure.txt", GC_MRG, GC_ARG,
                   str(resultdir) + "/all.Integron_structure.annotated.txt")
    # filter blast results and annotate ORF on class 1 integrons
    GC_ARG = dict()
    GC_MRG = dict()
    filter_blast_list(str(resultdir) + "/all.ORFs.ClassI.fasta.MRG.blast.txt", 80, 90, GC_MRG)
    filter_blast_list(str(resultdir) + "/all.ORFs.ClassI.fasta.ARG.blast.txt", 90, 80, GC_ARG)
    # output gene cassette annotation on class 1 integrons
    ORF_annotation(str(resultdir) + "/all.Integron_structure.ClassI.txt", GC_MRG, GC_ARG,
                   str(resultdir) + "/all.Integron_structure.ClassI.annotated.txt")
else:
    GC_ARG = dict()
    GC_MRG = dict()
    # output gene cassette annotation on class 1 integrons
    ORF_annotation(str(resultdir) + "/all.Integron_structure.ClassI.txt", GC_MRG, GC_ARG,
                   str(resultdir) + "/all.Integron_structure.ClassI.annotated.txt")
    GC_ARG = dict()
    GC_MRG = dict()
    # output gene cassette annotation on integrons of other classes
    ORF_annotation(str(resultdir) + "/all.Integron_structure.txt", GC_MRG, GC_ARG,
                   str(resultdir) + "/all.Integron_structure.annotated.txt")
# format integron structure
# load integron annotation files
list_file = [str(resultdir) + "/all.Integron_structure.ClassI.annotated.txt",
            str(resultdir) + "/all.Integron_structure.annotated.txt"]
for file_name in list_file:
    filedir, file_name = os.path.split(file_name)
    f1 = open(os.path.join(str(resultdir), str(file_name).replace('.annotated.txt',
                           '.annotated.formated.txt')), 'a')
    integron = 0
    lastintegron = ''
    ARG_MRG_number = 0
    Output = 0
    for line in open(os.path.join(str(resultdir),file_name),'r'):
        try:
            integronclass = str(line).split('\t')[0]
            targetfile = str(line).split('\t')[1].split(args.f)[0]
            integronname = str(line).split('\t')[1].split(':')[0]+':'+str(line).split('\t')[1].split(':')[1]
            ORFname=str(line).split('\t')[1]
            Tag = str(line).split('\t')[2]
            Locus = str(line).split('\t')[3]
            annotation = []
            if lastintegron != integronname:
                # output a new integron
                if integron != 0:
                    f1.write('\t' + str(ARG_MRG_number) + '\t')
                    ARG_MRG_number = 0
                    Output = 1
                    if args.tx != 'None':
                        # optimal output of taxonomy information
                        f1.write(str(taxon.get(targetfile, '')) + '\n')
                    else:
                        f1.write('\n')
                # output the first integron element of this integron
                integron += 1
                f1.write(str(integronclass) + '\t' + str(targetfile) + '\t' + str(integronname)+'\t')
            # continuously output the integron elements of this integron
            annotation.append(ORFname)
            annotation.append(Tag)
            annotation.append(Locus)
            for col in [4, 5, 6, 7, 8]:
                if str(line).split('\t')[int(col)].split('\n')[0].split('\r')[0] != 'None':
                    annotation.append(str(line).split('\t')[int(col)].split('\n')[0].split('\r')[0])
                    if col in [4 , 7] and Tag != 'SulI':
                        #avoid duplicate counting
                        ARG_MRG_number += 1
            f1.write(str(annotation) + ';')
            lastintegron = integronname
        except IndexError:
            pass
    if Output == 1:
    # for not empty integron result
    # output the last integron element
        f1.write('\t' + str(ARG_MRG_number) + '\t')
        if args.tx != 'None':
            # optimal output of taxonomy information
            f1.write(str(taxon.get(targetfile, '')) + '\n')
        else:
            f1.write('\n')
    f1.close()
flog.close()
