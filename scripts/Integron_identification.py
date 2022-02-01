import os
import glob
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input', metavar='input')
parser.add_argument("-f",
                    help="file type or filename extension of your sequences",
                    type=str, default='.fa',metavar='.fa or .fasta or .fna')
parser.add_argument("-o",
                    help="Optional: to provide your own CDS or ORFs files, \
                                        please input the file type or filename extension of your CDS or ORFs",
                    type=str, default='.faa', metavar='.faa')
parser.add_argument("--r",
                    help="output directory or folder of your integron searching results",
                    type=str, default='I-VIP_result',metavar='I-VIP_result')
parser.add_argument('--d',
                    default=4000, action='store', type=int, metavar='4000',
                    help='Optional: set the distance cutoff for two cassettes to be clustered together\
                                         (default is 4000)')
parser.add_argument('--c',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='Optional: set the e-value cutoff for attC search (default is 1.0)')
parser.add_argument('--ot',
                    help="Optional: to provide your own CDS or ORFs files, \
                                         please input the orf type or the method you used to extract orfs,\
                                         eg: 1 for genbank parsing; 2 for prodigal prediction",
                    metavar="a file containing target files and their specific ORF format",
                    action='store', default='Temp/all.orf.length', type=str)
parser.add_argument('--m',
                    help="Optional: set the search strategies for Module A \
                    (1: global search by Module A1; 2: local search by Module A2), \
                    (default \'1\' for global search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=1, type=int)


################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.i)
distance = int(args.d)
searchdir = args.r + "/output"
search_path = os.path.abspath(searchdir)
list_fasta = glob.glob(os.path.join(input_path,'*'+args.f))
resultdir = args.r + "/result"
try:
    os.mkdir(resultdir)
except OSError:
    pass
Cutoff_attc = float(args.c)
# load IntI1 list
IntI1=[]
for line in open(os.path.join('database/' + 'IntI1_list.txt'), 'r'):
    IntI1.append(str(line).split('\r')[0].split('\n')[0])
# record Integron_identification process
flog = open('Integron_identification.log', 'w')


################################################### Function #########################################################


def int_blast_list(file, list):
    # load integrase and sul1 search results
    # into the ORFs
    try:
        for line in open(file, 'r'):
            # ORFs on the target sequence (Contig)
            AA = str(line).split('\t')[0]
            if list.get(AA, 'None') == 'None':
                # store the (identity hit-length subject-ID) of the search results
                list.setdefault(AA,
                                '{0} {1} {2}'.format(str(line).split('\t')[2], str(line).split('\t')[3],
                                                     str(line).split('\t')[1]))
            elif float(str(line).split('\t')[2]) > float(str(list[AA]).split(' ')[0]):
                # for duplicate hit on the same ORF, reserve those with the highest identity
                list[AA] = '{0} {1} {2}'.format(str(line).split('\t')[2], str(line).split('\t')[3],
                                                str(line).split('\t')[1])
    except IOError:
        pass


def formatattc(Attc_list, filename, recordid):
    # format Attc_list into {contig (filename + '_' + record.id) : [attC]}
    Templist = dict()
    for key in Attc_list:
        if filename in key and recordid in key:
            if filename + '_' + recordid not in Templist:
                Templist.setdefault(filename + '_' + recordid,[])
            # avoid duplicate attC
            for attc in Attc_list[key]:
                if attc not in Templist[filename + '_' + recordid]:
                    Templist[filename + '_' + recordid].append(attc)
    return Templist


def int_attc_list(file, list):
    # load attC search results
    # into the Contig (target sequence)
    for line in open(file, 'r'):
        if float(str(line).split('\t')[15]) <= Cutoff_attc:
            # filter out attC search results by e-value (Cutoff_attc)
            if list.get(input_file+'_'+str(line).split('\t')[0], 'None') == 'None':
                # load the first attC on one sequence
                # store the (e-value locus1 locus2) of the search results
                list.setdefault(input_file+'_'+str(line).split('\t')[0],
                                [[str(line).split('\t')[15], str(line).split('\t')[7],
                                  str(line).split('\t')[8]]])
            else:
                # load other attC on the same sequence
                # avoid duplicate attC
                if [str(line).split('\t')[15], str(line).split('\t')[7],
                     str(line).split('\t')[8]] not in list[input_file+'_'+str(line).split('\t')[0]]:
                    list[input_file+'_'+str(line).split('\t')[0]].append(
                        [str(line).split('\t')[15], str(line).split('\t')[7],
                         str(line).split('\t')[8]])


def ORFannotation(otherORFs):
    # ORFs_name = {ORF, [locus1, locus2]}
    # return annotation information
    # Annotation = [averge(locus1, locus2), ORF + ' ' + element_type + ' ' + search_information, element_type]
    # search_information = (identity hit-length subject-ID)
    Annotation = []
    # whether ORFs are annotated as sul1
    if SulI_list.get(otherORFs, '-') != '-':
        if not Annotation:
            # add average loci to Annotation
            Annotation.append((float(ORFs_name.get(otherORFs, 'None')[0]) + float(
                ORFs_name.get(otherORFs, 'None')[1])) / 2)
            # add annotation information to Annotation
            Annotation.append(str(otherORFs) + ' ' + 'SulI' + ' ' + SulI_list.get(otherORFs, '-'))
            # add element type to Annotation
            Annotation.append('SulI')
    # whether ORFs are annotated as integrases and intI1
    if Integrase_list.get(otherORFs, '-') != '-':
        Tag='Integrase'
        if Integrase_list.get(otherORFs, '-').split(' ')[2] in IntI1:
            Tag='IntI1'
        if not Annotation:
            # add average loci to Annotation
            Annotation.append((float(ORFs_name.get(otherORFs, 'None')[0]) + float(
                ORFs_name.get(otherORFs, 'None')[1])) / 2)
            # add annotation information to Annotation
            Annotation.append(str(otherORFs) + ' ' + Tag + ' ' + Integrase_list.get(otherORFs, '-'))
            # add element type to Annotation
            Annotation.append(Tag)
    return Annotation


def FindORF(combine, key, contig, site1, site2, cutoff = 4000):
    # add attC on one integron
    # format combine list into {Integron_name (key), [Site, Annotation, Inttype]}
    # Contigs_name = {contig, [ORFs]}
    # ORFs_name = {ORF, [locus1, locus2]}
    if Contigs_name.get(contig, 'None')!='None':
        for otherORFs in Contigs_name.get(contig, 'None'):
            if ORFs_name[otherORFs] != 'None':
                # a valid ORF
                # get ORF annotation
                Annotation = ORFannotation(otherORFs)
                if Annotation:
                    # for sul1 or integrases
                    sitemid = float(Annotation[0])
                    if str(sitemid) not in combine[key][0]:
                        # add a new integron element
                        if (site1 - cutoff) <= sitemid <= (site2 + cutoff):
                            # within the distance cutoff > merge into one integron
                            if str(Annotation[2]) not in combine[key][2]:
                                combine[key][2] = str(combine[key][2]) + ' ' + str(Annotation[2])
                                # add element type information
                            # compare site of integron elements
                            if sitemid > site2:
                                # add afterwards
                                combine[key][0].append(str(sitemid))
                                # add annotation informaiton of ORF and search result
                                combine[key][1].append(Annotation[1])
                                # change the integron site2
                                site2 = sitemid
                            elif sitemid < site1:
                                # add forwards
                                combine[key][0] = [str(sitemid)] + combine[key][0]
                                # add annotation informaiton of ORF and search result
                                combine[key][1] = [Annotation[1]] + combine[key][1]
                                # change the integron site1
                                site1 = sitemid
                            else:
                                # insert in the middle
                                # do not change the integron sites
                                location = 0
                                if float(combine[key][0][0]) <= float(combine[key][0][-1]):
                                    # loci order from the smallest to the largest
                                    for element in combine[key][0]:
                                        if float(element) < float(sitemid):
                                            location += 1
                                        else:
                                            break
                                    combine[key][0].insert(location, str(sitemid))
                                    combine[key][1].insert(location, Annotation[1])
                                else:
                                    # loci order from the largest to the smallest
                                    for element in combine[key][0]:
                                        if float(element) > float(sitemid):
                                            location += 1
                                        else:
                                            break
                                    combine[key][0].insert(location, str(sitemid))
                                    combine[key][1].insert(location, Annotation[1])
                            # delete used ORF
                            ORFs_name[otherORFs] = 'None'
                            Contigs_name.get(contig, 'None').remove(otherORFs)
    # return new integron loci
    return [site1, site2]


def Findattc(combine, key, contig, site1, site2, cutoff=4000):
    # add attC on one integron
    # format combine list into {Integron_name (key), [Site, Annotation, Inttype]}
    for attC in TempattC.get(contig, '-'):
        if attC[-1] != '#':
            # add valid attC
            try:
                sitemid = (float(attC[1]) + float(
                    attC[2])) / 2
                if str(sitemid) not in combine[key][0]:
                    # add a new integron element
                    if (site1 - cutoff) <= sitemid <= (site2 + cutoff):
                        # within the distance cutoff > merge into one integron
                        if 'attC' not in combine[key][2]:
                            combine[key][2] = str(combine[key][2]) + ' ' + 'attC'
                            # add element type (attC) information
                        # compare site of integron elements
                        if sitemid > site2:
                            # add afterwards
                            combine[key][0].append(str(sitemid))
                            # add annotation informaiton of attC and e-value
                            combine[key][1].append('attC' + ' ' + attC[0])
                            # change the integron site2
                            site2 = max(float(attC[1]), float(
                                attC[2]))
                        elif sitemid < site1:
                            # add forwards
                            combine[key][0] = [str(sitemid)] + combine[key][0]
                            # add informaiton of attC and e-value
                            combine[key][1] = ['attC' + ' ' + attC[0]] + combine[key][1]
                            # change the integron site1
                            site1 = min(float(attC[1]), float(
                                attC[2]))
                        else:
                            # insert in the middle
                            # do not change the integron sites
                            location = 0
                            if float(combine[key][0][0]) <= float(combine[key][0][-1]):
                                # loci order from the smallest to the largest
                                for element in combine[key][0]:
                                    if float(element) < float(sitemid):
                                        location += 1
                                    else:
                                        break
                                combine[key][0].insert(location, str(sitemid))
                                combine[key][1].insert(location, 'attC' + ' ' + attC[0])
                            else:
                                # loci order from the largest to the smallest
                                for element in combine[key][0]:
                                    if float(element) > float(sitemid):
                                        location += 1
                                    else:
                                        break
                                combine[key][0].insert(location, str(sitemid))
                                combine[key][1].insert(location, 'attC' + ' ' + attC[0])
                        # delete used attC
                        attC.append('#')
            except IndexError:
                print('Not valid attC format for: ' + str([contig, attC]))
                flog.write('Not valid attC format for: ' + str([contig, attC])+'\n')
    #return new integron loci
    return [site1, site2]


def integron_identify(combine, target, i, attC, cutoff):
    # main function for integron identification
    # start with the first attC, set loci
    site1 = min(float(attC[1]),
                float(attC[2]))
    site2 = max(float(attC[1]),
                float(attC[2]))
    # set up the first integron on target sequence (with the name of str(target) + ':' + str(i))
    # information contains [loci, integron element information, integron element type]
    combine.setdefault(str(target) + ':' + str(i),
                       [[str((site1 + site2) / 2)],
                        ['attC' + ' ' + attC[0]], 'attC'])
    # delete attC that is used
    attC.append('#')
    # add the first attC to the candidate integron and return the integron loci (Newsite)
    Newsite = Findattc(combine, str(target) + ':' + str(i), target, float(site1), float(site2), cutoff)
    # add attC sites until no attC was in the range of Newsite
    while (float(site1) + float(site2)) / 2 != (float(Newsite[0]) + float(Newsite[1])) / 2:
        site1 = float(Newsite[0])
        site2 = float(Newsite[1])
        Newsite = Findattc(combine, str(target) + ':' + str(i), target, float(site1), float(site2), cutoff)
    if bool(Integrase_list) or bool(SulI_list):
        # add sulI and integrase
        Newsite = FindORF(combine, str(target) + ':' + str(i), target, float(Newsite[0]), float(Newsite[1]),
                          cutoff)
    # Contigs_name store a list of ORF names for the target sequence
    if str(Contigs_name.get(target, 'None')) != 'None':
        Contigs_name.get(target, 'None').reverse()
        # reverse contig ORFs to add integrase and sulI to the other end
        if bool(Integrase_list) or bool(SulI_list):
            Newsite = FindORF(combine, str(target) + ':' + str(i), target, float(Newsite[0]), float(Newsite[1]),
                              cutoff)
    # after adding integrase and sul1 (integron range changes)
    # repeat the previous process to exhaustedly add integron elements
    while (float(site1) + float(site2)) / 2 != (float(Newsite[0]) +float(Newsite[1])) / 2:
        site1 = float(Newsite[0])
        site2 = float(Newsite[1])
        Newsite = Findattc(combine, str(target) + ':' + str(i), target, float(site1), float(site2), cutoff)
        # add integrase and sul1 genes
        if bool(Integrase_list) or bool(SulI_list):
            Newsite = FindORF(combine, str(target) + ':' + str(i), target, float(Newsite[0]), float(Newsite[1]),
                              cutoff)
        # Contigs_name store a list of ORF names for the target sequence
        if str(Contigs_name.get(target, 'None')) != 'None':
            # reverse the ORF list
            Contigs_name.get(target, 'None').reverse()
            if bool(Integrase_list) or bool(SulI_list):
                Newsite = FindORF(combine, str(target) + ':' + str(i), target, float(Newsite[0]), float(Newsite[1]),
                                  cutoff)


def Combine(target, lable, cutoff=4000):
    # main function for identify integrons on one target sequence
    # store integron structure in combine list
    combine = dict()
    # combine = {Integron_name, [Site, Annotation, Inttype]}
    # Integron_name = Contig + ':' + i
    if lable == 'Contig':
        i = 1
        integron_identify(combine, target, i, TempattC.get(target, '-')[0], cutoff)
        # the identificaiton of first integron on target sequence is finished
        # to identify the other integrons on the same target sequence
        # i+=1 and repeat the previous process
        for otherattC in TempattC.get(target, '-'):
            if otherattC[-1] != '#':
                # attC not added on identified integrons
                i += 1
                integron_identify(combine, target, i, otherattC, cutoff)
    return combine


def write_integron(Name, Site, Annotation, Inttype, No_integron):
    # output integron into a format of
    # 'Integron_Type' + '\t' + 'Contig_ID' + '\t' + 'Gene_Cassatte_Number' + '\t' +
    #            'Gene_Locus' + '\t' + 'Gene_Annotation' + '\n'
    f1 = open(os.path.join(resultdir, str(input_file).split(str(args.f))[0] + '.Integron.txt'), 'a')
    # count the number of attC on the integron
    NumofattC = int(str(";".join([str(x) for x in Annotation])).count('attC'))
    if NumofattC >= 2:
        # integron classification
        # Type A
        if ('IntI1' in Inttype) and ('attC' in Inttype) \
                and ('SulI' in Inttype):
            f1.write('A' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                     ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+'\n')
            f1.close()
        # Type B
        elif ('IntI1' in Inttype) and ('attC' in Inttype) \
                and ('SulI' not in Inttype):
            f1.write('B' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                     ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
            f1.close()
        # Type C
        elif ('Integrase' in Inttype) and ('attC' in Inttype):
            f1.write('C' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                     ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
            f1.close()
        # Type D
        elif ('Integrase' not in Inttype) and ('IntI1' not in Inttype) and ('attC' in Inttype) \
                and ('SulI' in Inttype):
            f1.write('D' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                     ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
            f1.close()
        # Type E
        elif ('Integrase' not in Inttype) and ('IntI1' not in Inttype) and ('attC' in Inttype) \
                and ('SulI' not in Inttype):
            f1.write('E' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                     ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
            f1.close()
    return No_integron


################################################### Programme #######################################################
# load ORF information of loci
# ORFs_name = {ORF, [locus1, locus2]}
ORFs_name=dict()
for line in open(args.ot,'r'):
    ORFs_name.setdefault(str(line).split('\t')[0],
                         [float(str(line).split('\t')[-3]),
                          float(str(line).split('\t')[-2])])
# identify integrons in each target sequence file
for file_name in list_fasta:
    print('Identifying integrons in ' + str(file_name))
    flog.write('Identifying integrons in ' + str(file_name)+'\n')
    input_path, input_file = os.path.split(file_name)
    No_integron = 0
    try:
        # load attC results
        # Attc_list = {Contig, [attc infor list]}
        Attc_name = os.path.join(search_path, str(input_file) + '.Z.max.attc.hits.txt2.txt')
        Attc_list = dict()
        int_attc_list(Attc_name, Attc_list)
        try:
            # load sul1 results
            # SulI_list = {ORF_name, sul1 infor}
            SulI_name = os.path.join(search_path, str(input_file) + '.sul1.txt')
            SulI_list = dict()
            int_blast_list(SulI_name, SulI_list)
        except IOError:
            print('No sul1 found for ' + str(input_file))
            flog.write('No sul1 found for ' + str(input_file)+'\n')
        try:
            # load integrase results
            # Integrase_list = {ORF_name, Integrase infor}
            IntI_name = os.path.join(search_path, str(input_file) + '.int.txt')
            Integrase_list = dict()
            int_blast_list(IntI_name, Integrase_list)
        except IOError:
            print('No integrase found for ' + str(input_file))
            flog.write('No integrase found ' + str(input_file)+'\n')
        # output integrons identified
        f1 = open(
            os.path.join(resultdir, os.path.splitext(str(input_file))[0] + '.Integron.txt'),
            'a')
        Lable = 'Integron_Type' + '\t' + 'Integron_Name' + '\t' + 'Number_Of_AttC' + '\t' + \
                'Integron_Element_Locus' + '\t' + 'Integron_Element_Type_Annotation' + '\n'
        f1.write(Lable)
        f1.close()
        # for each sequence in target files, search for integron
        for record in SeqIO.parse(open(os.path.join(input_path,input_file)), "fasta"):
            if bool(Attc_list):
                TempattC = formatattc(Attc_list, input_file, record.id)
                if TempattC.get(input_file + '_' + record.id, '-') != '-':
                    Contigs_name = dict()
                    Contig = input_file + '_' + record.id
                    # set up ORFs for each sequences (contigs)
                    # ORFs_name = {ORF, [locus1, locus2]}
                    for AA in ORFs_name:
                        if Contig in str(AA).replace(args.o,args.f):
                            # Contigs_name = {Contig, [ORF list]}
                            if Contigs_name.get(Contig, 'None') == 'None':
                                Contigs_name.setdefault(Contig, [AA])
                            else:
                                Contigs_name[Contig].append(AA)
                    # main step to identify integrons
                    result = Combine(Contig, 'Contig', distance)
                    for key in result:
                        # output each integron
                        # result = {Integron_name, [Site, Annotation, Inttype]}
                        No_integron = write_integron(key, result[key][0], result[key][1], result[key][2], No_integron)
        print ('Identified %s integrons in %s' % (No_integron, file_name))
    except IOError:
        print('No attC found in ' + str(input_file))
        flog.write('No attC found in ' + str(input_file)+'\n')
print('Finished identifying integrons!')
flog.write('Finished identifying integrons!\n')
flog.close()