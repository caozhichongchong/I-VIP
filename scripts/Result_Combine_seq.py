import os
import glob
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-input",
                    help="input path", type=str, default='.',metavar='.')
parser.add_argument("--fasta",
                    help="format of fasta sequences", type=str, default='.fa', metavar='.fa')
parser.add_argument('--resultdir',
                    default="result", action='store', type=str, metavar='result',
                    help="Set the result directory for integron identification (default: result)")
parser.add_argument('--searchdir',
                    default="output", action='store', type=str, metavar='output',
                    help="Set the directory of attC and SulI searching results (default: output)")
parser.add_argument('--distance',
                    default=4000, action='store', type=int, metavar='4000',
                    help='Distance cutoff for two cassettes to be clustered together (default is 4kb)')
parser.add_argument('--attc_evalue',
                    default=1.0, action='store', type=float, metavar='1.0',
                    help='E-value cutoff for attC site identification (default is 1.0)')
parser.add_argument('--AA_type',
                    help="The AA type, eg: 0 for genbank parsing; 1 for prodigal prediction",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=0, type=int)
parser.add_argument("--orf",
                    help="format of ORFs", type=str, default='.faa', metavar='.faa')
################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.fasta)
distance = int(args.distance)
search_path = os.path.abspath(args.searchdir)
in_dir, input_file = os.path.split(args.input)
in_dir = os.path.abspath(in_dir)
list_fasta = glob.glob(os.path.join(in_dir,'*'+args.fasta))
try:
    os.mkdir(args.resultdir)
except OSError:
    pass
Cutoff_attc = float(args.attc_evalue)
Cutoff_identity=80
Cutoff_hitlength=50

#IntI and SulI length
Length = dict()
for line in open(os.path.join(in_dir, 'database/' + 'Integrase.fasta.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))
for line in open(os.path.join(in_dir, 'database/' + 'sul1_sequences_in_SARG.txt.length.txt'), 'rb'):
    Length.setdefault(str(line).split('\t')[0], float(str(line).split('\t')[2]))

#IntI1 list
IntI1=[]
for line in open(os.path.join(in_dir, 'database/' + 'IntI1_list.txt'), 'rb'):
    IntI1.append(str(line).split('\r')[0].split('\n')[0])

################################################### Function #########################################################
def int_blast_list(file, list,Cutoff_identity,Cutoff_hitlength):
    try:
        for line in open(file, 'rb'):
            if float(str(line).split('\t')[2])>=Cutoff_identity:
                if float(str(line).split('\t')[3])>=Cutoff_hitlength*float(Length.get(str(line).split('\t')[1]))/100:
                    Contig = str(line).split('\t')[0]
                    if list.get(Contig, 'None') == 'None':
                        list.setdefault(Contig,
                                        '{0} {1} {2}'.format(str(line).split('\t')[2], str(line).split('\t')[3],
                                                             str(line).split('\t')[1]))
                    elif float(str(line).split('\t')[2]) > float(str(list[Contig]).split(' ')[0]):
                        list[Contig] = '{0} {1} {2}'.format(str(line).split('\t')[2], str(line).split('\t')[3],
                                                                              str(line).split('\t')[1])
    except IOError:
        pass

def int_attc_list(file, list):
    for line in open(file, 'rb'):
        if float(str(line).split('\t')[15])<=Cutoff_attc:
            if list.get(str(line).split('\t')[0], 'None') == 'None':
                list.setdefault(str(line).split('\t')[0],
                                [[str(line).split('\t')[15], str(line).split('\t')[7],
                                  str(line).split('\t')[8]]])
            elif all(str(line).split('\t')[7] not in key and str(line).split('\t')[8] not in key for key in
                         list[str(line).split('\t')[0]]):
                list[str(line).split('\t')[0]].append(
                    [str(line).split('\t')[15], str(line).split('\t')[7],
                     str(line).split('\t')[8]])


def ORFannotation(otherORFs):
    Annotation = []
    if SulI_list.get(otherORFs, '-') != '-':
        if not Annotation:
            Annotation.append([(float(ORFs_name.get(otherORFs, 'None')[0]) + float(
                ORFs_name.get(otherORFs, 'None')[1])) / 2])
            Annotation.append([str(otherORFs) + ' ' + 'SulI' + ' ' + SulI_list.get(otherORFs, '-')])
            Annotation.append(['SulI'])
        else:
            Annotation[1].append(str(otherORFs) + ' ' + 'SulI' + ' ' + SulI_list.get(otherORFs, '-'))
            Annotation[2].append('SulI')
    if Integrase_list.get(otherORFs, '-') != '-':
        Tag='Integrase'
        if Integrase_list.get(otherORFs, '-').split(' ')[2] in IntI1:
            Tag='IntI1'
        if not Annotation:
            Annotation.append([(float(ORFs_name.get(otherORFs, 'None')[0]) + float(
                ORFs_name.get(otherORFs, 'None')[1])) / 2])
            Annotation.append([str(otherORFs) + ' ' + Tag + ' ' + Integrase_list.get(otherORFs, '-')])
            Annotation.append([Tag])
        else:
            Annotation[1].append(str(otherORFs) + ' ' + Tag + ' ' + Integrase_list.get(otherORFs, '-'))
            Annotation[2].append(Tag)

    return Annotation



def FindORF(combine, key, contig, site1, site2, cutoff=4000):
    Contigs_name.get(contig, 'None')
    if Contigs_name.get(contig, 'None')!='None':
        for otherORFs in Contigs_name.get(contig, 'None'):
            if ORFs_name[otherORFs] != 'None':
                Annotation = ORFannotation(otherORFs)
                if Annotation:
                    sitemid = float(Annotation[0][0])
                    if str(sitemid) not in combine[key][0]:
                        if (site1 - cutoff) <= sitemid <= (site2 + cutoff):
                            combine[key][2] = str(combine[key][2]) + ' ' + str(Annotation[2])
                            if sitemid > site2:
                                combine[key][0].append(str(sitemid))
                                combine[key][1].append(Annotation[1])
                                site2 = sitemid
                            elif sitemid < site1:
                                combine[key][0] = [str(sitemid)] + combine[key][0]
                                combine[key][1] = [Annotation[1]] + combine[key][1]
                                site1 = sitemid
                            else:
                                location = 0
                                if float(combine[key][0][0]) <= float(combine[key][0][-1]):
                                    for element in combine[key][0]:
                                        if float(element) < float(sitemid):
                                            location += 1
                                        else:
                                            break
                                    combine[key][0].insert(location, str(sitemid))
                                    combine[key][1].insert(location, Annotation[1])
                                else:
                                    for element in combine[key][0]:
                                        if float(element) > float(sitemid):
                                            location += 1
                                        else:
                                            break
                                    combine[key][0].insert(location, str(sitemid))
                                    combine[key][1].insert(location, Annotation[1])
                            ORFs_name[otherORFs] = 'None'
                            Contigs_name.get(contig, 'None').remove(otherORFs)
    return [site1, site2]


def Findattc(combine, key, contig, site1, site2, cutoff=4000):
    for attC in Attc_list.get(contig, '-'):
        if attC[-1] != '#':
            try:
                sitemid = (float(attC[1]) + float(
                    attC[2])) / 2
                if str(sitemid) not in combine[key][0]:
                    if (site1 - cutoff) <= sitemid <= (site2 + cutoff):
                        combine[key][2] = str(combine[key][2]) + ' ' + 'attC'
                        if sitemid > site2:
                            combine[key][0].append(str(sitemid))
                            combine[key][1].append('attC' + ' ' + attC[0])
                            site2 = max(float(attC[1]), float(
                                attC[2]))
                        elif sitemid < site1:
                            combine[key][0] = [str(sitemid)] + combine[key][0]
                            combine[key][1] = ['attC' + ' ' + attC[0]] + combine[key][1]
                            site1 = min(float(attC[1]), float(
                                attC[2]))
                        else:
                            location = 0
                            if float(combine[key][0][0]) <= float(combine[key][0][-1]):
                                for element in combine[key][0]:
                                    if float(element) < float(sitemid):
                                        location += 1
                                    else:
                                        break
                                combine[key][0].insert(location, str(sitemid))
                                combine[key][1].insert(location, 'attC' + ' ' + attC[0])
                            else:
                                for element in combine[key][0]:
                                    if float(element) > float(sitemid):
                                        location += 1
                                    else:
                                        break
                                combine[key][0].insert(location, str(sitemid))
                                combine[key][1].insert(location, 'attC' + ' ' + attC[0])
                        attC.append('#')
            except IndexError:
                print contig
                print attC
    return [site1, site2]



def Combine(target, lable, cutoff=4000):
    combine = dict()
    if lable == 'Contig':
        temp3 = target
        i = 1
        site1 = min(float(Attc_list.get(temp3, '-')[0][1]),
                    float(Attc_list.get(temp3, '-')[0][2]))
        site2 = max(float(Attc_list.get(temp3, '-')[0][1]),
                    float(Attc_list.get(temp3, '-')[0][2]))
        combine.setdefault(str(temp3) + ':' + str(i), [[str((site1 + site2) / 2)], ['attC' + ' ' + \
                                                                                    Attc_list.get(temp3, '-')[0][0]],
                                                       'attC'])
        Attc_list.get(temp3, '-')[0].append('#')
        Newsite =Findattc(combine, str(temp3) + ':' + str(i), temp3, float(site1), float(site2), cutoff)
        #combine attC sites until no attC was in the range
        while (float(site1) + float(site2)) / 2 != (float(Newsite[0]) + \
                                                            float(Newsite[1])) / 2:
            site1 = float(Newsite[0])
            site2 = float(Newsite[1])
            Newsite =Findattc(combine, str(temp3) + ':' + str(i), temp3, float(site1), float(site2), cutoff)
        if bool(Integrase_list) or bool(SulI_list): #combine sulI and integrase
            Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]), float(Newsite[1]),
                              cutoff)
        if str(Contigs_name.get(temp3, 'None'))!='None':
            Contigs_name.get(temp3, 'None').reverse()
            # reverse contig ORFs to add integrase and SulI to the other end
            if bool(Integrase_list) or bool(SulI_list):
                Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]), float(Newsite[1]),
                                  cutoff)
        # combine attC sites until no attC was in the range
        while (float(site1) + float(site2)) / 2 != (float(Newsite[0]) + \
                                                            float(Newsite[1])) / 2:
            site1 = float(Newsite[0])
            site2 = float(Newsite[1])
            Newsite = Findattc(combine, str(temp3) + ':' + str(i), temp3, float(site1), float(site2), cutoff)
            if bool(Integrase_list) or bool(SulI_list):
                Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]), float(Newsite[1]),
                                  cutoff)
            if str(Contigs_name.get(temp3, 'None')) != 'None':
                Contigs_name.get(temp3, 'None').reverse()
                if bool(Integrase_list) or bool(SulI_list):
                    Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]), float(Newsite[1]),
                                      cutoff)
        for otherattC in Attc_list.get(temp3, '-'): #for more than on integrons on a contig
            if otherattC[-1] != '#':
                i += 1
                site1 = min(float(otherattC[1]),
                            float(otherattC[2]))
                site2 = max(float(otherattC[1]),
                            float(otherattC[2]))
                combine.setdefault(str(temp3) + ':' + str(i), [[str((site1 + site2) / 2)], ['attC' + ' ' + \
                                                                                            Attc_list.get(temp3, '-')[
                                                                                                0][0]],
                                                               'attC'])
                Attc_list.get(temp3, '-')[0].append('#')
                Newsite =Findattc(combine, str(temp3) + ':' + str(i), temp3, float(site1), float(site2), cutoff)
                # combine attC sites until no attC was in the range
                while (float(site1) + float(site2)) / 2 != (float(Newsite[0]) + \
                                                                    float(Newsite[1])) / 2:
                    site1 = float(Newsite[0])
                    site2 = float(Newsite[1])
                    Newsite = Findattc(combine, str(temp3) + ':' + str(i), temp3, float(site1), float(site2), cutoff)
                if bool(Integrase_list) or bool(SulI_list):
                    Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]), float(Newsite[1]),
                                      cutoff)
                if str(Contigs_name.get(temp3, 'None')) != 'None':
                    Contigs_name.get(temp3, 'None').reverse()
                    # reverse contig ORFs to add integrase and SulI to the other end
                    if bool(Integrase_list) or bool(SulI_list):
                        Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]),
                                          float(Newsite[1]),
                                          cutoff)
                # combine attC sites until no attC was in the range
                while (float(site1) + float(site2)) / 2 != (float(Newsite[0]) + \
                                                                    float(Newsite[1])) / 2:
                    site1 = float(Newsite[0])
                    site2 = float(Newsite[1])
                    Newsite = Findattc(combine, str(temp3) + ':' + str(i), temp3, float(site1), float(site2), cutoff)
                    if bool(Integrase_list) or bool(SulI_list):
                        Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]),
                                          float(Newsite[1]),
                                          cutoff)
                    if str(Contigs_name.get(temp3, 'None')) != 'None':
                        Contigs_name.get(temp3, 'None').reverse()
                        if bool(Integrase_list) or bool(SulI_list):
                            Newsite = FindORF(combine, str(temp3) + ':' + str(i), temp3, float(Newsite[0]),
                                              float(Newsite[1]),
                                              cutoff)
    return combine


def write_integron(Name, Site, Annotation, Inttype):
    f1 = open(os.path.join(in_dir + '/' + args.resultdir, os.path.splitext(str(input_file))[0] + '.Integron.txt'), 'a')
    NumofattC = int(str(";".join([str(x) for x in Annotation])).count('attC'))
    if ('IntI1' in Inttype) and ('attC' in Inttype) \
            and ('SulI' in Inttype) and NumofattC >= 1:
        f1.write('A' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                 ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+'\n')
        f1.close()
    elif ('IntI1' in Inttype) and ('attC' in Inttype) \
            and ('SulI' not in Inttype) and NumofattC >= 1:
        f1.write('B' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                 ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
        f1.close()
    elif ('Integrase' in Inttype) and ('attC' in Inttype) \
            and NumofattC >= 1:
        f1.write('C' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                 ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
        f1.close()
    elif ('Integrase' not in Inttype) and ('IntI1' not in Inttype) and ('attC' in Inttype) \
            and ('SulI' in Inttype) and NumofattC >= 1:
        f1.write('D' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                 ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
        f1.close()
    elif ('Integrase' not in Inttype) and ('IntI1' not in Inttype) and ('attC' in Inttype) \
            and ('SulI' not in Inttype) and NumofattC >= 2:
        f1.write('E' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
                 ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
        f1.close()
    #elif ('Integrase' in Inttype) and ('attC' not in Inttype) and ('SulI' in Inttype):
    #    f1.write('A-' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
    #             ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
    #    f1.close()
    #elif ('Integrase' in Inttype) and ('attC' not in Inttype) and ('SulI' not in Inttype):
    #    f1.write('B-' + '\t' + str(Name) + '\t' + str(NumofattC) + '\t' +
    #             ";".join([str(x) for x in Site]) + '\t' + ";".join([str(x) for x in Annotation]) + '\t'+ str(Inttype)+ '\n')
    #    f1.close()


################################################### Programme #######################################################


ORFs_name=dict()
for line in open(os.path.join(in_dir, 'all.orfs.aa.length'),'rb'):
    if args.AA_type ==0:
        ORFs_name.setdefault(str(line).split('\t')[0],
                                     [float(str(line).split('\t')[1].split(' ')[-2]),
                                      float(str(line).split('\t')[1].split(' ')[-1])])
    elif args.AA_type ==1:
        ORFs_name.setdefault(str(line).split('\t')[0],
                                     [float(str(line).split('\t')[1].split(' # ')[1]),
                                      float(str(line).split('\t')[1].split(' # ')[2])])

for file_name in list_fasta:
    print 'Searching ' + str(file_name)
    in_dir, input_file = os.path.split(file_name)
    try:
        Attc_name = os.path.join(search_path, str(input_file) + '.Z.max.attc.hits.txt2.txt')
        Attc_list = dict()
        int_attc_list(Attc_name, Attc_list)
        Contigs_name = dict()
        Contigslist=[]
        for AA in ORFs_name:
            if input_file in str(AA):
                Contigname=str(AA)
                if args.orf + '_' in Contigname:
                    Contigname = Contigname.split(args.orf + '_')[1]
                if args.fasta + '_' in Contigname:
                    Contigname = Contigname.split(args.fasta + '_')[1]
                temp=Contigname.split('_')
                if len(temp)<2:
                    Contigname = Contigname.split('_')[0]
                else:
                    Contigname=Contigname.split('_')[0]+'_'+"_".join(Contigname.split('_')[1:-1])
                Contigslist.append(Contigname)
                if Contigs_name.get(Contigname, 'None') == 'None':
                    Contigs_name.setdefault(Contigname, [AA])
                else:
                    Contigs_name[Contigname].append(AA)
        Contigslist=list(Contigslist)
        try:
            SulI_name = os.path.join(search_path, str(input_file) + '.sul1.txt')
            SulI_list = dict()
            int_blast_list(SulI_name, SulI_list, 30, 50)
        except IOError:
            print 'No sulI for ' + str(input_file)
        try:
            IntI_name = os.path.join(search_path, str(input_file) + '.Int.txt')
            Integrase_list = dict()
            int_blast_list(IntI_name, Integrase_list, Cutoff_identity, Cutoff_hitlength)
        except IOError:
            print 'No Int for ' + str(input_file)
        Lable = 'Integron_Type' + '\t' + 'Contig_ID' + '\t' + 'Gene_Cassatte_Number' + '\t' + \
                'Gene_Locus' + '\t' + 'Gene_Annotation' + '\n'
        f1 = open(
            os.path.join(args.resultdir, os.path.splitext(str(input_file))[0] + '.Integron.txt'),
            'a')
        f1.write(Lable)
        f1.close()
        for Contigname in Contigslist:
            if bool(Attc_list) and Attc_list.get(Contigname, '-') != '-':
                result = Combine(Contigname, 'Contig', distance)
                for key in result:
                    write_integron(key, result[key][0], result[key][1], result[key][2])
    except IOError:
        print 'No attC for ' + str(input_file)
    print 'Cleaning Up...\n'
    del file_name
print 'Finished combining Integrons!'