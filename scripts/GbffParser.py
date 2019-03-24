from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import glob
import argparse


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input directory or folder of your sequences",
                    type=str, default='input',metavar='input')
parser.add_argument("-f",
                    help="file type or filename extension of your sequences.\n \
                         To input genbank file, set \"-f .gbff\" or \"-f .gbff.gz\"",
                    type=str, default='.fa',metavar='.fa or .fasta or .fna')


################################################## Definition ########################################################
args = parser.parse_args()


################################################### Programme #######################################################
print "GbffParser.py was developed by Li-Guan Li and An-Ni Zhang"
print "To extract the genome sequences and CDS sequences from both complete and draft genomes"
print "Script is available to extract complete genome sequence and CDS sequences:"
print "https://github.com/LiguanLi/ARG_MRG_Cooccurrence/blob/master/GBK_CDSextract.py"
# set up gbff files
in_dir = os.path.abspath(args.i)
if args.f == '.gbff.gz':
    # gunzip .gbff.gz
    cmd1 = 'for file in ' + str(in_dir) + '/*.gbff.gz\ndo\ngunzip \"${file}\"\ndone'
    os.system(cmd1)
    input_extension = '.gbff'
else:
    input_extension = args.f
# fetch all gbff files
list_gbff = glob.glob(os.path.join(in_dir,'*'+input_extension))
# main programme
if list_gbff == []:
    print 'Fing no genbank files (.gbff or .gbff.gz) in ' + str(in_dir)
else:
    Filenumber = 0
    for file_name1 in list_gbff:
        gb = SeqIO.parse(os.path.join(in_dir, file_name1), "genbank")
        in_dir, file_name1 = os.path.split(file_name1)
        tmp1 = file_name1.split("_")[0]
        tmp2 = file_name1.split("_")[1]
        assem_name = "_".join([tmp1, tmp2])
        print "Processing GenBank chr record %s" % assem_name
        Filenumber += 1
        try:
            os.mkdir('Chr/Chr' + str(int(Filenumber / 10000)))
        except OSError:
            pass
        for record in gb:
            for feature in record.features:
                if feature.type == "source":
                    # record taxaid
                    taxaid = feature.qualifiers.get('db_xref', 'None')
            # output genome sequence
            output_handle1 = open(os.path.join(in_dir, assem_name + ".fa"), "a")
            output_handle1.write(
                ">%s_%s %s %s\n%s\n" % (assem_name, record.id, record.description, taxaid[-1], record.seq))
            # fetch CDS sequence
            AA_true = 0
            for feature in record.features:
                if feature.type == "CDS":
                    AA_true = 1
                    ID = feature.qualifiers.get('db_xref', 'None')[0]
                    desc = feature.qualifiers.get('protein_id', 'None')[0]
                    product = feature.qualifiers.get('product', 'None')[0]
                    locus = feature.qualifiers.get('locus_tag', 'None')[0]
                    type = feature.type
                    start = feature.location.start.position
                    end = feature.location.end.position
                    try:
                        assert len(feature.qualifiers['translation']) == 1
                        aa_seq = feature.qualifiers['translation'][0]
                    except KeyError:
                        print desc, 'no amni acids found!'
                        aa_seq = ''
                    if AA_true == 1:
                        # output CDS sequence
                        output_handle3 = open(
                            os.path.join(in_dir, assem_name + ".faa"),
                            "a")
                        output_handle3.write(
                            ">%s_%s_%s %s %s %s %s %s %s\n%s\n" %
                            (assem_name, record.id, desc, locus, ID, product, type, start, end, aa_seq))
                        output_handle3.close()
                        AA_true = 1
            output_handle1.close()
        print 'Retrieving whole genome sequences!'
        print 'Done!'
