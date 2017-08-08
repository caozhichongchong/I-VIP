# Integron Visualization and Identification Pipeline (I-VIP)
#Copyright: LG209, Environmental Biotechnology Lab, The University of Hong Kong
#Author: An Ni ZHANG, email: caozhichongchong@gmail.com
#Citation

#I-VIP is a well-organized and user-friendly pipeline for integron identification, extraction, annotation and visualization
#Requirement: database folder containing all attC, integrase and sulI databases
#Requirement: scripts folder containing all scripts
#Requirement: both sequences and orfs of input data
#Attention: It's recommended to remove "_" and ":" from input data name
#Put I-VIP.py, scripts folder, database folder under the same directory of your input dataset(s)
python I-VIP.py -fasta .fa -orf .faa --gbff None --orftype 1  --module 1 --usearch usearch --quick 1

#usage: I-VIP.py [-h] [-fasta .fa] [-orf .faa] [--gbff gbff] [--orftype 0 or 1]
#                [--distance 4000] [--cutoff 1.0] [--module 0 or 1]
#                [--quick mid] [--usearch version of usearch] [--thread 15]

#optional arguments:
#  -h, --help            show this help message and exit
#  -fasta .fa            format of fasta sequences
#  -orf .faa             format of ORFs
#  --gbff gbff           Set the directory of Genbank files (default: gbff) No
#                        gbff input: None
#  --orftype 0 or 1      The orf type, eg: 0 for genbank parsing; 1 for
#                        prodigal prediction
#  --distance 4000       Distance cutoff for two cassettes to be clustered
#                        together (default is 4000)
#  --cutoff 1.0          E-value cutoff for attC site identification (default
#                        is 1.0)
#  --module 0 or 1       Alternative attC search strategies for Module A (0:
#                        global search by Module A1; 1: local search by Module
#                        A2), (default 0)
#  --quick mid           attC search and filtering setting (from least strict
#                        to most strict: 'max', 'nohmm', 'mid', '', 'rfam'),
#                        (default 'mid')
#  --usearch version of usearch
#                        Use usearch ahead of blastp for integrase and sulI
#                        alignment, 'None' for not using usearch, 'usearch' or
#                        other version of usearchfor using usearch (default:
#                        'None')
#  --thread 15           The thread number assigned for running I-VIP (default
#                        15)
