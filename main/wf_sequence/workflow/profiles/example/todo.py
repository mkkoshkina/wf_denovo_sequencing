'''
EXAMPLE:
SEQS.append(DenovoSequencing(project="TET3", 
                            MS_data= '/mnt/raid0/immunotet/data/immunopeptidomics/TET3', 
                            genome='GRCh38.p14'
                            ))
SEQS.append(DenovoSequencing(project="INI287", 
                            MS_data= '/mnt/raid0/data_for_denovo_IP/INI287', 
                            genome='GRCh38.p14',
                            rnaseq= '/mnt/raid0/data_for_denovo_IP/INI287/unique_RNA_PDX159Aligned_TE.fa'
                            ))
'''

SEQS.append(DenovoSequencing(project=<name of the project>, 
                            MS_data= <path to the MS data>, 
                            genome=<name of the genome that you have indexed with IndexGenome.jar on setup>,
                            model=<name of the Casanovo model you want to use>, #optional, default is casanovo_nontryptic.ckpt
                            rnaseq= <path to the RNAseq data>, # optional, Fasta file containing extra transcript sequences (will be 3-frame translated!)!
                            extra= <path to the DNAseq data>, #optional, Fasta file containing extra sequences (DNA will be 6-frame translated!)!
                            var=<path to the var file> #optional, File containing variants (each line); could also be in ${input%.csv}.var, in which case you do not have to specify this parameter!
                            ))