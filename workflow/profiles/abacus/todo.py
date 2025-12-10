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

SEQS.append(DenovoSequencing(project="INI305", 
                             MS_data="/mnt/beegfs/home/mveyssie/stagein/home/mveyssie/Documents/Projects/immuno_rhabdoid/data", 
                            genome="GRCh38.p14",
                            #model=<name of the Casanovo model you want to use>, #optional, default is casanovo_nontryptic.ckpt
                            rnaseq= "/mnt/beegfs/home/mveyssie/stagein/home/mveyssie/Documents/Projects/immuno_rhabdoid/data/PDX209_RNA_TE.fa", # optional, Fasta file containing extra transcript sequences (will be 3-frame translated!)!
                            #extra= <path to the DNAseq data>, #optional, Fasta file containing extra sequences (DNA will be 6-frame translated!)!
                            #var=<path to the var file> #optional, File containing variants (each line); could also be in ${input%.csv}.var, in which case you do not have to specify this parameter!
                            ))
