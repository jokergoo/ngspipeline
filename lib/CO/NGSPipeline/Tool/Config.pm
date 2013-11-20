package CO::NGSPipeline::Tool::Config;

#########################################
# common tools
#########################################

use strict;
use File::Spec;
use CO::Utils;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw($BWA
                 $CNYBWA
				 $PICARD_BIN_DIR
				 $TRIMPAIR_BIN_DIR
				 $BOWTIE_INDEX
				 $BOWTIE2_INDEX
				 $GENCODE_BOWTIE_INDEX
				 $GENCODE_BOWTIE2_INDEX
				 $GENCODE_GTF
				 $GENOME_HG19
				 $GENOME_1KG
                 );

our $BWA = 'bwa-0.6.2-tpx';
our $CNYBWA = 'cnybwa-0.6.2';
our $PICARD_BIN_DIR = '/ibios/tbi_cluster/11.4/x86_64/picard/picard-tools';

our $TRIMPAIR_BIN_DIR = '/ibios/co02/guz/program/general/pairTrim';


our $BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1_hg19_chr/hg19_1-22_X_Y_M';
our $BOWTIE2_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie2_hg19_chr/hg19_1-22_X_Y_M';

# gencode
our $GENCODE_BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1/gencode17/gencode.v17.annotation';
our $GENCODE_BOWTIE2_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie2/gencode17/gencode.v17.annotation';
#our $GENCODE_GTF = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode17/gencode.v17.annotation.gtf';
#our $GENCODE_GTF = '/ibios/raid6/users/herrmanc/Projects/LSC_Simon/output/merge_GTF/GENCODE_Broad.gtf';
our $GENCODE_GTF = '/ibios/co02/guz/project/hipo16/analysis/gencode/gencode_broad_merged.gtf';

# genome
our $GENOME_HG19 = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/hg19_chr/hg19_karyotypically_sorted.fa';
our $GENOME_1KG = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa';


1;