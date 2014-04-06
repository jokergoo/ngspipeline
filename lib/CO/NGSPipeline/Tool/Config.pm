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
                 $FASTQC
                 $CUTADAPT
                 $SAMTOOLS
				 $PICARD_BIN_DIR
				 $TRIMPAIR_BIN_DIR
				 $BOWTIE_INDEX
				 $BOWTIE2_INDEX
				 $GENCODE_BOWTIE_INDEX
				 $GENCODE_BOWTIE2_INDEX
				 $GENCODE_GTF
				 $GENOME_HG19
				 $GENOME_1KG
				 $GENCODE_DEXSEQ_GTF
                 );

our $PICARD_BIN_DIR = '/ibios/tbi_cluster/11.4/x86_64/picard/picard-tools';

our $TRIMPAIR_BIN_DIR = '/ibios/co02/guz/program/general/pairTrim';


our $BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1_hg19_chr/hg19_1-22_X_Y_M';
our $BOWTIE2_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie2_hg19_chr/hg19_1-22_X_Y_M';

# gencode
our $GENCODE_BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1/gencode17/gencode.v17.annotation';
our $GENCODE_BOWTIE2_INDEX = '/icgc/lsdf/mb/analysis/guz/gencode_index/gencode.v17.annotation';
our $GENCODE_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_broad_merged.gtf';
our $GENCODE_DEXSEQ_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_v17_broad_for_dexseq.gff';

#our $GENCODE_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_v19_merged.gtf';
#our $GENCODE_DEXSEQ_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_v19_for_dexseq.gff';


# genome
our $GENOME_HG19 = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/hg19_chr/hg19_karyotypically_sorted.fa';
our $GENOME_1KG = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa';


############################
## absolute path for binary programs
our $FASTQC = "/ibios/tbi_cluster/11.4/x86_64/bin/fastqc";
our $BWA = '/ibios/tbi_cluster/11.4/x86_64/bin/bwa-0.6.2-tpx';
our $CNYBWA = 'cnybwa-0.6.2';
our $CUTADAPT = '/ibios/tbi_cluster/11.4/x86_64/bin/cutadapt';
our $SAMTOOLS = '/ibios/tbi_cluster/11.4/x86_64/bin/samtools';
1;