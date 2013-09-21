package CO::NGSPipeline::Tool::RNAseq::Config;

use strict;
use File::Spec;
use CO::Utils;
use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw($BOWTIE_INDEX
                 $BOWTIE2_INDEX
                 $STAR_GENOME
                 $GENCODE_BOWTIE_INDEX
                 $GENCODE_BOWTIE2_INDEX
                 $GENCODE_GTF
                 $GENCODE_GTF_STAR
                 $GENOME_1KG
                 $GENOME_HG19
                 $GSNAP_GENOME_DIR
                 $GSNAP_GENOME
                 $GSNAP_IIT
                 $FUSIONMAP_BIN_DIR
                 $FUSIONMAP_DIR
                 $DEFUSE_DIR
                 $DEFUSE_CONFIG_FILE
                 $FUSIONHUNTER_BIN_DIR
                 $FUSIONHUNTER_DIR);

our $BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1_hg19_chr/hg19_1-22_X_Y_M';
our $BOWTIE2_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie2_hg19_chr/hg19_1-22_X_Y_M';

our $STAR_GENOME = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/STAR/genome';

# gencode
our $GENCODE_BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1/gencode17/gencode.v17.annotation';
our $GENCODE_BOWTIE2_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie2/gencode17/gencode.v17.annotation';
our $GENCODE_GTF = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode17/gencode.v17.annotation.gtf';

# genome
our $GENOME_HG19 = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/hg19_chr/hg19_karyotypically_sorted.fa';
our $GENOME_1KG = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa';

# gsnap using 1000genome as reference genome
our $GSNAP_GENOME_DIR = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/GNSAP/genome/';
our $GSNAP_GENOME = '1KGRef';
our $GSNAP_IIT = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/gencode_v17_nochr.iit';

###################################################
# gene fusion settings
###################################################

our $FUSIONMAP_BIN_DIR = '/ibios/co02/guz/program/gene_fusion/FusionMap/FusionMap_2013-07-30/bin';
our $FUSIONMAP_DIR = '/ibios/co02/guz/program/gene_fusion/FusionMap';
our $DEFUSE_DIR = '/ibios/co02/guz/program/gene_fusion/deFuse/defuse-0.6.1';
our $DEFUSE_CONFIG_FILE = '/ibios/co02/guz/program/gene_fusion/deFuse/defuse-0.6.1/scripts/config.txt';
our $FUSIONHUNTER_BIN_DIR = '/ibios/co02/guz/program/gene_fusion/FusionHunter/FusionHunter-v1.4/bin';
our $FUSIONHUNTER_DIR = '/ibios/co02/guz/program/gene_fusion/FusionHunter';

1;
