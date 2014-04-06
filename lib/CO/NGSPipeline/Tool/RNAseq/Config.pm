package CO::NGSPipeline::Tool::RNAseq::Config;

use strict;
use File::Spec;
use CO::Utils;
use File::Basename;
use CO::NGSPipeline::Tool::Config;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw($BOWTIE_INDEX
                 $BOWTIE2_INDEX
                 $STAR_GENOME
                 $GENCODE_BOWTIE_INDEX
                 $GENCODE_BOWTIE2_INDEX
                 $GENCODE_GTF
                 $GENCODE_GTF_STAR
                 $GENCODE_DEXSEQ_GTF
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
                 $FUSIONHUNTER_DIR
                 $STAR
                 $TOPHAT2);

our $BOWTIE_INDEX = $CO::NGSPipeline::Tool::Config::BOWTIE_INDEX;
our $BOWTIE2_INDEX = $CO::NGSPipeline::Tool::Config::BOWTIE2_INDEX;

our $STAR_GENOME = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/STAR_Gencode_Broad/genome/';
#our $STAR_GENOME = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/STAR_gencode19/';

# gencode
our $GENCODE_BOWTIE_INDEX = $CO::NGSPipeline::Tool::Config::GENCODE_BOWTIE_INDEX;
our $GENCODE_BOWTIE2_INDEX = $CO::NGSPipeline::Tool::Config::GENCODE_BOWTIE2_INDEX;
our $GENCODE_GTF = $CO::NGSPipeline::Tool::Config::GENCODE_GTF;
our $GENCODE_DEXSEQ_GTF = $CO::NGSPipeline::Tool::Config::GENCODE_DEXSEQ_GTF;

# genome
our $GENOME_HG19 = $CO::NGSPipeline::Tool::Config::GENOME_HG19;
our $GENOME_1KG = $CO::NGSPipeline::Tool::Config::GENOME_1KG;

# gsnap using 1000genome as reference genome
our $GSNAP_GENOME_DIR = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/GNSAP_hg19/genome/';
our $GSNAP_GENOME = 'hg19';
our $GSNAP_IIT = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/gencode_v17.iit';

###################################################
# gene fusion settings
###################################################

our $FUSIONMAP_BIN_DIR = '/ibios/co02/guz/program/gene_fusion/FusionMap/FusionMap_2013-07-30/bin';
our $FUSIONMAP_DIR = '/ibios/co02/guz/program/gene_fusion/FusionMap';
our $DEFUSE_DIR = '/ibios/co02/guz/program/gene_fusion/deFuse/defuse-0.6.1';
our $DEFUSE_CONFIG_FILE = '/ibios/co02/guz/program/gene_fusion/deFuse/defuse-0.6.1/scripts/config.txt';
our $FUSIONHUNTER_BIN_DIR = '/ibios/co02/guz/program/gene_fusion/FusionHunter/FusionHunter-v1.4/bin';
our $FUSIONHUNTER_DIR = '/ibios/co02/guz/program/gene_fusion/FusionHunter';


our $STAR = '/ibios/tbi_cluster/11.4/x86_64/bin/STAR-2.3.0e';
our $TOPHAT2 = '/ibios/tbi_cluster/11.4/x86_64/bin/tophat2';
1;
