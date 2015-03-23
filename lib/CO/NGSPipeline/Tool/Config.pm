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

		 $STAR_GENOME
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
                 $TOPHAT2

                 $BISMARK_BIN_DIR
                 $BSMAP_BIN_DIR
                 $BISSNP_BIN_DIR
                 $BISMARK_GENOME_DIR
		 $BISMARK_REF_GENOME
                 $BSMAP_GENOME_DIR
                 $BSMAP_REF_GENOME
                 $BISSNP_INTERVAL_FILE
                 $BISSNP_INDEL_1_FILE
                 $BISSNP_INDEL_2_FILE
                 $BISSNP_DBSNP_FILE
                 
                 $METHYLCTOOLS_BIN_DIR
                 $METHYLCTOOLS_REFERENCE_POS_FILE
                 $METHYLCTOOLS_GENOME_DIR
                 $METHYLCTOOLS_REF_GENOME
                 $METHYLCTOOLS_REF_GENOME_CONV
				 
		 $GENOME_DIR
		 $REF_GENOME

                $BSMAP_REF_GENOME_MM
                $REF_GENOME_MM
                $BISSNP_INTERVAL_FILE_MM
                $BISSNP_INDEL_1_FILE_MM
                $BISSNP_INDEL_2_FILE_MM
                $BISSNP_DBSNP_FILE_MM
                 );

our $PICARD_BIN_DIR = '/ibios/tbi_cluster/11.4/x86_64/picard/picard-tools';

our $TRIMPAIR_BIN_DIR = '/ibios/co02/guz/program/general/pairTrim';


our $BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1_hg19_chr/hg19_1-22_X_Y_M';
our $BOWTIE2_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie2_hg19_chr/hg19_1-22_X_Y_M';

# gencode
our $GENCODE_BOWTIE_INDEX = '/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bowtie/bowtie1/gencode17/gencode.v17.annotation';
our $GENCODE_BOWTIE2_INDEX = '/icgc/lsdf/mb/analysis/guz/gencode_index/gencode.v17.annotation';
#our $GENCODE_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_broad_merged.gtf';
#our $GENCODE_DEXSEQ_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_v17_broad_for_dexseq.gff';

our $GENCODE_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_v19_merged.gtf';
our $GENCODE_DEXSEQ_GTF = '/ibios/co02/guz/MethSuite/gencode/gencode_v19_for_dexseq.gff';

#our $GENCODE_GTF = '/ibios/raid6/users/falcone/Projects/lincRNA/data/annotations/gencode.v21_Broad_no_overlap.gtf';


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



#our $STAR_GENOME = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/STAR_Gencode_Broad/genome/';
our $STAR_GENOME = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/STAR_gencode19/';
#our $STAR_GENOME = '/icgc/lsdf/mb/analysis/guz/genome/RNAseq_genome/STAR_gencode21_broad/';

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


our $BISMARK_BIN_DIR = '/ibios/co02/guz/program/BSTools/bismark_bin';  # relative path to config.pm
our $BSMAP_BIN_DIR   = '/ibios/co02/guz/program/BSTools/bsmap_bin';
our $BISSNP_BIN_DIR  = '/ibios/co02/guz/program/BSTools/bissnp_bin';

# if other genome files are used, you should add lambda genome and the name for lambda genome should be 'lambda'

our $BISMARK_GENOME_DIR = '/icgc/lsdf/mb/analysis/guz/genome/WGBS_genome/genome_bismark_hg19';
our $BISMARK_REF_GENOME = 'hg19.fa';


our $BSMAP_GENOME_DIR = '/icgc/lsdf/mb/analysis/guz/genome/WGBS_genome/genome_protocol';  # add lambda genome
our $BSMAP_REF_GENOME = 'hg19.fa';
our $BSMAP_REF_GENOME_MM = 'mm9_lambda.fa';

our $GENOME_DIR = '/icgc/lsdf/mb/analysis/guz/genome/WGBS_genome/genome_protocol';
our $REF_GENOME = 'hg19.fa';
our $REF_GENOME_MM = 'mm9_lambda.fa';

our $BISSNP_INTERVAL_FILE = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/whole_genome_interval_list.hg19.bed';
our $BISSNP_INDEL_1_FILE  = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/1000G_phase1.indels.hg19.sort.vcf';
our $BISSNP_INDEL_2_FILE  = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf';
our $BISSNP_DBSNP_FILE    = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/dbsnp_135.hg19.sort.vcf';

our $BISSNP_INTERVAL_FILE_MM = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/whole_genome_interval_list.mm9.bed';
our $BISSNP_INDEL_1_FILE_MM  = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/mouse-20110602-callable-dinox-indels.annot.mm9.vcf';
our $BISSNP_INDEL_2_FILE_MM  = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/mouse-20110602-callable-dinox-indels.annot.mm9.vcf';
our $BISSNP_DBSNP_FILE_MM    = '/icgc/lsdf/mb/analysis/guz/genome/bissnp_files/mouse-20111102-snps-all.annotated.mm9.vcf';

our $METHYLCTOOLS_BIN_DIR            = '/ibios/co02/guz/program/BSTools/methylCtools_bin';
our $METHYLCTOOLS_REFERENCE_POS_FILE = 'hg19.reference.pos.gz';
our $METHYLCTOOLS_GENOME_DIR         = '/icgc/lsdf/mb/analysis/guz/genome/WGBS_genome/genome_protocol';
our $METHYLCTOOLS_REF_GENOME         = 'hg19.fa';
our $METHYLCTOOLS_REF_GENOME_CONV    = 'hg19.conv.fa';

1;