package CO::NGSPipeline::Tool::RNAseq::GeneFusion;

use strict;
use CO::NGSPipeline::Tool::RNAseq::Config;
use CO::Utils;
use File::Basename;

use base qw(CO::NGSPipeline::Tool::RNAseq::Common
            CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;

	my $self = {};
	
	return bless $self, $class;
}

sub fusionmap {
	my $self = shift;
	
	my %param = ( "fastq1" => undef,
	              "fastq2" => undef,
				  "sample_id" => "sample",
				  "delete_input" => 0,
				  @_);
	
	my $fastq1 = to_abs_path($param{fastq1});
	my $fastq2 = to_abs_path($param{fastq2});
	my $delete_input = $param{delete_input};
	my $sample_id = $param{sample_id};
	
	my $pm = $self->get_pipeline_maker;
	open OUT, ">$pm->{dir}/secontrol.txt";
	print OUT <<CONFIG;
<Files>
$pm->{dir}/fastq1
$pm->{dir}/fastq2

<Options>
PairedEnd=True
RnaMode=True
Use32BitMode=False
ThreadNumber=8 //Possible values: 1-100. Default value=1
FileFormat=FASTQ // Possible values: FASTQ, QSEQ, FASTA. Default value=FASTQ
MinimalFusionAlignmentLength=25 //Possible values: 15-50. Default value=25 (alpha)
FusionReportCutoff=1 // Possible values: 1-1000. Default value=1 (beta)
NonCanonicalSpliceJunctionPenalty=4 //Possible values: 01-. Default value = 2 (G)
MinimalHit=2 // Minimal distinct read; Possible values: 1-10000, Default value =2 
MinimalRescuedReadNumber=1 // Minimal rescued read number. Default value = 1
OutputFusionReads=True // Possible values: True, False. Default value = True
FilterBy=DefaultList //advanced filtering using default black list from FusionMap, set to None to avoid automatic downloading
MonoPath=/ibios/tbi_cluster/11.4/x86_64/bin/mono

<Output>
TempPath=$pm->{tmp_dir}
OutputPath=$pm->{dir}/results
OutputName=$sample_id

CONFIG
	close OUT;
	
	if($fastq1 =~/\.gz$/) {
		$pm->add_command("zcat -c $fastq1 > $pm->{dir}/fastq1", 0);
	} else {
		$pm->add_command("ln -s $fastq1 $pm->{dir}/fastq1", 0);
	}
	if($fastq2 =~/\.gz$/) {
		$pm->add_command("zcat -c $fastq2 > $pm->{dir}/fastq2", 0);
	} else {
		$pm->add_command("ln -s $fastq2 $pm->{dir}/fastq2", 0);
	}
	
	$pm->add_command("mono $FUSIONMAP_BIN_DIR/FusionMap.exe --semap $FUSIONMAP_DIR Human.B37 RefGene $pm->{dir}/secontrol.txt");
	$pm->add_command("rm $pm->{dir}/fastq1", 0);
	$pm->add_command("rm $pm->{dir}/fastq2", 0);
	$pm->del_file($fastq1, $fastq2) if($delete_input);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_fusionmap",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "10GB",
										walltime => "150:00:00"});

	return($qid);

}

sub defuse {
	my $self = shift;
	
	my %param = ( "fastq1_arrayref" => undef,
	              "fastq2_arrayref" => undef,
	              "sample_id" => "sample",
				  "delete_input" => 0,
				  @_);
	
	my $fastq1_arrayref = $param{fastq1_arrayref};
	my $fastq2_arrayref = $param{fastq2_arrayref};
	my $delete_input = $param{delete_input};
	my $sample_id = $param{sample_id};
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->del_file("$pm->{dir}/fastq1", "$pm->{dir}/fastq2");
	$pm->add_command("mkfifo $pm->{dir}/fastq1 $pm->{dir}/fastq2", 0);
	
	my $cmd1 = "perl /home/guz/project/development/ngs_pipeline/zcat.pl ";
	my $cmd2 = "perl /home/guz/project/development/ngs_pipeline/zcat.pl ";
	for(my $i = 0; $i < scalar(@$fastq1_arrayref); $i ++) {
		$cmd1 .= " ".to_abs_path($fastq1_arrayref->[$i]);
		$cmd2 .= " ".to_abs_path($fastq2_arrayref->[$i]);
	}
	$cmd1 .= " > $pm->{dir}/fastq1 &";
	$cmd2 .= " > $pm->{dir}/fastq2 &";
	$pm->add_command($cmd1, 0);
	$pm->add_command($cmd2, 0);
		
	$pm->add_command("perl $DEFUSE_DIR/scripts/defuse.pl -c $DEFUSE_CONFIG_FILE -1 $pm->{dir}/fastq1 -2 $pm->{dir}/fastq2 -o $pm->{dir}/results -p 16");
	$pm->add_command("rm $pm->{dir}/fastq1", 0);
	$pm->add_command("rm $pm->{dir}/fastq2", 0);
	
	$pm->add_command("perl /ibios/co02/guz/program/gene_fusion/deFuse/add_anno_to_defuse.pl $pm->{dir}/results/results.filtered.tsv > $pm->{dir}/$sample_id.gene_fusion.annotated.txt");
	$pm->add_command("cp $pm->{dir}/results/results.filtered.tsv $pm->{dir}/results.filtered.tsv", 0);
	$pm->add_command("chmod 755 $pm->{dir}/$sample_id.gene_fusion.annotated.txt");
	$pm->add_command("chmod 755 $pm->{dir}/results.filtered.tsv");
	$pm->add_command("rm -rf $pm->{dir}/results", 0);
	
	if($delete_input) {
		$pm->del_file(@$fastq1_arrayref);
		$pm->del_file(@$fastq2_arrayref);
	}
	
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_defuse",
							 "-l" => { nodes => "1:ppn=16:lsdf", 
									    mem => "50GB",
										walltime => "150:00:00"});

	return($qid);

}

sub tophatfusion {
	my $self = shift;
	
	my %param = ( "fastq1" => undef,
	              "fastq2" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $fastq1 = to_abs_path($param{fastq1});
	my $fastq2 = to_abs_path($param{fastq2});
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
		
	$pm->add_command("$TOPHAT2 -o $pm->{dir} -p 8 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM $BOWTIE_INDEX $fastq1 $fastq2");
	$pm->del_file($fastq1, $fastq2) if($delete_input);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_tophatfusion",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "10GB",
										walltime => "150:00:00"});

	return($qid);

}


sub fusionhunter {
	my $self = shift;
	
	my %param = ( "fastq1" => undef,
	              "fastq2" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $fastq1 = to_abs_path($param{fastq1});
	my $fastq2 = to_abs_path($param{fastq2});
	my $delete_input = $param{delete_input};

	my $pm = $self->get_pipeline_maker;
	-e "$pm->{dir}/results" ? 1: mkdir("$pm->{dir}/results");
	
	open OUT, ">$pm->{dir}/results/FusionHunter.cfg";
	print OUT <<CONFIG;
<Files>
########################
# version1.4
# modification 2012/06/09
##########################
# configuration file for FusionHunter
# '#' indicates comments
# '=' is in need for each configuration line
# edit this file in your working directory 
##########################

# hg19 annotation packages are added, replace all hg18 links with hg19 if you are using hg19

##########################
# Basic options, you must get them changed for your data!! 
##########################

#set 1 if running on hg18 or hg19; otherwise 0
IF_HUMAN = 1

# left part of read pairs, should be in fastq format, and named as XXXX/1
L=$pm->{dir}/fastq1

# right part of read pairs should be in fastq format, and named as XXXX/2
R=$pm->{dir}/fastq2

# reference dir/name, reference should be in fasta format. Only 'major' chromosomes shall be included in reference genome. Undetermined scaffolds must be excluded since they might lead to spurious fusion outputs. 
Reference=/icgc/lsdf/mb/analysis/guz/gene_fusion/TopHatFusion/bowtie_index/hg19.fa

# the directory containing Bowtie index/basename of Bowtie index, built from the Reference file you provided.
# NO '/' in the end
#

BowtieIdx=/icgc/lsdf/mb/analysis/guz/gene_fusion/TopHatFusion/bowtie_index/hg19

# the directory and name of gene annotation list, we suggest UCSC annotation, these AnnotationFiles are included in FusionHunter package. For non-human species, please download the annotation from UCSC table browser with first 10 columns# as the GenePred table format, and last column should be the gene name. This track is mandatory in FusionHunter.
Gene_annotation = /icgc/lsdf/mb/analysis/guz/gene_fusion/FusionHunter/AnnotationFiles_hg19/hg19.ucscKnownGene

# directory and file name repeats region annotation, these AnnotationFiles are included in FusionHunter package for hg18/hg19. Optional for non-human species in FusionHunter, leave it blank if not available.
Repeats = /icgc/lsdf/mb/analysis/guz/gene_fusion/FusionHunter/AnnotationFiles_hg19/hg19.repeats

# the directory and file name of self alignment regions, these AnnotationFiles are included in FusionHunter package for hg18/hg19. Optional for non-human species in FusionHunter, leave it blank if not available.
SelfAlign = /icgc/lsdf/mb/analysis/guz/gene_fusion/FusionHunter/AnnotationFiles_hg19/hg19.chain.pairs

# the directory and file name of human EST database, these AnnotationFiles are included in FusionHunter package for hg18/hg19. Optional for non-human species in FusionHunter, leave it blank if not available.
EST = /icgc/lsdf/mb/analysis/guz/gene_fusion/FusionHunter/AnnotationFiles_hg19/hg19.SpliceEST

# size of the segmented reads, we strongly suggest it should not be longer than half of full read length e.g. <=25 if your RNA-seq read is 50bp
segment_size = 25

# number of cores for bowtie and other parallel processes, for sake of speed, we suggest you use as many cores as possible
CORE = 8

#  min number of paired-end reads that support a fusion (used in regionPairsList)
#  the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get, and shall take longer time
PAIRNUM = 2

# min number of junction spanning reads to support a fusion
# the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get
#
# TO BE NOTED: if you set the MINSPAN = 1, in order to reduce false positives, any candidate junction supported 
# by only 1 spanning read would be discarded unless the fusion junction point is exactly on annotated exon boundary. This process is embeded in FusionHunter.
MINSPAN = 1

# min size of the maximum base coverage on either side of the junction
# the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get
MINOVLP = 8

# Size of exact match for each junction flanking tile, should not be larger than MINOVLP
# the larger this number, the more specificity you get; the smaller this number, the more sensitivity you get, and shall take longer time
# We strongly suggest it should be >=3
TILE = 4

# total mismatch
# the smaller this number, the more specificity you get; the larger this number, the more sensitivity you get, and shall take longer time
# We strongly suggest it should be <=4
MISMATCH = 2

#########################
# Advanced options, you may change them for your preference
# changes are optional
########################

# number of multi hits for partial reads
M1 = 20

# number of reads to keep for partial reads
K1 = 8

# number of multi hits for full reads
M2 = 1

# number of reads to keep for full reads
K2 = 1

# max allowed repeat proportion of a read (used in reduceBwt)
REAPTOVLP = 0.6

# number of chains to overlap with a read (used in reduceBwt)
CHAINNUM = 20

# max allowed repeat proportion of a read (more stringent, used in leftRightOvlp)
RPTOVLP = 0.2

# max allowed alignment proportion between a pair of reads (used in leftRightOvlp)
CHAINOVP = 0.2

# distance to self-chain boundary (used in postLeftRightOvlp)
CHAINDIS = 200000

# proportion of a overlaps with a region (used in regionPairs)
READOVLP = 0.8

# screen out repetitive regions of reference when doing gapped alignment
MASK = 1

# output of fusion by FusionHunter
fusion_output =$pm->{dir}/FusionHunter.fusion

# output of readthrough by FusionHunter
readthrough_output = $pm->{dir}/FusionHunter.readthrough

CONFIG
	close OUT;
	
	
	if($fastq1 =~/\.gz$/) {
		$pm->add_command("zcat -c $fastq1 > $pm->{dir}/fastq1", 0);
	} else {
		$pm->add_command("ln -s $fastq1 $pm->{dir}/fastq1", 0);
	}
	if($fastq2 =~/\.gz$/) {
		$pm->add_command("zcat -c $fastq2 > $pm->{dir}/fastq2", 0);
	} else {
		$pm->add_command("ln -s $fastq2 $pm->{dir}/fastq2", 0);
	}
	
	$pm->add_command("cd $pm->{dir}/results", 0);
 	$pm->add_command("perl /icgc/lsdf/mb/analysis/guz/gene_fusion/FusionHunter/FusionHunter-v1.4/bin/FusionHunter.pl FusionHunter.cfg");
	$pm->add_command("rm $pm->{dir}/fastq1", 0);
	$pm->add_command("rm $pm->{dir}/fastq2", 0);
	$pm->del_file($fastq1, $fastq2) if($delete_input);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_fusionhunter",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "20GB",
										walltime => "150:00:00"});

	return($qid);

}

1;
