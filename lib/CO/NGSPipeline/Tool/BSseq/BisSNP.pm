package CO::NGSPipeline::Tool::BSseq::BisSNP;

use strict;
use CO::NGSPipeline::Tool::BSseq::Config;
use CO::Utils;

use base qw(CO::NGSPipeline::Tool::BSseq::Common
            CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
	return bless $self, $class;
}

sub call_methylation {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	              "genome" => "$GENOME_DIR/$REF_GENOME",
				  "delete_input" => 0,
				  @_);
	
	my $bam = to_abs_path($param{bam});
	my $genome = to_abs_path($param{genome});
	my $delete_input = $param{delete_input};
	
	my $prefix = $bam;
	$prefix =~s/\.bam$//;
	
	my $pm = $self->get_pipeline_maker;
	my $common_name = $pm->get_job_name;
	my $qid;
	
	# indel realignment
	$pm->set_job_name($common_name."_add_RG");
	$pm->add_command("JAVA_OPTIONS=-Xmx10G picard.sh AddOrReplaceReadGroups INPUT=$bam OUTPUT=$prefix.withRG.bam RGID=readGroup_name RGLB=readGroup_name RGPL=illumina RGPU=run RGSM=sample_name SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_add_RG",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "15GB",
										walltime => "10:00:00"});

	
	# it only use ~3.5 cores
	$pm->set_job_name($common_name."_BisulfiteRealignerTargetCreator");
	$pm->set_job_dependency($qid);
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -I $prefix.withRG.bam -T BisulfiteRealignerTargetCreator -L $BISSNP_INTERVAL_FILE -known $BISSNP_INDEL_1_FILE -known $BISSNP_INDEL_2_FILE -o $prefix.indel_target_interval.intervals -nt 4");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteRealignerTargetCreator",
							 "-l" => { nodes => "1:ppn=4:lsdf", 
									    mem => "15GB",
										walltime => "40:00:00"});
	
	$pm->set_job_name($common_name."_BisulfiteIndelRealigner");
	$pm->set_job_dependency($qid);	
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -I $prefix.withRG.bam -T BisulfiteIndelRealigner -targetIntervals $prefix.indel_target_interval.intervals -known $BISSNP_INDEL_1_FILE -known $BISSNP_INDEL_2_FILE -cigar -o $prefix.realigned.bam");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteIndelRealigner",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "15GB",
										walltime => "6:00:00"});
	# it only use ~6 cores
	$pm->set_job_name($common_name."_BisulfiteCountCovariates");
	$pm->set_job_dependency($qid);
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -I $prefix.realigned.bam -T BisulfiteCountCovariates -knownSites $BISSNP_DBSNP_FILE -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -recalFile $prefix.recalFile_before.csv -nt 6");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteCountCovariates",
							 "-l" => { nodes => "1:ppn=6:lsdf", 
									    mem => "15GB",
										walltime => "20:00:00"});
										
	$pm->set_job_name($common_name."_BisulfiteTableRecalibration");
	$pm->set_job_dependency($qid);
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -I $prefix.realigned.bam -o $prefix.realigned.recal.bam -T BisulfiteTableRecalibration -recalFile $prefix.recalFile_before.csv -maxQ 40");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteTableRecalibration",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "15GB",
										walltime => "15:00:00"});
	
	$pm->set_job_name($common_name."_BisulfiteCountCovariates_after");
	$pm->set_job_dependency($qid);	
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -I $prefix.realigned.recal.bam -T BisulfiteCountCovariates -knownSites $BISSNP_DBSNP_FILE -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -recalFile $prefix.recalFile_after.csv -nt 6");
	
	$pm->add_command("java -Xmx4g -jar $BISSNP_BIN_DIR/BisulfiteAnalyzeCovariates.jar -recalFile $prefix.recalFile_before.csv -outputDir $pm->{dir}/bissnp_recal_before -ignoreQ 5 --max_quality_score 40");
	
	$pm->add_command("java -Xmx4g -jar $BISSNP_BIN_DIR/BisulfiteAnalyzeCovariates.jar -recalFile $prefix.recalFile_after.csv -outputDir $pm->{dir}/bissnp_recal_after -ignoreQ 5 --max_quality_score 40");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteCountCovariates_after",
							 "-l" => { nodes => "1:ppn=6:lsdf", 
									    mem => "15GB",
										walltime => "24:00:00"});
	
	# genotyping
	# note the -stand_call_conf option
	$pm->set_job_name($common_name."_BisulfiteGenotyper");
	$pm->set_job_dependency($qid);
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -T BisulfiteGenotyper -I $prefix.realigned.recal.bam -D $BISSNP_DBSNP_FILE -vfn1 $prefix.cpg.raw.vcf -vfn2 $prefix.snp.raw.vcf -stand_call_conf 10 -stand_emit_conf 0 -mmq 30 -mbq 0 -nt 12");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteGenotyper",
							 "-l" => { nodes => "1:ppn=12:lsdf", 
									    mem => "15GB",
										walltime => "200:00:00"});
										
	$pm->set_job_name($common_name."_BisulfiteVCFsort");
	$pm->set_job_dependency($qid);
	$pm->add_command("perl $BISSNP_BIN_DIR/sortByRefAndCor.pl --k 1 --c 2 --tmp $pm->{tmp_dir} $prefix.cpg.raw.vcf $genome.fai | grep 'lambda' -v > $prefix.cpg.sorted.vcf");
	$pm->add_command("perl $BISSNP_BIN_DIR/sortByRefAndCor.pl --k 1 --c 2 --tmp $pm->{tmp_dir} $prefix.snp.raw.vcf $genome.fai | grep 'lambda' -v > $prefix.snp.sorted.vcf");
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_BisulfiteVCFsort",
							 "-l" => { nodes => "1:ppn=3:lsdf", 
									    mem => "15GB",
										walltime => "100:00:00"});									
										
	# post process
	$pm->set_job_name($common_name."_VCFpostprocess");
	$pm->set_job_dependency($qid);
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -T VCFpostprocess -oldVcf $prefix.snp.sorted.vcf -newVcf $prefix.snp.filtered.vcf -snpVcf $prefix.snp.sorted.vcf -o $prefix.snp.sorted.filter.summary.txt");
	
	$pm->add_command("java -Xmx10g -jar $BISSNP_BIN_DIR/BisSNP.jar -R $genome -T VCFpostprocess -oldVcf $prefix.cpg.sorted.vcf -newVcf $prefix.cpg.filtered.vcf -snpVcf $prefix.cpg.sorted.vcf -o $prefix.cpg.sorted.filter.summary.txt");
	
	$pm->add_command("perl $BISSNP_BIN_DIR/vcf2wig.pl $prefix.cpg.filtered.vcf CG");
	$pm->add_command("perl $BISSNP_BIN_DIR/vcf2bedGraph.pl $prefix.cpg.filtered.vcf CG");
	$pm->add_command("perl $BISSNP_BIN_DIR/vcf2bed.pl $prefix.cpg.filtered.vcf CG");
	$pm->add_command("perl $BISSNP_BIN_DIR/vcf2coverage.pl $prefix.cpg.filtered.vcf CG");
	
	$pm->del_file("$prefix.withRG.bam");
	$pm->del_file("$prefix.withRG.bai");
	$pm->del_file("$prefix.indel_target_interval.intervals");
	$pm->del_file("$prefix.realigned.bam");
	$pm->del_file("$prefix.realigned.bai");
	$pm->del_file("$prefix.recalFile_before.csv");
	$pm->del_file("$prefix.recalFile_after.csv");
	$pm->del_file("$prefix.realigned.recal.bam");
	$pm->del_file("$prefix.realigned.recal.bai");
	$pm->del_file("$prefix.cpg.raw.vcf");
	$pm->del_file("$prefix.snp.raw.vcf");
	
	$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bissnp_VCFpostprocess",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "15GB",
										walltime => "10:00:00"});

	return($qid);

}

1;
