package CO::NGSPipeline::Pipeline::WGS;

use strict;
use File::Basename;
use base qw/CO::NGSPipeline/;

use CO::NGSPipeline::Tool::Config;

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
	return bless $self, $class;
}

sub run {
	
	my $pipeline = shift;
	
	my $pm = $pipeline->get_pipeline_maker;
	unless(defined($pm)) {
		die "pipeline maker should be defined in a pipeline.\n";
	}
	
	my %param = (sample_id => undef,
	             r1 => undef,
				 r2 => undef,
				 library => undef,
				 dir => undef,
				 @_);
				 
	my $sample_id = $param{sample_id};
	my $r1 = $param{r1};
	my $r2 = $param{r2};
	my $library = $param{library};
	my $dir = $param{dir};
	
	my $qid = {};
	$qid->{flagstat} = [];
	$qid->{fastqc_r1} = [];
	$qid->{fastqc_r2} = [];
	my $sam_sort_list = [];
	
	for(my $i = 0; $i < scalar(@$r1); $i ++) {
		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
		
	
		
		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_wgs_fastqc_r1_$i");
		$qid->{fastqc_r1}->[$i] = $pipeline->wgs->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);
		
		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_wgs_fastqc_r2_$i");
			$qid->{fastqc_r2}->[$i] = $pipeline->wgs->fastqc(
				fastq      => $r2_fastq,
				output_dir => "$pm->{dir}/fastqc_r2_$i"
			);
		}
	
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_wgs_trimmed_$i");
		$qid->{trim} = $pipeline->wgs->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output2 => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
		);

				
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_wgs_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->wgs->fastqc(
			fastq      => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r1_trimmed_$i",
		);

		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_wgs_fastqc_r2_trimmed_$i");
			$pm->set_job_dependency($qid->{trim});
			$qid->{fastqc_r2_trimmed} = $pipeline->wgs->fastqc(
				fastq      => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i",
			);
		}
	
		###################################################################
		# alignment
		###################################################################
		$pm->set_job_name("$sample_id"."_wgs_alignment_r1_$i");
		$pm->set_job_dependency($qid->{fastqc_r1_trimmed});
		$qid->{alignment_r1} = $pipeline->wgs->bwa_aln(
			fastq        => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			genome       => $GENOME_MM10,
			output       => "$pm->{dir}/$sample_id.r1.$i.sai",
			delete_input => 0,
		);

		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_wgs_alignment_r2_$i");
			$pm->set_job_dependency($qid->{fastqc_r2_trimmed});
			$qid->{alignment_r2} = $pipeline->wgs->bwa_aln(
				fastq        => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				genome       => $GENOME_MM10,
				output       => "$pm->{dir}/$sample_id.r2.$i.sai",
				delete_input => 0,
			);
		}

		##################################################################
		# sampe/samse
		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_wgs_sampe_$i");
			$pm->set_job_dependency($qid->{alignment_r1}, $qid->{alignment_r2});
			$qid->{samse_or_pe} = $pipeline->wgs->sampe(
				aln1 => "$pm->{dir}/$sample_id.r1.$i.sai",
				aln2 => "$pm->{dir}/$sample_id.r2.$i.sai",
				fastq1 => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
				fastq2 => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				genome => $GENOME_MM10,
				output => "$pm->{dir}/$sample_id.$i.sorted.bam",
				delete_input => 0,
				);
		} else {
			$pm->set_job_name("$sample_id"."_wgs_samse_$i");
			$pm->set_job_dependency($qid->{alignment_r1}, $qid->{alignment_r2});
			$qid->{samse_or_pe} = $pipeline->wgs->samse(
				aln1 => "$pm->{dir}/$sample_id.r1.$i.sai",
				fastq1 => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
				genome => $GENOME_MM10,
				output => "$pm->{dir}/$sample_id.$i.sorted.bam",
				delete_input => 0,
				);
		}
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_wgs_flagstat_$i");
		$pm->set_job_dependency($qid->{samse_or_pe});
		$qid->{flagstat}->[$i] = $pipeline->wgs->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.sorted.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0
		);

		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.sorted.bam";
	}

	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_wgs_merge_and_mkdup");
	$pm->set_job_dependency(@{$qid->{flagstat}});
	$qid->{remove_duplicate} = $pipeline->wgs->merge_nodup(
		sam_list     => $sam_sort_list,
		output       => "$pm->{dir}/$sample_id.mkdup.bam",
		library      => $library,
		REMOVE_DUPLICATES => 'FALSE',
		delete_input => 0
	);
	

	####################################################################
	# flagstat
	####################################################################
	$pm->set_job_name("$sample_id"."_wgs_flagstat_final");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$pipeline->wgs->samtools_flagstat(
		sam          => "$pm->{dir}/$sample_id.mkdup.bam",
		output       => "$pm->{dir}/$sample_id.mkdup.flagstat",
		delete_input => 0
	);	

	####################################################################
	# coverage
	####################################################################
	$pm->set_job_name("$sample_id"."_wgs_coverage");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$pipeline->wgs->coverage(
		sam          => "$pm->{dir}/$sample_id.mkdup.bam",
		sample_id   => "$sample_id",
		ungap        => "/icgc/ngs_share/assemblies/mm10/stats/GRCm38mm10.fa.chrLenOnlyACGT_realChromosomes.tab",
		output       => "$pm->{dir}/$sample_id.mkdup.coverage_qc",
		delete_input => 0
	);	
}

1;
