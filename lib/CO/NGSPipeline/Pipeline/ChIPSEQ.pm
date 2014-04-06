package CO::NGSPipeline::Pipeline::ChIPSEQ;

use strict;
use File::Basename;
use base qw/CO::NGSPipeline/;

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
	
	my $use_convey = rand(1) > 0.8 ? 1 : 0;
	for(my $i = 0; $i < scalar(@$r1); $i ++) {
		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
		
		
		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_chipseq_fastqc_r1_$i");
		$qid->{fastqc_r1}->[$i] = $pipeline->bsmap->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);
		
		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_chipseq_fastqc_r2_$i");
			$qid->{fastqc_r2}->[$i] = $pipeline->bsmap->fastqc(
				fastq      => $r2_fastq,
				output_dir => "$pm->{dir}/fastqc_r2_$i"
			);
		}
		
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_chipseq_trimmed_$i");
		$qid->{trim} = $pipeline->bsmap->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output2 => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
		);

				
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_bsmap_chipseq_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->bsmap->fastqc(
			fastq      => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r1_trimmed_$i",
		);

		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_chipseq_fastqc_r2_trimmed_$i");
			$pm->set_job_dependency($qid->{trim});
			$qid->{fastqc_r2_trimmed} = $pipeline->bsmap->fastqc(
				fastq      => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i",
			);
		}
		
		###################################################################
		# alignment
		###################################################################
		$pm->set_job_name("$sample_id"."_chipseq_alignment_r1_$i");
		$pm->set_job_dependency($qid->{fastqc_r1_trimmed});
		$qid->{alignment_r1} = $pipeline->bsmap->bwa_aln(
			fastq        => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			genome       => "$GENOME_HG19",
			output       => "$pm->{dir}/$sample_id.r1.$i.sai",
			delete_input => 0,
			use_convey   => $use_convey,
		);
		
		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_chipseq_alignment_r2_$i");
			$pm->set_job_dependency($qid->{fastqc_r2_trimmed});
			$qid->{alignment_r2} = $pipeline->bsmap->bwa_aln(
				fastq        => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				genome       => "$GENOME_HG19",
				output       => "$pm->{dir}/$sample_id.r2.$i.sai",
				delete_input => 0,
				use_convey   => $use_convey,
			);	
		}
		
		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_chipseq_sampe_$i");
			$pm->set_job_dependency($qid->{alignment_r1}, $qid->{alignment_r2});
			$qid->{sampe} = $pipeline->bsmap->sampe(
				aln1         => "$pm->{dir}/$sample_id.r1.$i.sai",
				aln2         => "$pm->{dir}/$sample_id.r2.$i.sai",
				fastq1       => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
				fastq2       => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				genome       => "$GENOME_HG19",
				output       => "$pm->{dir}/$sample_id.$i.sorted.bam",
				delete_input => 1,
			);	
		} else {
			$pm->set_job_name("$sample_id"."_chipseq_samse_$i");
			$pm->set_job_dependency($qid->{alignment_r1});
			$qid->{sampe} = $pipeline->bsmap->samse(
				aln1         => "$pm->{dir}/$sample_id.r1.$i.sai",
				fastq1       => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
				genome       => "$GENOME_HG19",
				output       => "$pm->{dir}/$sample_id.$i.sorted.bam",
				delete_input => 1,
			);	
		}
		
		
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_chipseq_flagstat_$i");
		$pm->set_job_dependency($qid->{sampe});
		$qid->{flagstat}->[$i] = $pipeline->bsmap->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.sorted.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0
		);

		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.sorted.bam";
	}

	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_chipseq_merge_and_nodup");
	$pm->set_job_dependency(@{$qid->{flagstat}});
	$qid->{remove_duplicate} = $pipeline->bsmap->merge_nodup(
		sam_list     => $sam_sort_list,
		output       => "$pm->{dir}/$sample_id.nodup.bam",
		library      => $library,
		delete_input => 1
	);
									   
	########################################################################
	# insert size
	########################################################################
	$pm->set_job_name("$sample_id"."_chipseq_insertsize");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{insertsize} = $pipeline->bsmap->picard_insertsize(
		sam => "$pm->{dir}/$sample_id.nodup.bam"
	);

	
}

1;
