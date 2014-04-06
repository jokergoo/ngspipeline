package CO::NGSPipeline::Tool::Common;

##############################################################################
# provide command method for each pipeline. It is a base module for all specific
# pipelines.

use strict;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;
use List::Vectorize;
use File::Basename;

use base qw/CO::NGSPipeline/;

sub fastqc {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	my %param = ( "fastq" => undef,
	              "output_dir" => $pm->{dir},
	              @_);

	my $fastq      = to_abs_path( $param{fastq} );
	my $output_dir = to_abs_path( $param{output_dir} );

	if(! -e $output_dir) {
		$pm->add_command("mkdir -p $output_dir", 0);
	}
	$pm->add_command("$FASTQC -o $output_dir $fastq");
	my $f = basename($fastq);
	$f =~s/\.(fastq|fq)(\.gz)?$/_fastqc.zip/i;
	$pm->check_filesize("$output_dir/$f", 10*1024); # 10K
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_fastqc",
					          "-l" => { nodes => "1:ppn=1:lsdf",
							            mem => "1GB",
							            walltime => "30:00:00"},);
	return($qid);
	
}

=item C<$self-E<gt>trim(HASH)>

Trim FALSTQ files. FASTQ file can be either gziped or not

  fastq1        path of read 1
  fastq2        path of read 2
  output1       path of output 1
  output2       path of output 2
  output_dir    dir for the trimmed data
  delete_input  whether delete input files
 
=cut
sub trim {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	my %param = ( "fastq1" => undef,
	              "fastq2" => undef,
				  "output1" => undef,
				  "output2" => undef,
				  "polya" => 0,
				  "delete_input" => 0,
				  @_);
	
	my $fastq1  = to_abs_path( $param{fastq1} );
	my $fastq2  = to_abs_path( $param{fastq2} );
	my $output1 = to_abs_path( $param{output1} );
	my $output2 = to_abs_path( $param{output2} );
	my $polya = $param{polya};  # now it is not supported
	
	my $delete_input = $param{delete_input};
	
	my $qid;
	if($fastq2) {
		$pm->add_command("perl $TRIMPAIR_BIN_DIR/trim.pl --fastq1=$fastq1 --fastq2=$fastq2 --output1=$output1 --output2=$output2 --tmp=$pm->{tmp_dir}");
		$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_trim",
								  "-l" => { nodes => "1:ppn=3:lsdf",
											mem => "1GB",
											walltime => "50:00:00"},);
	} else {
		$pm->add_command("$CUTADAPT $fastq1 --quality-base 33 --quality_cutoff 20 --adapter AGATCGGAAGAGC ---minimum-length 20 | gzip -c > $output1");
		$qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_trim",
								  "-l" => { nodes => "1:ppn=1:lsdf",
											mem => "1GB",
											walltime => "50:00:00"},);
	}
	return($qid);
}

=item C<$self-E<gt>sort_sam(HASH)>

Sort SAM or BAM files

  same          sam file or bam file
  output        path of output
  delete_input  whether delete input files
 
=cut

# sort_sam would be used by pipe
sub sort_sam {
	my $self = shift;
	
	my %param = ( "sam" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  "sort_by" => "coordinate",
				  "add_index" => 0,
				  @_);
	
	my $sam    = to_abs_path( $param{sam} );
	my $output = to_abs_path( $param{output} );
	my $delete_input = $param{delete_input};
	my $sort_by = $param{sort_by};
	my $add_index = $param{add_index};
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("JAVA_OPTIONS=-Xmx50G picard.sh SortSam INPUT=$sam OUTPUT=$output SORT_ORDER=$sort_by TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=50000000");
	if($output =~/\.bam$/ and $add_index) {
		$pm->add_command("$SAMTOOLS index $output");
	}
	$pm->check_filesize($output);
	$pm->del_file($sam) if($delete_input);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_sort_sam",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "50GB",
										walltime => "20:00:00"});
	return($qid);
}

=item C<$self-E<gt>samtools_view(HASH)>

convert between SAM and BAM

  input         path of input
  output        path of output
  delete_input  whether delete input files
 
=cut
sub samtools_view {
	my $self = shift;
	
	my %param = ( "input" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $input    = to_abs_path( $param{input} );
	my $output = to_abs_path( $param{output} );
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	if($input =~/\.sam$/ and $output =~/\.bam/) {
		$pm->add_command("$SAMTOOLS view -Sbh $input -o $output");
	} elsif($input =~/\.bam$/ and $output =~/\.sam/) {
		$pm->add_command("$SAMTOOLS view -h $input -o $output");
	} else {
		die "Wrong extended file name.\n";
	}
	
	
	$pm->del_file($input) if($delete_input);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_samtools_view",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "20:00:00"});
	return($qid);

}

=item C<$self-E<gt>merge_nodup(HASH)>

Merge and remove duplications from SAM/BAM files

  sam_list      list of SAM/BAM files, array reference
  output        path of output
  same_library  whether these multiple lanes from same library
  delete_input  whether delete input files
 
=cut
sub merge_nodup {
	my $self = shift;
	
	my %param = ( "sam_list" => [],
	               "output" => undef,
				   "library" => undef,
				   "delete_input" => 0,
				   "sort_by" => "coordinate",
				   "REMOVE_DUPLICATES" => 'TRUE',
				   @_);
	
	
	my $sam_list = $param{sam_list};
	$sam_list = sapply($sam_list, \&to_abs_path);
	my $output = to_abs_path( $param{output} );
	my $library = $param{library}; $library = defined($library) ? $library : rep(1, len($sam_list));
	my $delete_input = $param{delete_input};
	my $sort_by = $param{sort_by};
	my $REMOVE_DUPLICATES = $param{REMOVE_DUPLICATES};
	
	$sam_list->[0] =~/\.(sam|bam)$/i;
	my $suffix = $1;
	
	my $pm = $self->get_pipeline_maker;
	
	if(scalar(@$sam_list) == 1) {
		my $sam_file;
		my $sam_sort_file;
		my $sam_nodup_file;
		my $sam_nodup_metric_file;
	
		$sam_file = $sam_list->[0];
		$sam_nodup_file = $output;
		$sam_nodup_metric_file = $output; $sam_nodup_metric_file =~s/\.(sam|bam)$/.mkdup.metrics/;
		$pm->add_command("JAVA_OPTIONS=-Xmx50G picard.sh MarkDuplicates INPUT=$sam_file OUTPUT=$output METRICS_FILE=$sam_nodup_metric_file TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=$REMOVE_DUPLICATES ASSUME_SORTED=TRUE CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=50000000");
		$pm->del_file($sam_file) if($delete_input);
		$pm->del_file("$sam_file.bai") if($delete_input);
	} else {
		# merge and remove duplicate in each library
		my $sam_nodup_file = tapply([0..$#$library], $library, sub {
			my @index = @_;
			my $library_subset_name = subset($library, \@index);
			$library_subset_name = $library_subset_name->[0];
			my $library_bam = [];
			
			my $sam_file;
			my $sam_sort_file;
			my $sam_nodup_file;
			my $sam_nodup_metric_file;
			
			if(scalar(@index) == 1) {
				$sam_file = $sam_list->[ $index[0] ];
				$sam_nodup_file = "$output";
				$sam_nodup_file =~s/\.(sam|bam)$/.library_$library_subset_name.$1/;
				$sam_nodup_metric_file = $sam_nodup_file; $sam_nodup_metric_file =~s/\.(sam|bam)$/.mkdup.metrics/;
				$pm->add_command("JAVA_OPTIONS=-Xmx50G picard.sh MarkDuplicates INPUT=$sam_file OUTPUT=$sam_nodup_file METRICS_FILE=$sam_nodup_metric_file TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=$REMOVE_DUPLICATES ASSUME_SORTED=TRUE CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=50000000");
				$pm->del_file($sam_file) if($delete_input);
				$pm->del_file("$sam_file.bai") if($delete_input);
				
			} else {
				my $input_str;
				for(my $i = 0; $i < scalar(@index); $i ++) {
					$input_str .= "INPUT=$sam_list->[ $index[$i] ] ";
				}
				$sam_sort_file = dirname($sam_list->[0])."/_tmp_$library_subset_name.".int(rand(999999)).".$suffix";
				$pm->add_command("JAVA_OPTIONS=-Xmx16G picard.sh MergeSamFiles $input_str OUTPUT=$sam_sort_file SORT_ORDER=$sort_by TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");
				$sam_nodup_file = "$output";
				$sam_nodup_file =~s/\.(sam|bam)$/.library_$library_subset_name.$1/;
				$sam_nodup_metric_file = $sam_nodup_file; $sam_nodup_metric_file =~s/\.(sam|bam)$/.mkdup.metrics/;
				$pm->add_command("JAVA_OPTIONS=-Xmx50G picard.sh MarkDuplicates INPUT=$sam_sort_file OUTPUT=$sam_nodup_file METRICS_FILE=$sam_nodup_metric_file TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=$REMOVE_DUPLICATES ASSUME_SORTED=TRUE CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=50000000");
				for(my $i = 0; $i < scalar(@index); $i ++) {
					$pm->del_file($sam_list->[ $index[$i] ]) if($delete_input);
					$pm->del_file("$sam_list->[ $index[$i] ].bai") if($delete_input);
				}
				$pm->del_file("$sam_sort_file");
			}
			
			return $sam_nodup_file;
		});
		$sam_nodup_file = [values %$sam_nodup_file];
		
		# finally merge samples from different libraries
		if(len($sam_nodup_file) == 1) {
			$pm->add_command("mv $sam_nodup_file->[0] $output", 0);
			my $sam_nodup_metric_file = $sam_nodup_file->[0]; $sam_nodup_metric_file =~s/\.(sam|bam)$/.mkdup.metrics/;
			my $output_metric_file = $output; $output_metric_file =~s/\.(sam|bam)$/.mkdup.metrics/;
			$pm->add_command("mv $sam_nodup_metric_file $output_metric_file", 0);
			my $sam_nodup_file_bai = $sam_nodup_file->[0]; $sam_nodup_file_bai=~s/\.bam/.bai/;
			my $output_bai = $output; $output_bai =~s/\.bam/.bai/;
			$pm->add_command("mv $sam_nodup_file_bai $output_bai", 0);
		} else {
			my $input_str;
			for(my $i = 0; $i < scalar(@$sam_nodup_file); $i ++) {
				$input_str .= "INPUT=$sam_nodup_file->[$i] ";
			}
			$pm->add_command("JAVA_OPTIONS=-Xmx16G picard.sh MergeSamFiles $input_str OUTPUT=$output SORT_ORDER=$sort_by TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");
			$pm->add_command("$SAMTOOLS index $output");
			for(my $i = 0; $i < scalar(@$sam_nodup_file); $i ++) {
				$pm->del_file($sam_nodup_file->[$i]);
				my $bai = $sam_nodup_file->[$i]; $bai =~s/\.bam/.bai/;
				$pm->del_file($bai);
			}
		}
	}
	$pm->check_filesize($output);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_merge_nodup",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "55GB",
										walltime => "200:00:00"});

	return($qid);
}

=item C<$self-E<gt>bwa_aln(HASH)>

bwa alignment

  fastq         fastq file
  genome        genome
  output        output file
  delete_input  whether delete input files
 
=cut
sub bwa_aln {
	my $self = shift;
	
	my %param = ( "fastq" => undef,
	              "genome" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  "use_convey" => 0,
				  @_);
	
	my $fastq = to_abs_path($param{fastq});
	my $genome = to_abs_path($param{genome});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	my $use_convey = $param{use_convey};
	
	my $pm = $self->get_pipeline_maker;
	
	if($use_convey) {
		$pm->add_command("$CNYBWA aln -t 12 -q 20 $genome $fastq > $output");
		$pm->del_file($fastq) if($delete_input);
		$pm->check_filesize($output); # 1M
		my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bwa_align",
								"-l" => { nodes => "1:ppn=12:lsdf", 
											mem => "10GB",
											walltime => "150:00:00"},
								"-q" => "convey",
								"-S" => "/bin/bash");

		return($qid);
	} else {
		$pm->add_command("$BWA aln -t 8 -q 20 $genome $fastq > $output");
		$pm->del_file($fastq) if($delete_input);
		$pm->check_filesize($output); # 1M
		my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bwa_align",
								"-l" => { nodes => "1:ppn=8:lsdf", 
											mem => "10GB",
											walltime => "150:00:00"},
								);
		return($qid);
	}
}


=item C<$self-E<gt>sampe(HASH)>

bwa pair-end

  aln1          alignment for read 1
  aln2          alignment for read 2
  fastq1        path of read 1
  fastq2        path of read 2
  genome        genome
  output        output
  delete_input  whether delete input files
 
=cut
sub sampe {
	my $self = shift;
	
	my %param = ( "aln1" => undef,
	              "aln2" => undef,
	              "fastq1" => undef,
	              "fastq2" => undef,
	              "genome" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $aln1 = to_abs_path($param{aln1});
	my $aln2 = to_abs_path($param{aln2});
	my $fastq1 = to_abs_path($param{fastq1});
	my $fastq2 = to_abs_path($param{fastq2});
	my $genome = to_abs_path($param{genome});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	
	if(!$fastq2) {
		return($self->samse(aln1 => $aln1,
		                    fastq1 => $fastq1,
							genome => $genome,
							output => $output,
							delete_input => $delete_input));
	}
	
	my $pm = $self->get_pipeline_maker;
	
	my $r = time().rand();
	$pm->add_command("$BWA sampe $genome $aln1 $aln2 $fastq1 $fastq2 | mbuffer -q -m 2G -l /dev/null | samtools view -uSbh - | mbuffer -q -m 2G -l /dev/null | samtools sort - $pm->{dir}/$r");
	$pm->add_command("mv $pm->{dir}/$r.bam $output", 0);
	$pm->add_command("$SAMTOOLS index $output");
	$pm->del_file($aln1, $aln2, $fastq1, $fastq2) if($delete_input);
	$pm->check_filesize($output); # 1M
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_sampe",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "150:00:00"});

	return($qid);
}

sub samse {
	my $self = shift;
	
	my %param = ( "aln1" => undef,
	              "fastq1" => undef,
	              "genome" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $aln1 = to_abs_path($param{aln1});
	my $fastq1 = to_abs_path($param{fastq1});
	my $genome = to_abs_path($param{genome});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	my $r = time().rand();
	$pm->add_command("$BWA sampe $genome $aln1 $fastq1 | mbuffer -q -m 2G -l /dev/null | samtools view -uSbh - | mbuffer -q -m 2G -l /dev/null | samtools sort - $pm->{dir}/$r");
	$pm->add_command("mv $pm->{dir}/$r.bam $output", 0);
	$pm->add_command("$SAMTOOLS index $output");
	$pm->del_file($aln1, $fastq1) if($delete_input);
	$pm->check_filesize($output); # 1M
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_samse",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "150:00:00"});

	return($qid);
}


sub merge_sam {
	my $self = shift;
	
	my %param = ( "sam_list" => [],
	               "output" => undef,
				   "delete_input" => 0,
				   @_);
	
	my $sam_list = $param{sam_list};
	$sam_list = sapply($sam_list, \&to_abs_path);
	my $output = to_abs_path( $param{output} );
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	if(len($sam_list) == 1) {
		$pm->add_command("mv $sam_list->[0] $output", 0);
		$pm->add_command("if [ -f $sam_list->[0].bai ];then\n mv $sam_list->[0].bai $output.bai\nfi\n", 0);

		my $tmp = $sam_list->[0];
		$tmp =~s/\.[bs]am$/bai/;
		$pm->add_command("if [ -f $tmp ];then\n mv $tmp $output.bai\nfi\n", 0);

	} else {
		my $input_str;
		for(my $i = 0; $i < scalar(@$sam_list); $i ++) {
			$input_str .= "INPUT=$sam_list->[$i] ";
		}

		$pm->add_command("JAVA_OPTIONS=-Xmx16G picard.sh MergeSamFiles $input_str OUTPUT=$output SORT_ORDER=coordinate TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");
		for(my $i = 0; $i < scalar(@$sam_list); $i ++) {
			$pm->del_file($sam_list->[$i]) if($delete_input);
		}
	}
	
	$pm->check_filesize($output); # 1M
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_merge_sam",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "20GB",
										walltime => "100:00:00"});

	return($qid);
}


sub samtools_flagstat {
	my $self = shift;
	
	my %param = ( "sam" => undef,
				  "output" => undef,
				  @_);
	
	my $sam    = to_abs_path( $param{sam} );
	my $output = to_abs_path( $param{output} );
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("$SAMTOOLS flagstat $sam > $output");
	
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_flagstat",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "30:00:00"});
	return($qid);
}

sub picard_metrics {
	my $self = shift;
	
	my %param = ( "sam" => undef,
				"genome" => undef,
				  @_);
	
	my $sam    = to_abs_path( $param{sam} );
	my $genome    = to_abs_path( $param{genome} );
	my $output_aln = $sam.".aln.metrics";
	my $output_gcbias = $sam.".gcbias.metrics";
	my $output_insertsize = $sam.".insertsiz.metrics";
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("JAVA_OPTIONS=-Xmx10G picard.sh CollectAlignmentSummaryMetrics IS_BISULFITE_SEQUENCED=true INPUT=$sam OUTPUT=$output_aln REFERENCE_SEQUENCE=$genome TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");
	$pm->add_command("JAVA_OPTIONS=-Xmx10G picard.sh CollectGcBiasMetrics REFERENCE_SEQUENCE=$genome INPUT=$sam OUTPUT=$output_gcbias CHART_OUTPUT=$output_gcbias.pdf IS_BISULFITE_SEQUENCED=true TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");
	$pm->add_command("JAVA_OPTIONS=-Xmx10G picard.sh CollectInsertSizeMetrics HISTOGRAM_FILE=$output_insertsize.pdf INPUT=$sam OUTPUT=$output_insertsize REFERENCE_SEQUENCE=$genome TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_picard_metrics",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "20:00:00"});
	return($qid);
}


sub picard_insertsize {
	my $self = shift;
	
	my %param = ( "sam" => undef,
				  @_);
	
	my $sam    = to_abs_path( $param{sam} );
	my $output_insertsize = $sam.".insertsiz.metrics";
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("JAVA_OPTIONS=-Xmx10G picard.sh CollectInsertSizeMetrics HISTOGRAM_FILE=$output_insertsize.pdf INPUT=$sam OUTPUT=$output_insertsize TMP_DIR=$pm->{tmp_dir} VALIDATION_STRINGENCY=SILENT");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_picard_metrics",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "2GB",
										walltime => "20:00:00"});
	return($qid);
}


1;
