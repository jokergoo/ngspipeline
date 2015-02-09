package CO::NGSPipeline::Tool::RNAseq::Common;

##############################################################################
# provide command method for each pipeline. It is a base module for all specific
# pipelines.

use strict;
use CO::NGSPipeline::Tool::RNAseq::Config;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;
use File::Basename;

use base qw/CO::NGSPipeline/;

sub rnaseqqc {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	              "sample_id" => undef,
	              @_);
	
	my $bam    = to_abs_path( $param{bam} );
	my $sample_id = $param{sample_id};
	
	my $pm = $self->get_pipeline_maker;
	
	my $bam_base = basename($bam);
	my $bam_rg = "$bam_base.RG.bam";
	my $bam_reorder = "$bam_base.reorder.bam";
	
	open OUT, ">$pm->{dir}/sample_list";
	print OUT "Sample ID	Bam File	Notes\n";
	print OUT "$sample_id\t$bam_reorder\tNo Note\n";
	close OUT;
	
	$pm->add_command("picard.sh AddOrReplaceReadGroups INPUT=$bam OUTPUT=$bam_rg RGID=readGroup_name RGLB=readGroup_name RGPL=illumina RGPU=run RGSM=sample_name SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=$pm->{tmp_dir}");
	$pm->add_command("picard.sh ReorderSam INPUT=$bam_rg OUTPUT=$bam_reorder CREATE_INDEX=true R=$GENOME_HG19 TMP_DIR=$pm->{tmp_dir}");
	$pm->del_file("$bam_rg", "$bam_base.RG.bai");
	$pm->add_command("java -jar /home/guz/soft/GenePatternServer/taskLib/RNASeQC.2.0/RNAseqMetrics.jar -s $pm->{dir}/sample_list -t $GENCODE_GTF -r $GENOME_HG19 -n 1000 -o $pm->{dir}/rnaseqqc");
	$pm->del_file("$bam_reorder", "$bam_base.reorder.bai");
	
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_rnaseqqc",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "20GB",
										walltime => "40:00:00"});
	return($qid);

}

sub rpkm {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	              "strand" => 0,
	              @_);
	
	my $bam    = to_abs_path( $param{bam} );
	my $strand = $param{"strand"};
	
	my $bam_prefix = $bam;
	$bam_prefix =~s/\.bam$//;
	
	my $stand_opt = $strand ? "-S" : "";
	
	my $pm = $self->get_pipeline_maker;
	
	open OUT, ">$pm->{tmp_dir}/exonCoverage.sh";
	print OUT <<SH;
samtools view  -bu -q 1 $bam | coverageBed $stand_opt -split -abam stdin -b /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/RefSeq/RefSeq_Nov15_2011_from_annovar_Exons_plain.bed.gz > $bam_prefix.Exons_RPKM.bed.tmp

mv $bam_prefix.Exons_RPKM.bed.tmp $bam_prefix.Exons_RPKM.bed
total=\$(awk '{a+=\$(NF-3);}END{print a}' $bam_prefix.Exons_RPKM.bed)
echo -e "#chrom\tchromStart\tchromEnd\tname\tstrand\texonNr\tlength\treads\tbases_covered\tlength\tcoverage\tRPKM" > $bam_prefix.Exons_RPKM.bed.tmp
perl /icgc/ngs_share/ngsPipelines/tools/RPKM.pl $bam_prefix.Exons_RPKM.bed \$total | sort -k1,1d -k2,2n >> $bam_prefix.Exons_RPKM.bed.tmp
mv $bam_prefix.Exons_RPKM.bed.tmp $bam_prefix.Exons_RPKM.bed
bgzip -c -d /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/RefSeq/RefSeq_Nov15_2011_from_annovar_Genes_plain.bed.gz | perl /icgc/ngs_share/ngsPipelines/ExonCoverage/countcoverageBedRPKM.pl $bam_prefix.Exons_RPKM.bed - \$total | awk 'NR==1; NR > 1 {print \$0 | "sort -k1,1d -k2,2n"}' > $bam_prefix.Genes_RPKM.bed.tmp
mv $bam_prefix.Genes_RPKM.bed.tmp $bam_prefix.Genes_RPKM.bed

SH
	close OUT;
	
	$pm->add_command("sh $pm->{tmp_dir}/exonCoverage.sh");
	
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_RPKM",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "10:00:00"});
	return($qid);
}

sub counting {
	my $self = shift;
	
	my %param = ( "bam" => undef,
		          "sample_id" => undef,
	              @_);
	
	my $bam    = to_abs_path( $param{bam} );
	my $sample_id = $param{sample_id};
	
	my $pm = $self->get_pipeline_maker;

	$pm->add_command("
S1=`awk 'NR==8{print \$5}' $pm->{dir}/rnaseqqc/$sample_id/$sample_id.metrics.txt`
S2=`awk 'NR==8{print \$6}' $pm->{dir}/rnaseqqc/$sample_id/$sample_id.metrics.txt`

if [ `echo \"\$S1<30\" | bc` -eq 1 ]
then
	STRAND='reverse'
elif [ `echo \"\$S1>70\" | bc` -eq 1 ]
then
	STRAND='yes'
else
	STRAND='no'
fi

echo \"strand=\$STRAND\"
", 0);
	
	my $dexseq_script_dir = "/ibios/tbi_cluster/11.4/x86_64/R/R-3.0.0/lib64/R/library/DEXSeq/python_scripts";
 	
	$pm->add_command("
$SAMTOOLS view -F 0x0404 $bam | htseq-count -s \$STRAND -t exon -m intersection-nonempty - $GENCODE_GTF > $pm->{dir}/$sample_id.mkdup.exon.nodup.count
$SAMTOOLS view -F 0x0004 $bam | htseq-count -s \$STRAND -t exon -m intersection-nonempty - $GENCODE_GTF > $pm->{dir}/$sample_id.mkdup.exon.count

if [ \"\$STRAND\" != 'no' ]
then
	if [ \"\$STRAND\" == 'yes' ]
	then
		ANTISENSE='reverse'
	else
		ANTISENSE='yes'
	fi
	
	$SAMTOOLS view -F 0x0404 $bam | htseq-count -s \$ANTISENSE -t exon -m intersection-nonempty - $GENCODE_GTF > $pm->{dir}/$sample_id.mkdup.exon.nodup.strand.count
	$SAMTOOLS view -F 0x0004 $bam | htseq-count -s \$ANTISENSE -t exon -m intersection-nonempty - $GENCODE_GTF > $pm->{dir}/$sample_id.mkdup.exon.strand.count
fi

# $SAMTOOLS view -F 0x0404 $bam | python $dexseq_script_dir/dexseq_count.py -p yes -r name -s \$STRAND $GENCODE_DEXSEQ_GTF - $pm->{dir}/$sample_id.mkdup.dexseq.nodup.count
# $SAMTOOLS view -F 0x0004 $bam | python $dexseq_script_dir/dexseq_count.py -p yes -r name -s \$STRAND $GENCODE_DEXSEQ_GTF - $pm->{dir}/$sample_id.mkdup.dexseq.count

# if [ \"\$STRAND\" != 'no' ]
# then
#	if [ \"\$STRAND\" == 'yes' ]
#	then
#		ANTISENSE='reverse'
#	else
#		ANTISENSE='yes'
#	fi
#	
#	$SAMTOOLS view -F 0x0404 $bam | python $dexseq_script_dir/dexseq_count.py -p yes -r name -s \$ANTISENSE $GENCODE_DEXSEQ_GTF - $pm->{dir}/$sample_id.mkdup.dexseq.nodup.strand.count
#	$SAMTOOLS view -F 0x0004 $bam | python $dexseq_script_dir/dexseq_count.py -p yes -r name -s \$ANTISENSE $GENCODE_DEXSEQ_GTF - $pm->{dir}/$sample_id.mkdup.dexseq.strand.count
# fi

", 0);
 	
 	#$pm->del_file("$bam");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_counting",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "1GB",
										walltime => "100:00:00"},);
	return($qid);
}

1;
