package CO::NGSPipeline;

use strict;

sub set_pipeline_maker {
	my $pipeline = shift;
	my $pm = shift;
	
	$pipeline->{_pipeline_maker} = $pm;
	
	return $pipeline;
}

sub get_pipeline_maker {
	my $pipeline = shift;
	
	return $pipeline->{_pipeline_maker};
}



############################################################################
# tools section
############################################################################
use CO::NGSPipeline::Tool::BSseq::BSMAP;
use CO::NGSPipeline::Tool::BSseq::Bismark;
use CO::NGSPipeline::Tool::BSseq::Bsmooth;
use CO::NGSPipeline::Tool::BSseq::BisSNP;
use CO::NGSPipeline::Tool::BSseq::methylCtools;
use CO::NGSPipeline::Tool::RNAseq::GeneFusion;
use CO::NGSPipeline::Tool::RNAseq::TopHat;
use CO::NGSPipeline::Tool::RNAseq::STAR;
use CO::NGSPipeline::Tool::RNAseq::GSNAP;

sub bsmap {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::BSseq::BSMAP->new->set_pipeline_maker($pm);
}

sub bismark {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::BSseq::Bismark->new->set_pipeline_maker($pm);
}

sub bsmooth {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::BSseq::Bsmooth->new->set_pipeline_maker($pm);
}

sub bissnp {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::BSseq::BisSNP->new->set_pipeline_maker($pm);
}

sub methylctools {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::BSseq::methylCtools->new->set_pipeline_maker($pm);
}


sub tophat {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::RNAseq::TopHat->new->set_pipeline_maker($pm);
}

sub star {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::RNAseq::STAR->new->set_pipeline_maker($pm);
}

sub gsnap {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::RNAseq::GSNAP->new->set_pipeline_maker($pm);
}

sub genefusion {
	my $self = shift;
	my $pm = $self->get_pipeline_maker;
	
	return CO::NGSPipeline::Tool::RNAseq::GeneFusion->new->set_pipeline_maker($pm);
}


1;
