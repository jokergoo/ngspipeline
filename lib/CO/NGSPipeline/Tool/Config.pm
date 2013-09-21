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
				 $PICARD_BIN_DIR
				 $TRIMPAIR_BIN_DIR
                 );

our $BWA = "bwa-0.6.2-tpx";
our $CNYBWA = "cnybwa-0.6.2";
our $PICARD_BIN_DIR = "/ibios/tbi_cluster/11.4/x86_64/picard/picard-tools";

our $TRIMPAIR_BIN_DIR = "/ibios/co02/guz/program/general/pairTrim";

1;