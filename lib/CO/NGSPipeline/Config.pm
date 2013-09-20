package CO::NGSPipeline::Config;

#########################################
# global configurations
#########################################

use strict;
use File::Spec;
use CO::NGSPipeline::Utils;
use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw($BWA
                 $CNYBWA
				 $PICARD_BIN_DIR
				 $EMAIL
                 );

our $BWA = "bwa-0.6.2-tpx";
our $CNYBWA = "cnybwa-0.6.2";
our $PICARD_BIN_DIR = "/ibios/tbi_cluster/11.4/x86_64/picard/picard-tools";
our $EMAIL = "z.gu\@dkfz.de";

our $TRIMPAIR_BIN = "/ibios/co02/guz/program/general/pairTrim/trim.pl";

1;