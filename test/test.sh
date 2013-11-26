perl /home/guz/perllib/ngs/WGBS_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_wgbs --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/whole_genome_bisulfite_sequencing/analysis_bismark --tool bismark
perl /home/guz/perllib/ngs/WGBS_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_wgbs --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/whole_genome_bisulfite_sequencing/analysis_bsmap --tool bsmap
perl /home/guz/perllib/ngs/WGBS_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_wgbs --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/whole_genome_bisulfite_sequencing/analysis_methylctools --tool methylctools


perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_gf --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_defuse --tool defuse
perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_gf --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_fusionmap --tool fusionmap
perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_gf --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_fusionhunter --tool fusionhunter


perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_rnaseq --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_star --tool star
perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_rnaseq --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_star_nodup --tool star --nodup
perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_rnaseq --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_tophat --tool tophat
perl /home/guz/perllib/ngs/RNASeq_pipeline.pl --list /home/guz/perllib/ngs/list/list_test_rnaseq --dir /icgc/lsdf/mb/analysis/guz/pipeline_test/rna_sequencing/analysis_tophat_nodup --tool tophat --nodup
