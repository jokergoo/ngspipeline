# Light weight NGS pipeline maker.

A pipeline is in fact a list of commands. More correctly, a pipeline combines a
list of independent steps or jobs in which each jobs are composed of closely
related commands (such as `samtools view`, `samtools sort` and `samtools index`
are always embed in one step). With knowing the dependency between jobs, the 
complete dependency tree can be built.

For PBS system, a pipeline maker should at least do following things:

- put commands into shell scripts.
- set PBS parameters such as memory, time requested, number of CPUs.
- solve the dependency between jobs.

A more stronger pipeline maker can support:

- monitor status of jobs (i.e. whether it is failed or successful).
- check whether jobs are successful finished.
- recovery from last failed jobs instead of start from beginning of the pipeline.

Here I developed a light weight pipeline maker. It allows quick implementing
a NGS pipeline. When it is fully tested, it can be immigrated to more stronger
pipeline administrators such as `Roddy`.

In this repository, there is a `CO::PipelineMaker` class which do the pipeline
maker stuff and `CO:::NGSPipeline` classes which implement tools and pipelines
for WGBS pipeline and RNAseq pipeline.

## `CO::PipelineMaker`

Mainly do following things:

- write commands into shell script.
- add job checking in shell script.
- add job flags
- handle job dependency
- submit by qsub

To use this class, first create a `CO::PipelineMaker` object which must contains
the path of the working directory for the pipeline.

    my $pm = CO::PipelineMaker->new("dir" => $dir);
  
Set the job name and dependency

    $pm->set_job_name("test_job");
    $pm->set_job_dependency($previous_pid1, $previous_pid2);
  
Add commands to this step

    $pm->add_command($cmd1);
    $pm->add_command($cmd2);
  
Check file size of some output files. This is important since in some NGS commands,
error occurs while the commands only finish without return a non-zeor exit code.

    $pm->check_filesize($output);

Finally, when all settings for current job are done, you can prepare the shell
script, configure PBS settings and submit the job. `run` will return the PID of
this job and you can use it to set next job's dependency.

    my $pid = $pm->run("-l" => { walltime => '10:00:00',
	                             memory   => '10G',
					             nodes    => '1:ppn=8' }
			          );

Normally, for a NGS step, I will embed above codes into a single subroutine. 
Examples can be found under `CO::NGSPipeline::Tool` classes.

## `CO::NGSPipeline` namespace

There are two classes in this namespace. `CO::NGSPipeline::Tool` defines
independent jobs and `CO::NGSPipeline::Pipeline` defines pipelines by integrating
methods from `CO::NGSPipeline::Tool`.

### `CO::NGSPipeline::Tool`

Under each level of different namespace, there is a `Common` class which defines
common methods which is base class of more specific classes. Also, there is a 
`Config` module which defines global variables.

#### `CO::NGSPipeline::Tool::Common`
- fastqc
- trim
- sort sam/bam
- samtools view
- merge and remove duplicates
- flagstat
- bwa aln
- bwa sampe
- merge sam/bam
- picard metric
- picard insertsize
	
#### `CO::NGSPipeline::Tool::Config`
- global variables for common methods

Currently I implemented methods for WGBS pipelines and RNAseq pipelines.

### `CO::NGSPipeline::Tool::BSseq`

Methods for WGBS pipeline
	 
#### `CO::NGSPipeline::Tool::BSseq::Common`
Common methods for WGBS pipeline, including:

- QC
- save methylation into R `Bsseq` object.
	 
#### `CO::NGSPipeline::Tool::BSseq::Config`
Global variables for WGBS pipeline
	 
#### `CO::NGSPipeline::Tool::BSseq::Bismark`
Methods for Bismark pipeline, including:

- alignment
- Bismark's methylation calling
- lambda conversion rate
	 
#### `CO::NGSPipeline::Tool::BSseq::BSMAP`
Methods for BSMAP pipeline, including:

- alignment
- BSMAP's methylation calling
- lambda conversion rate
	 
#### `CO::NGSPipeline::Tool::BSseq::methylCtools`
Methods for methylCtools pipeline, including:

- fqconv
- bconv
- bcall
- lambda conversion rate
	 
#### `CO::NGSPipeline::Tool::BSseq::BisSNP`

- steps for BisSNP methylation calling

#### `CO::NGSPipeline::Tool::BSseq::Bsmooth`

- save methylation as RData for downstream DMR calling

### `CO::NGSPipeline::Tool::RNAseq`
Methods for RNAseq pipeline
	 
#### `CO::NGSPipeline::Tool::RNAseq::Common`

- rnaseqqc
- rpkm
- counting
	 
#### `CO::NGSPipeline::Tool::RNAseq::Config`
Global variables for RNA seq pipeline
	 
#### `CO::NGSPipeline::Tool::RNAseq::GeneFusion`
Methods for gene fusion pipeline

- defuse
- fusionmap
- fusionhunter
- tophatfusion
	 
#### `CO::NGSPipeline::Tool::RNAseq::GSNAP`
Methods for GSNAP pipeline

- alignment
	 
#### `CO::NGSPipeline::Tool::RNAseq::STAR`
Methods for STAR pipeline

- alignment
	 
#### `CO::NGSPipeline::Tool::RNAseq::TopHat`
Methods for TopHat pipeline

- alignment
 
### `CO::NGSPipeline::Pipeline`
Integrated pipelines which is a collection of methods from `CO::NGSPipeline::Tool`

### `CO::NGSPipeline`
This class provides 'shortcut' methods to call real methods in `CO::NGSPipeline::Tool`
namespace. For example, we already had a method called `align` in `CO::NGSPipeline::Tool::BSseq::BSMAP`.
In order to use this methed in a pipeline, you do not need to initialize the object
and deal with PipelineMaker stuff. Just using the 'shortcut' method:

    $pipeline->bsmap->align(@param);
	
in which `$pipeline` is a pipeline object and should be initialized with
a pipeline maker object. `$pipeline->bsmap` is a shortcut method which will initialize 
a `CO::NGSPipeline::Tool::BSseq::BSMAP` object and attach the pipeline maker object, 
finally you can call `align` method on this BSMAP object.
  
### `CO::NGSPipeline::Report`
Scirpts for pipeline report, currently only reports for WGBS pipeline.
  
### `CO::NGSPipeline::Getopt`
All established pipelines need paired-end FastQ files, so we need pathes of FastQ
files and the sample names. Also, working directory as well as some running mode
are common for all pipelines. Therefore, this class take charge of command line
parameters, construct and print help messages, validate parameters and finally
returns validated and transformed variables.
  
## `CO::Utils`
Util methods for all `CO` classes
