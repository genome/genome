package Genome::Model::Tools::DetectVariants2::CopyCatSomaticWithBamWindow;

use warnings;
use strict;

use Cwd;
use Genome;
use Workflow::Simple;

my $DEFAULT_VERSION = '0.1';

class Genome::Model::Tools::DetectVariants2::CopyCatSomaticWithBamWindow{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces somatic copy-number calls from paired samples",
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[maxtmp>10000] span[hosts=1]',
        },
    ],
    has => [
        bamwindow_filter_to_chromosomes => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'Chromosomes that bamwindow will filter output to (generally for test purposes)',
        },
        #TODO: annotation_version is deprecated.  Remove it
        annotation_version => {
            is => 'Text',
            doc => 'version of the copycat annotation to use.  Deprecated',
            is_optional => 1,
        }
    ],
};


sub _detect_variants {
    my $self = shift;
    my $cnvs = $self->_temp_staging_directory."/cnvs.hq";

    ##parse input params string - expected format
    #--bamwindow-version 0.4 --bamwindow-params [-w 10000 -r -l -s]
    my $params = $self->params;
    my $bamwindow_version;
    my $bamwindow_params;

    if ($params =~ m/--bamwindow-version/) {
        $params =~ m/--bamwindow-version\s*(\d+\.?\d?)\s*/;
        $bamwindow_version = $1;
        $params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
        $bamwindow_params  = $1;
        $bamwindow_params =~ s/--bamwindow-params\s*//;
    }

    my $copycat_params;
    #--copycat-params [--per-read-length --per-library]
    if ($params =~ m/--copycat-params/) {
        $params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
        $bamwindow_params = $1;
        $params =~ m/--copycat-params\s*\{([^\}]+)\}/;
        $copycat_params = $1;
    }
    my $per_library = 0;
    my $per_read_length = 0;
    if($copycat_params =~ m/--per-library/) {
        $per_library = 1;
    }
    if($copycat_params =~ m/--per-read-length/) {
        $per_read_length = 1;
    }

    #annotation-version
    my $annotation_version;
    if($params =~ m/--annotation-version/){
        $params =~ m/--annotation-version\s*(\d+\.?\d?)\s*/;
        $annotation_version = $1;
    }

    my $samtools_strategy;
    if($params =~ m/--samtools-strategy/) {
        $params =~ m/--samtools-strategy\s*\{([^\}]+)\}/;
        $samtools_strategy = $1;
    }

    unless($samtools_strategy) {
        die $self->error_message("No --samtools-strategy found in params");
    }

    my %input;

    # Define a workflow from the static XML at the bottom of this module
    my $workflow = Workflow::Operation->create_from_xml(\*DATA);

    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    # Collect and set input parameters
    #for bam-window
    $input{bamwindow_version} = $bamwindow_version;
    $input{bamwindow_params} = $bamwindow_params;
    $input{tumor_bam_file} = $self->aligned_reads_input;
    $input{normal_bam_file} = $self->control_aligned_reads_input;
    $input{tumor_window_file} = join('/', $self->_temp_staging_directory, 'tumor_window_file');
    $input{normal_window_file} = join('/', $self->_temp_staging_directory, 'normal_window_file');
    my @bamwindow_filter_to_chromosomes = $self->bamwindow_filter_to_chromosomes;
    $input{bamwindow_filter_to_chromosomes} = \@bamwindow_filter_to_chromosomes;

    #for copycat
    $input{per_read_length} = $per_read_length;
    $input{per_library} = $per_library;
    $input{tumor_samtools_file} = $self->get_samtools_results(
        $self->aligned_reads_input,
        $self->aligned_reads_sample,
        $samtools_strategy,
        $self->output_directory . '/tumor-samtools-result',
    );
    $input{normal_samtools_file} = $self->get_samtools_results(
        $self->control_aligned_reads_input,
        $self->control_aligned_reads_sample,
        $samtools_strategy,
        $self->output_directory . '/normal-samtools-result',
    );
    $input{copycat_output_directory} = $self->_temp_staging_directory;

    my $annotation_sr = $self->_find_annotation_data(
        $self->reference_build,
        $annotation_version,
        $self->result_users
    );

    $input{annotation_data_id} = $annotation_sr->id;

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }
    $workflow->log_dir($log_dir);

    # Launch workflow
    $self->debug_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %input);

    # Collect and analyze results
    unless($result){
        if (@Workflow::Simple::ERROR){
            print Data::Dumper->Dumper(@Workflow::Simple::ERROR), "\n";
        }
        die $self->error_message("Workflow did not return correctly");
    }

    unless(-e $cnvs){
        system("touch $cnvs");
    }

    return 1;
}

sub _find_annotation_data {
    my $self = shift;
    my $reference_build = shift;
    my $annotation_version = shift;
    my $result_users = shift;

    my $annotation_reference = $reference_build;
    my $annotation_sr;
    until ($annotation_sr) {
        unless ($annotation_reference) {
            die $self->error_message(
                'No annotation data found for version %s and reference %s',
                $annotation_version,
                $reference_build->id,
            );
        }

        $annotation_sr = Genome::Model::Tools::CopyCat::AnnotationData->get_with_lock(
            reference_sequence => $annotation_reference,
            version            => $annotation_version,
            users              => $result_users,
        );
    } continue {
        $annotation_reference = $annotation_reference->derived_from;
    }

    if($annotation_sr->reference_sequence ne $reference_build) {
        $self->status_message(
            'No annotation data for reference %s.  Using data for reference %s from which it was derived.',
            $reference_build->id,
            $annotation_sr->reference_sequence->id,
        );
    }

    return $annotation_sr;
}

sub has_version {
    return 1; #FIXME implement this when this module is filled out
}

sub _sort_detector_output {
    return 1;
}

sub get_samtools_results{
    my $self = shift;
    my $bam_path = shift;
    my $sample_name = shift;
    my $samtools_strategy = shift;
    my $output_directory = shift;

    Genome::Sys->create_directory($output_directory) unless -d $output_directory;

    my $dispatcher = Genome::Model::Tools::DetectVariants2::Dispatcher->create(
        snv_detection_strategy => $samtools_strategy,
        aligned_reads_input => $bam_path,
        aligned_reads_sample => $sample_name,
        reference_build_id => $self->reference_build->id,
        output_directory => $output_directory,
        result_users => $self->result_users,
    );
    $dispatcher->execute or die $self->error_message('Failed to run dispatcher to get samtools result');

    my $result = $dispatcher->snv_result;
    unless($result) {
        die $self->error_message('Failed to find result after running dispatcher');
    }

    my $snv_hq = $result->path('snvs.hq');
    unless($snv_hq) {
        die $self->error_message('no snv file found on result');
    }

    $result->add_user(user => $self, label => 'uses_for_running_copycat');
    return $snv_hq;
}

sub _promote_staged_data {
    my $self = shift;
    my $output_directory = $self->SUPER::_promote_staged_data;
    Genome::Sys->remove_directory_tree($self->_temp_staging_directory);
    return $output_directory;
}

sub _create_temp_directories {
    my $self = shift;
    my $staging_tempdir = File::Temp->newdir(
            "staging-XXXXX",
            DIR     => $self->output_directory,
            CLEANUP => 0,
    );
    $self->_temp_staging_directory($staging_tempdir->dirname);
    $self->_temp_scratch_directory(Genome::Sys->create_temp_directory);
    return 1;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>
    <workflow name="CopyCatSomatic Detect Variants Module">


    <link fromOperation="input connector" fromProperty="normal_bam_file" toOperation="BamWindow Normal" toProperty="bam_file"/>
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Normal" toProperty="options" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Normal" toProperty="version" />
    <link fromOperation="input connector" fromProperty="normal_window_file" toOperation="BamWindow Normal" toProperty="output_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_filter_to_chromosomes" toOperation="BamWindow Normal" toProperty="filter_to_chromosomes" />


    <link fromOperation="input connector" fromProperty="tumor_bam_file" toOperation="BamWindow Tumor" toProperty="bam_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Tumor" toProperty="options" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Tumor" toProperty="version" />
    <link fromOperation="input connector" fromProperty="tumor_window_file" toOperation="BamWindow Tumor" toProperty="output_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_filter_to_chromosomes" toOperation="BamWindow Tumor" toProperty="filter_to_chromosomes" />


    <link fromOperation="BamWindow Normal" fromProperty="output_file" toOperation="CopyCat Somatic" toProperty="normal_window_file" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_file" toOperation="CopyCat Somatic" toProperty="tumor_window_file" />
    <link fromOperation="input connector" fromProperty="per_library" toOperation="CopyCat Somatic" toProperty="per_library" />
    <link fromOperation="input connector" fromProperty="per_read_length" toOperation="CopyCat Somatic" toProperty="per_read_length" />
    <link fromOperation="input connector" fromProperty="tumor_samtools_file" toOperation="CopyCat Somatic" toProperty="tumor_samtools_file" />
    <link fromOperation="input connector" fromProperty="normal_samtools_file" toOperation="CopyCat Somatic" toProperty="normal_samtools_file" />
    <link fromOperation="input connector" fromProperty="copycat_output_directory" toOperation="CopyCat Somatic" toProperty="output_directory" />
    <link fromOperation="input connector" fromProperty="annotation_data_id" toOperation="CopyCat Somatic" toProperty="annotation_data_id" />


    <link fromOperation="BamWindow Normal" fromProperty="output_file" toOperation="output connector" toProperty="bam_window_normal_output_file" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_file" toOperation="output connector" toProperty="bam_window_tumor_output_file" />
    <link fromOperation="CopyCat Somatic" fromProperty="output_directory" toOperation="output connector" toProperty="copycat_output_directory" />

    <operation name="BamWindow Normal">
    <operationtype commandClass="Genome::Model::Tools::BamWindow" typeClass="Workflow::OperationType::Command" />
    </operation>

    <operation name="BamWindow Tumor">
    <operationtype commandClass="Genome::Model::Tools::BamWindow" typeClass="Workflow::OperationType::Command" />
    </operation>

    <operation name="CopyCat Somatic">
    <operationtype commandClass="Genome::Model::Tools::CopyCat::Somatic" typeClass="Workflow::OperationType::Command" />
    </operation>


    <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>tumor_bam_file</inputproperty>
    <inputproperty>normal_bam_file</inputproperty>
    <inputproperty>tumor_window_file</inputproperty>
    <inputproperty>normal_window_file</inputproperty>
    <inputproperty>bamwindow_params</inputproperty>
    <inputproperty>bamwindow_version</inputproperty>
    <inputproperty>per_read_length</inputproperty>
    <inputproperty>per_library</inputproperty>
    <inputproperty>tumor_samtools_file</inputproperty>
    <inputproperty>normal_samtools_file</inputproperty>
    <inputproperty>copycat_output_directory</inputproperty>
    <inputproperty>annotation_data_id</inputproperty>
    <inputproperty isOptional="Y">bamwindow_filter_to_chromosomes</inputproperty>

    <outputproperty>bam_window_normal_output_file</outputproperty>
    <outputproperty>bam_window_tumor_output_file</outputproperty>
    <outputproperty>copycat_output_directory</outputproperty>
    </operationtype>

    </workflow>
