package Genome::ModelGroup::Command::RunPindel;

use warnings;
use strict;

use Genome;
use Workflow;
use File::Copy;
use Workflow::Simple;
use Cwd;

my $DEFAULT_VERSION = '0.5';
my $PINDEL_COMMAND = 'pindel_64';

class Genome::ModelGroup::Command::RunPindel {
    is => ['Genome::Command::Base'],
    doc => "Runs the pindel pipeline on the last complete build of a somatic model.",
    has => [
        model_group_id=> {
            is => 'String',
            is_optional=>0,
            doc=>"which model group to process",
        },
        output_directory=> {
            is=>'String',
            is_optional=>0,
            doc=>"directory to put results",
        },
        chromosome_list => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of chromosomes to run on.',
            default_value => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'],
        },
        chr_mem_usage => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of mem to request per chromosomes to run on.',
        },

   ],
    has_transient_optional => [
        _reference_build_id=> {
        },
        _workflow_result => {
            doc => 'Result of the workflow',
        },
        _indel_output_dir => {
            is => 'String',
            doc => 'The location of the indels.hq.bed file',
        },
        _chr_mem_usage => {
            doc => 'This is a hashref containing the amount of memory in MB to request for each chromosome job of pindel',
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW},
        },
    ],
};

sub execute {
    my $self = shift;
    # Obtain normal and tumor bams and check them. Either from somatic model id or from direct specification. 

    # Set default params
    my $mg = Genome::ModelGroup->get($self->model_group_id);
    my ($model, undef) = $mg->models;
    my $refbuild_id = $model->reference_sequence_build->id;
    $self->_reference_build_id($refbuild_id);
    unless($refbuild_id){
        die $self->error_message("Received no reference build id.");
    }
    print "refbuild_id = ".$refbuild_id."\n";

    my %input;

    # Define a workflow from the static XML at the bottom of this module
    my $workflow = Workflow::Operation->create_from_xml(\*DATA);
    
    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    my @chrom_list = $self->default_chromosomes;

    # Collect and set input parameters
    $input{chromosome_list} = \@chrom_list;
    $input{reference_build_id} = $refbuild_id;
    $input{model_group_id} = $self->model_group_id;
    $input{output_directory}  =  $self->output_directory;#$self->_temp_staging_directory;
    
    $self->_dump_workflow($workflow);
    $workflow->log_dir($self->output_directory);

    # Launch workflow
    $self->status_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %input);

    # Collect and analyze results
    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }
    $self->_workflow_result($result);

    return 1;
}

sub _dump_workflow {
    my $self = shift;
    my $workflow = shift;
    my $xml = $workflow->save_to_xml;
    my $xml_location = $self->output_directory."/workflow.xml";
    my $xml_file = Genome::Sys->open_file_for_writing($xml_location);
    print $xml_file $xml;
    $xml_file->close;
    #$workflow->as_png($self->output_directory."/workflow.png"); #currently commented out because blades do not all have the "dot" library to use graphviz
}

sub _create_temp_directories {
    my $self = shift;
    $self->_temp_staging_directory($self->output_directory);
    $self->_temp_scratch_directory($self->output_directory);
    return 1;

    return $self->SUPER::_create_temp_directories(@_);
}


sub _generate_standard_files {
    my $self = shift;
    my $staging_dir = $self->_temp_staging_directory;
    my $output_dir  = $self->output_directory;
    my @chrom_list = $self->default_chromosomes;
    my $test_chrom = $chrom_list[0];
    my $raw_output_file = $output_dir."/indels.hq";
    my @raw_inputs = map { $output_dir."/".$_."/indels.hq" } @chrom_list;
    my $cat_raw = Genome::Model::Tools::Cat->create( dest => $raw_output_file, source => \@raw_inputs);
    unless($cat_raw->execute){
        die $self->error_message("Cat command failed to execute.");
    }
    $self->SUPER::_generate_standard_files(@_);
    return 1;
}

sub _promote_staged_data {
    my $self = shift;
    return 1;
}

sub _sort_detector_output {
    my $self = shift;
    return 1;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my @versions = Genome::Model::Tools::Pindel::RunPindel->available_pindel_versions;

    for my $v (@versions){
        if($v eq $version){
            return 1;
        }
    }

    return 0;
}

sub chromosome_sort {
    # numeric sort if chrom starts with number
    # alphabetic sort if chrom starts with non-number
    # numeric before alphabetic
    # Not perfect, NT_1 will sort the same as NT_100
    my ($a_chrom) = $a =~ /^(\S+)/;
    my ($b_chrom) = $b =~ /^(\S+)/;
    ($a_chrom, $b_chrom) = (uc($a_chrom), uc($b_chrom));
    my $a_is_numeric = ($a_chrom =~ /^\d+$/ ? 1 : 0);
    my $b_is_numeric = ($b_chrom =~ /^\d+$/ ? 1 : 0);

    my $rv;
    if ($a_chrom eq $b_chrom) {
        $rv = 0;
    }
    elsif ($a_is_numeric && $b_is_numeric) {
        $rv = ($a_chrom <=> $b_chrom);
    }
    elsif ($a_is_numeric && !$b_is_numeric) {
        $rv = -1;
    }
    elsif ($b_is_numeric && !$a_is_numeric) {
        $rv = 1;
    }
    elsif ($a_chrom eq 'X') {
        $rv = -1;
    }
    elsif ($b_chrom eq 'X') {
        $rv = 1;
    }
    elsif ($a_chrom eq 'Y') {
        $rv = -1;
    }
    elsif ($b_chrom eq 'Y') {
        $rv = 1;
    }
    else {
        $rv = ($a_chrom cmp $b_chrom);
    }

    return $rv;
}

sub default_chromosomes {
    my $self = shift;
    my $refbuild = Genome::Model::Build::ReferenceSequence->get($self->_reference_build_id);
    die unless $refbuild;
    my $chromosome_array_ref = $refbuild->chromosome_array_ref;
    my @chromosome_list = sort chromosome_sort @$chromosome_array_ref;
    return @chromosome_list;
}

sub default_chromosomes_as_string {
    return join(',', $_[0]->default_chromosomes);
}

sub params_for_detector_result {
    my $self = shift;
    my ($params) = $self->SUPER::params_for_detector_result;

    $params->{chromosome_list} = $self->default_chromosomes_as_string;

    return $params;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Pindel Detect Variants Module">

  <link fromOperation="input connector" fromProperty="model_group_id" toOperation="Pindel" toProperty="model_group_id" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="Pindel" toProperty="output_directory" />
  <link fromOperation="input connector" fromProperty="chromosome_list" toOperation="Pindel" toProperty="chromosome" />
  <link fromOperation="input connector" fromProperty="reference_build_id" toOperation="Pindel" toProperty="reference_build_id" />

  <link fromOperation="Pindel" fromProperty="output_directory" toOperation="output connector" toProperty="output" />

  <operation name="Pindel" parallelBy="chromosome">
    <operationtype commandClass="Genome::Model::Tools::Pindel::RunPindelModelGroup" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="N">model_group_id</inputproperty>
    <inputproperty isOptional="Y">output_directory</inputproperty>
    <inputproperty isOptional="Y">chromosome_list</inputproperty>
    <inputproperty isOptional="Y">reference_build_id</inputproperty>

    <outputproperty>output</outputproperty>
  </operationtype>

</workflow>
