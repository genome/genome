package Genome::Model::Tools::DetectVariants2::MethRatio;

use warnings;
use strict;

use Genome;
use Workflow;
use File::Copy;
use Workflow::Simple;
use Cwd;

my $DEFAULT_VERSION = '2.6';

class Genome::Model::Tools::DetectVariants2::MethRatio {
    is => ['Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase'],
    doc => "Runs methyl counting script on a bsmap alignment model.",
    has => [
        chromosome_list => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of chromosomes to run on. This will run on a default set of chromosomes from the reference sequence if not set.',
        },
   ],
    has_transient_optional => [
        _snv_output_dir => {
            is => 'String',
            doc => 'The location of the snvs.hq file',
        },
    ],
};

sub _ensure_chromosome_list_set {
    my $self = shift;

    unless ($self->chromosome_list) {
        $self->chromosome_list($self->default_chromosomes);
    }
    return;
}

sub set_output {
    my $self = shift;

    unless ($self->snv_output) {
        $self->snv_output($self->_temp_staging_directory. '/snvs.hq');  #???
    }

    return;
}

sub add_bams_to_input {
    my ($self, $input) = @_;

    $input->{bam_file} = $self->aligned_reads_input;
    return;
}
sub get_reference {
    my $self = shift;

    my $refbuild_id = $self->reference_build_id;
    unless($refbuild_id){
        die $self->error_message("Received no reference build id.");
    }
    print "refbuild_id = ".$refbuild_id."\n";
    my $ref_seq_build = Genome::Model::Build->get($refbuild_id);
    return $ref_seq_build->full_consensus_path('fa');
}

sub workflow_xml {
    return \*DATA;
}

sub variant_type {
    return 'snv';
}

sub raw_output_file {
    my $self = shift;
    return $self->output_directory . "/" . $self->variant_type . "s.hq";
}

sub raw_inputs {
    my $self = shift;

    return map {$self->raw_input_for_chromosome($_)} @{$self->chromosome_list};
}

sub raw_input_for_chromosome {
    my ($self, $chromosome) = @_;
    return $self->output_directory . "/"
        .  Genome::Utility::Text::sanitize_string_for_filesystem($chromosome)
        . "/" . $self->variant_type . "s.hq";
}

sub _generate_standard_files {
    my $self = shift;
    my $staging_dir = $self->_temp_staging_directory;
    my $cat_raw = Genome::Model::Tools::Cat->create(dest => $self->raw_output_file,
        source => [$self->raw_inputs]);
    unless($cat_raw->execute){
        die $self->error_message("Cat command failed to execute.");
    }
    $self->SUPER::_generate_standard_files(@_);
    return 1;
}

sub _sort_detector_output {
    my $self = shift;
    return 1;
}

sub versions {
    return Genome::Model::Tools::Bsmap::MethRatioWorkflow->available_methratio_versions;
}

sub default_chromosome_list {
    my $self = shift;

    return $self->default_chromosomes;
}

sub default_chromosomes {
    my $self = shift;
    my $refbuild = Genome::Model::Build::ReferenceSequence->get($self->reference_build_id);
    die unless $refbuild;
    my $chromosome_array_ref = $refbuild->chromosome_array_ref;
    return $self->sort_chromosomes($refbuild->chromosome_array_ref);
}

sub default_chromosomes_as_string {
    return join(',', $_[0]->default_chromosomes);
}

sub chromosome_list_as_string {
    my $self = shift;
    my $chromosome_list = $self->chromosome_list;
    return join(',', @$chromosome_list);
}

sub params_for_detector_result {
    my $self = shift;
    my ($params) = $self->SUPER::params_for_detector_result;

    if ($self->chromosome_list) {
        $params->{chromosome_list} = $self->chromosome_list_as_string;
    } else {
        $params->{chromosome_list} = $self->default_chromosomes_as_string;
    }

    return $params;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="MethRatio Detect Variants Module">

  <link fromOperation="input connector" fromProperty="bam_file" toOperation="MethRatio" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="MethRatio" toProperty="output_directory" />
  <link fromOperation="input connector" fromProperty="chromosome_list" toOperation="MethRatio" toProperty="chromosome" />
  <link fromOperation="input connector" fromProperty="reference" toOperation="MethRatio" toProperty="reference" />
  <link fromOperation="input connector" fromProperty="version" toOperation="MethRatio" toProperty="version" />


  <link fromOperation="MethRatio" fromProperty="output_directory" toOperation="output connector" toProperty="output" />

  <operation name="MethRatio" parallelBy="chromosome">
    <operationtype commandClass="Genome::Model::Tools::Bsmap::MethRatioWorkflow" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="Y">bam_file</inputproperty>
    <inputproperty isOptional="Y">output_directory</inputproperty>
    <inputproperty isOptional="Y">chromosome_list</inputproperty>
    <inputproperty isOptional="Y">version</inputproperty>
    <inputproperty isOptional="Y">reference</inputproperty>

    <outputproperty>output</outputproperty>
  </operationtype>

</workflow>
