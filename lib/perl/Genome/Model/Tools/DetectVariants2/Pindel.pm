package Genome::Model::Tools::DetectVariants2::Pindel;

use warnings;
use strict;

use Genome;
use Workflow;
use File::Copy;
use Workflow::Simple;
use Cwd;
use Genome::Utility::Text;

my $DEFAULT_VERSION = '0.2';
my $PINDEL_COMMAND = 'pindel_64';

class Genome::Model::Tools::DetectVariants2::Pindel {
    is => ['Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase'],
    doc => "Runs the pindel pipeline on the last complete build of a somatic model.",
    has => [
        chromosome_list => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of chromosomes to run on.',
        },
        chr_mem_usage => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of mem to request per chromosomes to run on.',
        },
   ],
    has_transient_optional => [
        _indel_output_dir => {
            is => 'String',
            doc => 'The location of the indels.hq.bed file',
        },
        _chr_mem_usage => {
            doc => 'This is a hashref containing the amount of memory in MB to request for each chromosome job of pindel',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => "-M 25000000 -R 'select[mem>25000] rusage[mem=25000]'",
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

    unless ($self->indel_bed_output) {
        $self->indel_bed_output($self->_temp_staging_directory. '/indels.hq.bed');
    }

    return;
}

sub add_bams_to_input {
    my ($self, $input) = @_;

    $input->{tumor_bam} = $self->aligned_reads_input;
    $input->{normal_bam} = $self->control_aligned_reads_input if defined $self->control_aligned_reads_input;
    return;
}

sub get_reference {
    my $self = shift;

    my $refbuild_id = $self->reference_build_id;
    unless($refbuild_id){
        die $self->error_message("Received no reference build id.");
    }
    print "refbuild_id = ".$refbuild_id."\n";

    return $refbuild_id;
}

sub workflow_xml {
    return \*DATA;
}

sub variant_type {
    return 'indel';
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

sub versions {
    return Genome::Model::Tools::Pindel::RunPindel->available_pindel_versions;
}

sub default_chromosomes_as_string {
    return join(',', $_[0]->default_chromosomes);
}

sub params_for_detector_result {
    my $self = shift;
    my ($params) = $self->SUPER::params_for_detector_result;
    my $chrom_string = $self->default_chromosomes_as_string;
    if(length($chrom_string) >= Genome::SoftwareResult::Param->__meta__->property(property_name => 'value_id')->data_length) {
        #FIXME if the default is ever changed this breaks
        $chrom_string = 'all';

    }
    $params->{chromosome_list} = $chrom_string; #$self->default_chromosomes_as_string;

    return $params;
}

sub _sort_detector_output {
    my $self = shift;
    return 1;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Pindel Detect Variants Module">

  <link fromOperation="input connector" fromProperty="normal_bam" toOperation="Pindel" toProperty="control_aligned_reads_input" />
  <link fromOperation="input connector" fromProperty="tumor_bam" toOperation="Pindel" toProperty="aligned_reads_input" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="Pindel" toProperty="output_directory" />
  <link fromOperation="input connector" fromProperty="chromosome_list" toOperation="Pindel" toProperty="chromosome" />
  <link fromOperation="input connector" fromProperty="version" toOperation="Pindel" toProperty="version" />
  <link fromOperation="input connector" fromProperty="reference" toOperation="Pindel" toProperty="reference_build_id" />

  <link fromOperation="Pindel" fromProperty="output_directory" toOperation="output connector" toProperty="output" />

  <operation name="Pindel" parallelBy="chromosome">
    <operationtype commandClass="Genome::Model::Tools::Pindel::RunPindel" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="Y">normal_bam</inputproperty>
    <inputproperty isOptional="Y">tumor_bam</inputproperty>
    <inputproperty isOptional="Y">output_directory</inputproperty>
    <inputproperty isOptional="Y">version</inputproperty>
    <inputproperty isOptional="Y">chromosome_list</inputproperty>
    <inputproperty isOptional="Y">reference</inputproperty>

    <outputproperty>output</outputproperty>
  </operationtype>

</workflow>
