package Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase;

use strict;
use warnings;
use Genome;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(chromosome_sort);

class Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase {
    doc => "A base class for detectors that run workflows",
    is_abstract => 1,
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has => [
        chromosome_list => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of chromosomes to run on. This will run on a default set of chromosomes from the reference sequence if not set.',
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => $ENV{GENOME_LSF_QUEUE_DV2_WORKER},
        },
    ],
    has_transient_optional => [
        _workflow_result => {
            doc => 'Result of the workflow',
        },
    ],
};

sub sort_chromosomes {
    my $self = shift;
    my $chromosome_array_ref = shift;
    my @chromosome_list = sort {chromosome_sort()} @$chromosome_array_ref;
    return \@chromosome_list;
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
    elsif ($a_chrom eq 'MT') {
        $rv = -1;
    }
    elsif ($b_chrom eq 'MT') {
        $rv = 1;
    }
    else {
        $rv = ($a_chrom cmp $b_chrom);
    }

    return $rv;
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

sub _promote_staged_data {
    my $self = shift;
    return 1;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }

    for my $v ($self->versions){
        if($v eq $version){
            return 1;
        }
    }

    return 0;
}

sub _create_temp_directories {
    my $self = shift;
    $self->_temp_staging_directory($self->output_directory);
    $self->_temp_scratch_directory($self->output_directory);
    return 1;

    return $self->SUPER::_create_temp_directories(@_);
}

sub _detect_variants {
    my $self = shift;

    $self->_ensure_chromosome_list_set;

    $self->set_output;

    # Define a workflow from the static XML at the bottom of this module
    my $workflow = Workflow::Operation->create_from_xml($self->workflow_xml);

    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    my %input;
    $input{chromosome_list} = $self->chromosome_list;
    $input{reference} = $self->get_reference;
    $input{output_directory} = $self->output_directory;
    $input{version} = $self->version;

    $self->add_bams_to_input(\%input);

    $self->_dump_workflow($workflow);

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }
    $workflow->log_dir($log_dir);

    Genome::Sys->disconnect_default_handles;

    # Launch workflow
    $self->debug_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %input);

    # Collect and analyze results
    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }
    $self->_workflow_result($result);

    return 1;
}

sub default_chromosomes {
    my $self = shift;
    my $refbuild = Genome::Model::Build::ReferenceSequence->get($self->reference_build_id);
    die unless $refbuild;
    my $chromosome_array_ref = $refbuild->chromosome_array_ref;
    return $self->sort_chromosomes($refbuild->chromosome_array_ref);
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

sub _ensure_chromosome_list_set {
    my $self = shift;

    unless ($self->chromosome_list) {
        $self->chromosome_list($self->default_chromosomes);
    }
    return;
}

1;

