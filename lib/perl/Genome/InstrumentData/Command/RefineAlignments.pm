package Genome::InstrumentData::Command::RefineAlignments;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::Command::RefineAlignments {
    is => ['Command::V2'],
    has_input => [
        merge_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            doc => 'The merge result to be merged',
        },
        refiner_name => {
            is => 'Text',
            doc => 'Name of the refiner to use',
        },
        refiner_version => {
            is => 'Text',
            doc => 'Version of the refiner to use',
        },
        refiner_params => {
            is => 'Text',
            doc => 'Params for the refiner to use',
        },
        refiner_known_sites_ids => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'ID of the variant list to use for refinement',
        },
    ],
    has_optional_output => [
        result_id => {
            is => 'Text',
            doc => 'The result generated/found when running the command',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => &bsub_rusage,
        },
    ],
};

sub bsub_rusage {
    return "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[tmp=90000:mem=16000]' -M 16000000";
}

# sub shortcut {
    # my $self = shift;

    # #try to get using the lock in order to wait here in shortcut if another process is creating this alignment result
    # my $result = $self->_process_merged_alignment('get_with_lock');
    # unless($result) {
        # $self->debug_message('No existing alignment found.');
        # return;
    # }

    # $self->debug_message('Using existing alignment ' . $result->__display_name__);

    # my $refiner_result = $self->_process_refinement('shortcut', $result);
    # unless($refiner_result) {
        # $self->debug_message('No existing refinement found.');
        # return;
    # }

    # return ref($refiner_result)? $refiner_result : $result;
# }

sub execute {
    my $self = shift;

    my $refiner_result = $self->_process_refinement('execute', $self->merge_result);
    unless($refiner_result) {
        $self->error_message('Failed to generate refinement.');
        die $self->error_message;
    }

    return $refiner_result;
}

sub _process_refinement {
    my $self = shift;
    my $mode = shift;
    my $merged_result = shift;

    unless($self->refiner_name) {
        return 1;
    }

    my @known_sites_ids = $self->refiner_known_sites_ids;
    my @known_sites = Genome::Model::Build::ImportedVariationList->get(id => \@known_sites_ids);
    my %params = (
        version => $self->refiner_version,
        params => $self->refiner_params,
        known_sites => \@known_sites,
        bam_source => $merged_result
    );

    my $cmd_class_name = $self->_refiner_for_name($self->refiner_name);
    my $cmd = $cmd_class_name->create(%params);
    if ( not $cmd ) {
        $self->error_message("Failed to create refiner command $cmd_class_name with params ".Data::Dumper::Dumper(\%params));
        return;
    }
    my $result = eval { $cmd->$mode; };
    if($@) {
        $self->error_message($mode . ': ' . $@);
        return;
    }

    $self->result_id($result->id) if $result;
    $DB::single=1;
    return $result;
}

sub _refiner_for_name {
    my $self = shift;
    my $name = shift;

    $name =~ s/-/_/g;
    return 'Genome::InstrumentData::Command::RefineReads::' . Genome::Utility::Text::string_to_camel_case($name);
}

1;
