package Genome::InstrumentData::Command::MergeAlignments;

use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use Genome;

use Genome::Utility::Text;

class Genome::InstrumentData::Command::MergeAlignments {
    is => ['Command::V2'],
    has_input => [
        alignment_result_ids => {
            is => 'ARRAY',
            doc => 'The alignments to be merged',
        },
        merger_name => {
            is => 'Text',
            doc => 'The name of the merge program to use (e.g. "samtools")',
        },
        merger_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the merge program',
        },
        merger_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the merge program to use',
        },
        duplication_handler_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'The name of the program to use for marking or removing duplicate reads',
        },
        duplication_handler_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the dpulication handler',
        },
        duplication_handler_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the duplication handler to use',
        },
        refiner_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'Name of the refiner to use',
        },
        refiner_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the refiner to use',
        },
        refiner_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Params for the refiner to use',
        },
        refiner_known_sites_ids => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'ID of the variant list to use for refinement',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of Samtools to use when needed by mergers/deduplicators',
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

sub params_for_merged_alignment {
    my $self = shift;

    my @alignment_results = Genome::InstrumentData::AlignmentResult->get($self->alignment_result_ids);

    my $filters = [];
    my $segments = [];
    for my $result (@alignment_results) {
        if($result->filter_name) {
            push @$filters, join(':', $result->instrument_data_id, $result->filter_name);
        }

        if($result->instrument_data_segment_id) {
            push @$segments, join(':', $result->instrument_data_id, $result->instrument_data_segment_id, $result->instrument_data_segment_type);
        }
    }
    my @instrument_data_ids = uniq (map($_->instrument_data_id, @alignment_results));
    my %params = (
        #after merged alignments are refactored can replace most params with this:
        #alignment_result_id => [$self->alignment_result_ids],
        instrument_data_id => \@instrument_data_ids,
        reference_build_id => $alignment_results[0]->reference_build_id,
        annotation_build_id => $alignment_results[0]->annotation_build_id,
        aligner_name => $alignment_results[0]->aligner_name,
        aligner_version => $alignment_results[0]->aligner_version,
        aligner_params => $alignment_results[0]->aligner_params,
        picard_version => $alignment_results[0]->picard_version,
        force_fragment => $alignment_results[0]->force_fragment,
        trimmer_name => $alignment_results[0]->trimmer_name,
        trimmer_version => $alignment_results[0]->trimmer_version,
        trimmer_params => $alignment_results[0]->trimmer_params,

        merger_name => $self->merger_name,
        merger_params => $self->merger_params,
        merger_version => $self->merger_version,
        duplication_handler_name => $self->duplication_handler_name,
        duplication_handler_params => $self->duplication_handler_params,
        duplication_handler_version => $self->duplication_handler_version,
        samtools_version => $self->samtools_version || undef,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );

    if(scalar @$filters) {
        $params{filter_name} = $filters;
    }
    if(scalar @$segments) {
        $params{instrument_data_segment} = $segments;
    }

    return \%params;
}


sub shortcut {
    my $self = shift;

    #try to get using the lock in order to wait here in shortcut if another process is creating this alignment result
    my $result = $self->_process_merged_alignment('get_with_lock');
    unless($result) {
        $self->debug_message('No existing alignment found.');
        return;
    }

    $self->debug_message('Using existing alignment ' . $result->__display_name__);

    my $refiner_result = $self->_process_refinement('shortcut', $result);
    unless($refiner_result) {
        $self->debug_message('No existing refinement found.');
        return;
    }

    return ref($refiner_result)? $refiner_result : $result;
}

sub execute {
    my $self = shift;

    my $result = $self->_process_merged_alignment('get_or_create');
    unless($result) {
        $self->error_message('Failed to generate merged alignment.');
        die $self->error_message;
    }
    $self->debug_message('Generated merged alignment');

    my $refiner_result = $self->_process_refinement('execute', $result);
    unless($refiner_result) {
        $self->error_message('Failed to generate refinement.');
        die $self->error_message;
    }

    return ref($refiner_result)? $refiner_result : $result;
}

sub _process_merged_alignment {
    my $self = shift;
    my $mode = shift;

    my $params = $self->params_for_merged_alignment();
    unless($params) {
        $self->error_message('Could not get parameters for merger');
        return;
    }
    my $result = eval { Genome::InstrumentData::AlignmentResult::Merged->$mode(%$params) };
    if($@) {
        $self->error_message($mode . ': ' . $@);
        return;
    }

    $self->result_id($result->id) if $result;
    return $result;
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
    return $result;
}

sub _refiner_for_name {
    my $self = shift;
    my $name = shift;

    $name =~ s/-/_/g;
    return 'Genome::InstrumentData::Command::RefineReads::' . Genome::Utility::Text::string_to_camel_case($name);
}


1;
 
