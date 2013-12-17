package Genome::Model::Tools::DetectVariants2::Result::Combine::LqUnion;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Combine::LqUnion{
    is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
    has_input => [
        result_ids => {
            is => 'Text',
            is_many => 1,
            doc => 'ids of the results whose lq results to union',
        },
    ],
    has_param => [
        variant_type => {
            is => 'Text',
            valid_values => ['snv', 'indel'],
        },
    ],
    has => [
        results => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            calculate => q{ return Genome::Model::Tools::DetectVariants2::Result::Base->get([$self->result_ids]); },
            is_many => 1,
        },
    ],
    doc => 'Unions the lq files of a set of other results',
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless ($self);

    unless($self->_validate_inputs) {
        die $self->error_message('Failed to validate inputs.');
    }

    unless($self->_prepare_staging_directory) {
        die $self->error_message('Failed to prepare staging directory.');
    }

    unless($self->_combine_variants){
        die $self->error_message('Failed to combine variants');
    }

    unless($self->_validate_output) {
        die $self->error_message('Failed to validate output.');
    }

    unless ($self->_prepare_output_directory) {
        die $self->error_message('Failed to prepare output directory.');
    }

    unless($self->_promote_data) {
        die $self->error_message('Failed to promote data.');
    }

    unless($self->_reallocate_disk_allocation) {
        #this is suboptimal, but no need to die as our data should still be good and complete
        $self->warning_message('Failed to reallocate disk allocation.');
    }

    unless($self->_add_as_user_of_inputs) {
        die $self->error_message('Failed to add self as user of inputs.');
    }

    return $self;
}

sub _validate_inputs {
    my $self = shift;

    my @r = $self->results;
    my @r_ids = $self->result_ids;

    return (scalar(@r_ids) == scalar(@r));
}

sub compare_results {
    my ($first, $second) = @_;

    ($first->class cmp $second->class) or
    ($first->can('input_result_a') && $second->can('input_result_a')? &compare_results($first->input_result_a, $second->result_a) : 0) or
    ($first->can('input_result_b') && $second->can('input_result_b')? &compare_results($first->input_result_b, $second->result_b) : 0) or
    ($first->can('detector_name')? $first->detector_name : '') cmp ($second->can('detector_name')? $second->detector_name : '') or
    ($first->can('detector_version')? $first->detector_version : '') cmp ($second->can('detector_version')? $second->detector_version : '') or
    ($first->can('detector_params')? $first->detector_params : '') cmp ($second->can('detector_params')? $second->detector_params : '') or
    ($first->can('filter_name')? $first->filter_name : '') cmp ($second->can('filter_name')? $second->filter_name : '') or
    ($first->can('filter_version')? $first->filter_version : '') cmp ($second->can('filter_version')? $second->filter_version : '') or
    ($first->can('filter_params')? $first->filter_params : '') cmp ($second->can('filter_params')? $second->filter_params : '') or
    $first->id cmp $second->id #cry
}

sub _combine_variants {
    my $self = shift;
    my @results = $self->results;

    @results = sort { &compare_results($a, $b); } @results;

    my $type = $self->variant_type;
    my $file = $self->_file_for_type($type);
    my @files = map($_->path($file), @results);

    @files = grep(-e $_, @files);

    my $result_file = join('/', $self->temp_staging_directory, $self->_file_for_type($type));
    if(scalar(@files) == 0) {
        Genome::Sys->shellcmd(
            cmd => 'touch ' . $result_file,
            output_files => [$result_file],
            allow_zero_size_output_files => 1,
        );
    } elsif(scalar(@files) == 1) {
        Genome::Sys->create_symlink($files[0], $result_file);
    } else {
        #TODO This approach doesn't really work for cnv or sv... but to date we don't have multiple LQ files to combine.
        my $sort_command = Genome::Model::Tools::Joinx::Sort->create(
            input_files => \@files,
            merge_only => 1,
            output_file => $result_file,
        );

        unless ($sort_command->execute) {
            $self->error_message("Error executing lq sort command");
            die $self->error_message;
        }
    }

    return 1;
}

sub _validate_output {
    my $self = shift;

    #TODO wc -l all the files and make sure everything's there if we're paranoid
    return 1;
}

sub _add_as_user_of_inputs {
    my $self = shift;

    for my $r ($self->results) {
        $r->add_user(label => 'uses', user => $self);
    }

    return 1;
}

sub _file_for_type {
    my $self = shift;
    my $type = shift;

    my $file = $type . 's.lq';
    if($type eq 'snv' or $type eq 'indel') {
        $file .= '.bed';
    }

    return $file;
    return;
}

sub _needs_symlinks_followed_when_syncing { 1 };
sub _working_dir_prefix { 'dv2-lq-union-result' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    return join('/', 'build_merged_alignments', $self->id, 'dv2-lq-union-' . $staged_basename);
};

1;
