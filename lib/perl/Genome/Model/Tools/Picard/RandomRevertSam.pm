package Genome::Model::Tools::Picard::RandomRevertSam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::RandomRevertSam {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'A BAM or SAM file to revert.',
        },
        output_file => {
            is => 'String',
            doc => 'The reverted BAM or SAM file.',
        },
        probability => {
            is => 'Number',
            doc => 'The probability of keeping an alignment, between 0 and 1.  Default value: 0.0. This option can be set to \'null\' to clear the default value.',
            is_optional => 1,
        },
        random_seed => {
            is => 'Integer',
            doc => 'Random seed to use if reproducibilty is desired.  Default value: 1. This option can be set to \'null\' to clear the default value.',
            is_optional => 1,
        },
        sort_order => {
            is => 'String',
            is_optional => 1,
            doc => 'The sort order to create the reverted output file with.',
            default_value => 'queryname',
            valid_values => ['unsorted', 'queryname', 'coordinate'],
        },
        restore_original_qualities => {
            is => 'Boolean',
            doc => 'True to restore original qualities from the OQ field to the QUAL field if available.',
            default_value => 1,
            is_optional => 1,
        },
        remove_duplicate_information => {
            is => 'Boolean',
            doc => 'Remove duplicate read flags from all reads.',
            default_value => 1,
            is_optional => 1,
        },
        remove_alignment_information => {
            is => 'Boolean',
            doc => 'Remove all alignment information from the file.',
            default_value => 1,
            is_optional => 1,
        },
        attribute_to_clear => {
            is => 'String',
            is_many => 1,
            doc => 'When removing alignment information, the set of optional tags to remove. This option may be specified 0 or more times.',
            is_optional => 1,
        },
        sample_alias => {
            is => 'String',
            doc => 'The sample alias to use in the reverted output file. This will override the existing sample alias in the file and is used only if all the read groups in the input file have the same sample alias',
            is_optional => 1,
        },
        library_name => {
            is => 'String',
            is_optional => 1,
            doc => 'The library name to use in the reverted output file. This will override the existing sample alias in the file and is used only if all the read groups in the input file have the same sample alias.'
        },
    ],
};

sub help_brief {
    'Randomly Reverts SAM or BAM files to a previous state by removing certain types of information and/or substituting in the original quality scores when available.';
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;
    
    my $picard_version = $self->use_version;
    my $picard_dir = $self->picard_path;
    my $picard_jar_path = $picard_dir .'/sam-'. $picard_version .'.jar';
    my $sam_jar_path = $picard_dir .'/picard-'. $picard_version .'.jar';
    my $jar_path = $self->class->base_dir . '/RandomRevertSam.jar';
    my $cmd = $jar_path .' net.sf.picard.sam.RandomRevertSam INPUT='. $self->input_file .' OUTPUT='. $self->output_file .' SORT_ORDER='. $self->sort_order;
    my @boolean_attributes = qw/
                                    restore_original_qualities
                                    remove_duplicate_information
                                    remove_alignment_information
                                /;
    my $string = $self->resolve_boolean_attributes_string(\@boolean_attributes);
    $cmd .= $string;
    if (defined($self->sample_alias)) {
        $cmd .= ' SAMPLE_ALIAS='. $self->sample_alias;
    }
    if (defined($self->library_name)) {
        $cmd .= ' LIBRARY_NAME='. $self->library_name;
    }
    for my $attribute ($self->attribute_to_clear) {
        $cmd .= ' ATTRIBUTE_TO_CLEAR='. $attribute;
    }
    if (defined($self->probability)) {
        $cmd .= ' PROBABILITY='. $self->probability;
    }
    if (defined($self->random_seed)) {
        $cmd .= ' RANDOM_SEED='. $self->random_seed;
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

sub resolve_boolean_attributes_string {
    my $self = shift;
    my $attributes = shift;

    my $cmd = '';
    for my $attribute (@{$attributes}) {
        unless (defined($self->$attribute)) { next; }
        my $flag;
        if ($self->$attribute == 1) {
            $flag = 'true';
        } elsif ($self->$attribute == 0) {
            $flag = 'false';
        } else {
            $self->error_message('Failed to parse boolean attribute '. $attribute .' with value '. $self->$attribute);
            die($self->error_message);
        }
        my $uc_attribute = uc($attribute);
        $cmd .= ' '. $uc_attribute .'='. $flag;
    }
    return $cmd;
}

1;
