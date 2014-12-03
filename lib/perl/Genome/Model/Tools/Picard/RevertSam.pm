package Genome::Model::Tools::Picard::RevertSam;

use strict;
use warnings;

use Genome;
use Sort::strverscmp;

class Genome::Model::Tools::Picard::RevertSam {
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
        sort_order => {
            is => 'String',
            doc => 'The sort order to create the reverted output file with.',
            default_value => 'queryname',
            valid_values => ['unsorted', 'queryname', 'coordinate'],
            is_optional => 1,
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
        sanitize => {
            is => 'Boolean',
            doc => 'WARNING: This option is potentially destructive. If enabled will discard reads in order to produce a consistent output BAM. Reads discarded include (but are not limited to) paired reads with missing mates, duplicated records, records with mismatches in length of bases and qualities. This option can only be enabled if the output sort order is queryname and will always cause sorting to occur.',
            default_value => 0,
            is_optional => 1,
        },
        max_discard_fraction => {
            is => 'Float',
            doc => 'If SANITIZE=true and higher than MAX_DISCARD_FRACTION reads are discarded due to sanitization then the program will exit with an Exception instead of exiting cleanly. Output BAM will still be valid.',
            default_value => '0.01',
            is_optional => 1,
        },
        sample_alias => {
            is => 'String',
            doc => 'The sample alias to use in the reverted output file. This will override the existing sample alias in the file and is used only if all the read groups in the input file have the same sample alias',
            is_optional => 1,
        },
        library_name => {
            is => 'String',
            doc => 'The library name to use in the reverted output file. This will override the existing sample alias in the file and is used only if all the read groups in the input file have the same sample alias.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Reverts SAM or BAM files to a previous state by removing certain types of information and/or substituting in the original quality scores when available.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#BamIndexStats
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/RevertSam.jar';
    unless (-e $jar_path) {
        if ( strverscmp($self->use_version,'1.29') == -1 ) {
            die('Please define version 1.29 or greater.');
        } else {
            die('Missing jar path: '. $jar_path);
        }
    }
    my $cmd = $jar_path .' net.sf.picard.sam.RevertSam INPUT='. $self->input_file .' OUTPUT='. $self->output_file .' SORT_ORDER='. $self->sort_order;

    my @boolean_attributes = qw/
                                   restore_original_qualities
                                   remove_duplicate_information
                                   remove_alignment_information
                               /;

    # Sanitize option was first made available in Picard v1.108
    if ( strverscmp($self->use_version,'1.107') == 1 ) {
        push @boolean_attributes, 'sanitize';
        $cmd .= ' MAX_DISCARD_FRACTION='. $self->max_discard_fraction;
    } elsif ($self->sanitize) {
        die('Version '. $self->use_version .' does not support the sanitize option.');
    }

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
