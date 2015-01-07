package Genome::FeatureList::Command::Import;

use strict;
use warnings;

use feature qw(say);

use Genome;

class Genome::FeatureList::Command::Import {
    is => 'Command::V2',
    has_input => [
        file_path => {
            is => 'Text',
            doc => 'Path to a file or directory to import',
            shell_args_position => 1,
        },
        name => {
            is => 'Text',
            doc => 'A name for the new feature-list',
        },
        reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference for the coordinates in the supplied BED file',
        },
        content_type => {
            is => 'Text',
            default => 'targeted',
            valid_values => Genome::FeatureList->__meta__->property('content_type')->valid_values,
        },
    ],
    has_optional_input => [
        is_1_based => {
            is => 'Boolean',
            default => 0,
            doc => 'Indicates this is a "1-based" BED file instead of a true-BED (0-based) file',
        },
        source => {
            is => 'Text',
            doc => 'Where the BED file came from (e.g. Agilent, Nimblegen)',
        },
        description => {
            is => 'Text',
            doc => 'General description of the BED file',
        },
    ],
    has_transient_optional => [
        is_multitracked => {
            is => 'Boolean',
            default => 0,
            doc => 'determined during the import process',
        },
    ],
    has_optional_output => [
        new_feature_list => {
            is => 'Genome::FeatureList',
            doc => 'the newly imported feature list',
        },
    ],
    doc => 'command to import a new feature-list from a BED file',
};

sub help_brief {
    return 'Import a feature-list from a capture BED file';
}

sub help_synopsis {
    return <<EOS
genome feature-list import --name="New capture probe set" --reference=NCBI-human-build36 /path/to/capture.bed
EOS
;
}

sub help_detail {
    return <<EOHELP
This command provides an interface around `genome feature-list create` that performs additional BED file
sanitization and validation to ensure the resulting FeatureList is suitable for use as a "region of interest"
set in analysis pipelines.
EOHELP
;
}

sub execute {
    my $self = shift;

    my $file_path = $self->file_path;
    unless(-e $file_path) {
        die $self->error_message('File does not exist: %s', $file_path);
    }

    if(-d $file_path) {
        $file_path = $self->_resolve_bed_file_from_directory($file_path);
    }

    my $sanitized_bed_path = $self->validate_and_sanitize_bed($file_path);
    unless($sanitized_bed_path) {
        die $self->error_message('Unable to import BED file due to validation errors.');
    }

    my $md5 = Genome::Sys->md5sum($sanitized_bed_path);
    my $imported_feature_list = $self->find_existing_list($md5, $self->name);

    unless($imported_feature_list) {
        my $create_cmd = Genome::FeatureList::Command::Create->create(
            file_path => $sanitized_bed_path,
            reference => $self->reference,
            content_type => $self->content_type,
            format => Genome::FeatureList->_derive_format($self->is_1_based, $self->is_multitracked),
            description => ($self->description || 'created with `genome feature-list import`'),
            source => $self->source,
            name => $self->name,
        );
        $imported_feature_list = $create_cmd->execute;
        unless($imported_feature_list) {
            die $self->error_message('Failed to execute command to create FeatureList');
        }
    }

    say $imported_feature_list->id;
    $self->new_feature_list($imported_feature_list);

    return 1;
}

sub validate_bed_content_line {
    my $self = shift;
    my $line = shift;

    unless(exists $self->{_chromosome_hash}) {
        my $reference = $self->reference;

        my $chromosome_list = $reference->chromosome_array_ref;
        for my $chr (@$chromosome_list) {
            $self->{_chromosome_hash}{$chr} = 1;
        }
    }
    my $chromosomes = $self->{_chromosome_hash};

    my @fields = Genome::FeatureList->_parse_entry($line);
    unless(exists $chromosomes->{$fields[0]}) {
        die $self->error_message(
            'Chromosome %s not found in reference %s',
            $fields[0],
            $self->reference->name
        );
    }

    return 1;
}

sub validate_and_sanitize_bed {
    my $self = shift;
    my $input_bed_file = shift;
    my $input_bed_fh = Genome::Sys->open_file_for_reading($input_bed_file);

    my ($sanitized_bed_fh, $sanitized_bed_file) = Genome::Sys->create_temp_file;

    while(my $line = <$input_bed_fh>) {
        next if $line =~ /^browser/;
        $line =~ s/\015$//; #remove CRs if BED came from Windows

        if($line =~ /^track/) {
            $self->is_multitracked(1);
            my %track_attrs = Genome::FeatureList->parse_track_definition($line);
            my $name = Genome::FeatureList->_standardize_track_name($track_attrs{name});
            unless (grep { $_ eq $name } ('probes', 'targets')) {
                die $self->error_message(
                    'Unknown track name %s found in %s, line %s.',
                    $name,
                    $input_bed_file,
                    $input_bed_fh->input_line_number,
                );
            }
            $line = qq(track name="$name"\n);
        } else {
            $self->validate_bed_content_line($line);
        }

        $sanitized_bed_fh->print($line);
    }

    $input_bed_fh->close;
    $sanitized_bed_fh->close;

    return $sanitized_bed_file;
}

sub _resolve_bed_file_from_directory {
    my $self = shift;
    my $directory = shift;

    unless(-d $directory) {
        die $self->error_message('Could not read %s as a directory.', $directory);
    }

    my @bed_files = glob(join '/', $directory, '*.bed');

    unless(scalar(@bed_files)) {
        die $self->error_message('No BED files found in directory %s.', $directory);
    }

    if(scalar(@bed_files) == 1) {
        return $bed_files[0];
    }

    my @all_tracks = grep { $_ =~ '/[^/]+_AllTracks.bed$' } @bed_files;
    if(scalar(@all_tracks) == 1) {
        return $all_tracks[0];
    }

    die $self->error_message('Multiple candidate BED files found in directory %s. Please select one.', $directory);
}

sub find_existing_list {
    my $class = shift;
    my ($bed_md5, $name) = @_;

    my @existing_lists = Genome::FeatureList->get(file_content_hash => $bed_md5);

    unless(@existing_lists) {
        my $list_with_same_name = Genome::FeatureList->get(name => $name);
        if($list_with_same_name) {
            die $class->error_message(
                'Found existing feature-list with the same name but a different BED file: %s',
                $list_with_same_name->id,
            );
        } else {
            return;
        }
    }

    my $list_to_use;

    if(@existing_lists == 1) {
        $list_to_use = $existing_lists[0];
    } else {
        my @matches = grep { $_->name eq $name } @existing_lists;
        if(@matches == 1) {
            $list_to_use = $matches[0];
        }
    }

    unless($list_to_use) {
        $class->error_message(
            'Multiple matching existing FeatureLists found: %s',
            join(' ', map $_->id, @existing_lists)
        );
        die $class->error_message('To proceed make sure the --name parameter matches one of the existing lists.');
    }

    $class->status_message('Found existing list: %s', $list_to_use->__display_name__);
    return $list_to_use;
}

1;
