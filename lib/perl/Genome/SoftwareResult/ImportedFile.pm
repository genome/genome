package Genome::SoftwareResult::ImportedFile;

use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::ImportedFile {
    is => 'Genome::SoftwareResult::StageableSimple::SingleFile',
    has_input => [
        file_content_hash => {
            is => 'Text',
            doc => 'MD5 hash of the file that is contained in this result',
        },
    ],
    has_optional_metric => [
        source_file_path => {
            is => 'Path',
            doc => 'Path to the original file used to create this result',
        }
    ],
};

sub create {
    my ($class, %params) = @_;

    my $file = delete $params{source_file_path};
    $class->fatal_message('No source file path given to create software result!') if not defined $file;
    $file = Cwd::abs_path($file);
    $class->fatal_message('Source file path (%s) given to create software result does not exist!', $file) if not -s $file;

    my $md5 = Genome::Sys->md5sum($file);
    $class->fatal_message('No md5 for source file! %s', $file) if not $md5;

    return $class->SUPER::create(
        source_file_path => $file,
        file_content_hash => $md5,
    );
}

sub _run{
    my $self = shift;
    Genome::Sys->copy_file($self->source_file_path, $self->_temp_staging_file_path);
}

sub _file_name {
    my $self = shift;
    return File::Basename::basename($self->source_file_path);
}

1;
