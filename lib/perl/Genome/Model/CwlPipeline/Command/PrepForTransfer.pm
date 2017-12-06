package Genome::Model::CwlPipeline::Command::PrepForTransfer;

use strict;
use warnings;

use Genome;

class Genome::Model::CwlPipeline::Command::PrepForTransfer {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build::CwlPipeline',
            is_many => 1,
            doc => 'The build(s) to prepare for data transfer',
        },
        directory => {
            is => 'Text',
            doc => 'The diretory to prepare and write symlinks.',
        },
        md5sum => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Create a md5sum in the MANIFEST for each file.',
            default_value => 0,
        },
    ],
    doc => 'prepare a directory for transfer of CwlPipeline builds',
};

sub help_detail {
    return <<EOHELP
Prepare a directory for transfer of CwlPipeline build(s) files. Include a MANIFEST of files, subject names and optionally md5 checksums.
EOHELP
    ;
}

sub execute {
    my $self = shift;

    my @headers = qw/
                        subject.name
                        file.name
                    /;
    if ($self->md5sum) {
        push @headers, 'md5sum';
    }
    if (-d $self->directory) {
        $self->fatal_message('Unable to use existing directory: '. $self->directory);
    } else {
        Genome::Sys->create_directory($self->directory);
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->directory .'/MANIFEST',
        separator => "\t",
        headers => \@headers,
    );
    for my $build ($self->builds) {
        my $results_directory = $build->data_directory .'/results';
        my @file_paths = grep {-f $_} glob($results_directory .'/*');
        for my $file_path (@file_paths) {
            my ($file_name, $dir, $suffix) = File::Basename::fileparse($file_path);
            my $symlink_name = $build->id .'.'. $file_name;
            my $symlink_path = $self->directory .'/'. $symlink_name;
            my %data = (
                'subject.name' => $build->model->subject->name,
                'file.name' => $symlink_name,
            );
            if ($self->md5sum) {
                my $md5sum = Genome::Sys->md5sum($file_path);
                $data{md5sum} = $md5sum;
            }
            $writer->write_one(\%data);
            Genome::Sys->create_symlink($file_path,$symlink_path);
        }
    }
    $writer->output->close();

    return 1;
}


1;
