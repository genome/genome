package Genome::Model::Tools::Sam::Compare;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;
use Sys::Hostname;
use Genome::Utility::AsyncFileSystem qw(on_each_line);

class Genome::Model::Tools::Sam::Compare {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        file1 => {
            is  => 'String',
            doc => 'The first file in the comparison',
        },
        file2 => {
            is  => 'String',
            doc => 'The second file in the comparison',
        },
    ],
};

sub help_brief {
    'Tool to compare BAM or SAM files';
}

sub help_detail {
    return <<EOS
    Tool to compare BAM or SAM files.
EOS
}

sub execute {
    my $self = shift;
    
    my $output_file = Genome::Sys->create_temp_file_path();
    
    my $compare_cmd = Genome::Model::Tools::Picard::CompareSams->create(
        input_file_1 => $self->file1,
        input_file_2 => $self->file2,
        output_file  => $output_file,
        
        validation_stringency => 'SILENT',
    );
    
    unless($compare_cmd->execute()) {
        $self->error_message('Failed to execute command');
        die($self->error_message);
    }

    my $ret = 0;
    my $response;

    my $output_fh = Genome::Sys->open_file_for_reading($output_file);
    while (<$output_fh>) {
        if (m/Cannot compare/) {
            $ret = 0;
            last;
        }
        if (m/Differ\s*0$/) {
            $ret = 1;
        }
        $response .= $_;
    }
    close($output_fh);

    print $response;

    return $ret;

}

1;
