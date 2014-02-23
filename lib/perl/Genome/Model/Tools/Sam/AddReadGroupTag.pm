package Genome::Model::Tools::Sam::AddReadGroupTag;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::AddReadGroupTag {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        read_group_tag => {
            is  => 'String',
            doc => 'The value which will be added to every record of the input file in the RG tag and PG tag.',
        },
    ],
    has_optional => [
        input_file => {
            is  => 'String',
            doc => 'The SAM file to add a read group and program group tag to.',
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting file',
        },
        input_filehandle => {
            is => 'IO::File',
            doc => 'Pass in a filehandle to stream the SAM records in to add read groups',
        },
        output_filehandle => {
            is => 'IO::File',
            doc => 'Pass in a filehandle to stream read groups out of',
        },
        pass_sam_headers => {
            is => 'Boolean',
            doc => 'Should headers be passed through',
            default_value => 1
        }
    ],
};

sub help_brief {
    'Tool to add a read group tag to SAM files.';
}

sub help_detail {
    return <<EOS
    Tool to add a read group tag to SAM files.
EOS
}

sub execute {
    my $self = shift;

    my $input_fh;
    if ($self->input_file && -s $self->input_file) {
        $input_fh = Genome::Sys->open_file_for_reading($self->input_file);
    } elsif ($self->input_filehandle) {
        $input_fh = $self->input_filehandle;
    } else {
        $self->error_message("You must pass in either an input file or input filehandle.");
        return;
    }

    my $output_fh;
    if ($self->output_file) {
        if (-s $self->output_file )  {
           $self->error_message("The target file already exists at: " . $self->output_file . ". Please remove this file and rerun to generate a new merged file.");
           return;
        }

        $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    } elsif ($self->output_filehandle) {
        $output_fh = $self->output_filehandle;
    } else {
        $self->error_message("You must pass in either an output file or output filehandle.");
        return;
    }


    my $read_group_tag = $self->read_group_tag;

    $self->debug_message("Attempting to add read group tag $read_group_tag.");

    my $now = UR::Context->current->now;
    $self->debug_message(">>> Beginning add read tag at $now");

    my $rg_added_records = 0;

    while (my $line = $input_fh->getline) {
        my $first_char = substr($line, 0, 1);
        if ( $first_char ne '@') {
            $line =~ m/((?:\S*\s){11})(.*)$/;
            my $front = $1;
            my $back = $2;
            chomp $front;
            if (defined $back) {
                chomp $back;
            }
            if ($back eq "") {
                print $output_fh $front."\tRG:Z:$read_group_tag\tPG:Z:$read_group_tag\n";
            } else {
                print $output_fh $front."RG:Z:$read_group_tag\tPG:Z:$read_group_tag\t".$back."\n";
            }
            $rg_added_records++;
        } else {
            print $output_fh $line if ($self->pass_sam_headers);
        }
    }

    # don't close the file if we didn't open it!
    unless ($self->output_filehandle) {
        $output_fh->close;
    }
    $now = UR::Context->current->now;
    $self->debug_message("<<< Completed add read tag at $now.  Processed $rg_added_records SAM records.");

    return 1;
}


1;
