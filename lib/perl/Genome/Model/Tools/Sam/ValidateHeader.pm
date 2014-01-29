package Genome::Model::Tools::Sam::ValidateHeader;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::ValidateHeader {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        input_file => {
            is  => 'String',
            doc => 'The BAM file to validate a header on.',
        },
    ],
};

sub help_brief {
    'Tool to make sure a TCGA header exists.  Only verifies the existance of @PG record and @RG record';
}

sub help_detail {
    return <<EOS
    Tool to add a read group tag to SAM files.
EOS
}

sub execute {
    my $self = shift;

    my $input_file = $self->input_file;
    my $samtool = $self->samtools_path;
    my $get_header_cmd = "$samtool view -H $input_file";
    my @header_lines = `$get_header_cmd`;
    my @pg_result = grep (/\@PG/,@header_lines);
    my @rg_result = grep (/\@RG/,@header_lines);

    unless ( scalar(@pg_result) > 0 ) {
        $self->debug_message("Input file $input_file does not contain a \@PG record.");
        return;
    }
    unless ( scalar(@rg_result) > 0 ) {
        $self->debug_message("Input file $input_file does not contain a \@RG record.");
        return;
    }

    #check the first read to see if the PG and RG values are present
    my $get_read_cmd = "$samtool view $input_file | head -n1";
    my $get_read_result = `$get_read_cmd`;

    unless ($get_read_result =~ m/RG:Z/) {
        $self->debug_message("Input file $input_file does not contain a RG record.");
        return; 
    }

    unless ($get_read_result =~ m/PG:Z/) {
        $self->debug_message("Input file $input_file does not contain a PG record.");
        return; 
    }

    #$self->debug_message(join("\n",@header_lines));
    return 1;
}


1;
