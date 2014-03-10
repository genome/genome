package Genome::Model::Tools::Sam::SortBam;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::SortBam {
    is  => 'Genome::Model::Tools::Sam',
    has => [
	file_name => {
	    is  => 'String',
	    doc => 'Input File',
	},
	name_sort => {
	    is => 'Boolean',
	    is_optional => 1,
	    doc => 'Sort by read names instead of chrom coordinates'
	},
    output_file => {
        is  => 'String',
        doc => 'Output file, defaults to file_name.sorted.bam',
        is_optional => 1,
    },
    ],
};

sub help_brief {
    'Tool to sort BAM files';
}

sub help_detail {
    return <<EOS
    Tool to sort BAM files.
EOS
}



sub execute {
    my $self = shift;

    #samtools sort does dumb things with the file extension.  must snap off the output filename.

    my $real_output_file; 

    if ($self->output_file) {
        $real_output_file = $self->output_file;
        $real_output_file =~ s/\.bam$//;
    } else {
        $real_output_file = $self->file_name;
        $real_output_file =~ s/\.bam$//;
        $real_output_file .= ".sorted";	
    }

    $self->debug_message(sprintf("attempting to sort %s into file with prefix %s", $self->file_name, $real_output_file));

    my $sort_params = ($self->name_sort ? " -n " : "");
    $sort_params .= ' -m '. $self->maximum_memory;

    my $sam_sort_cmd = sprintf("%s sort %s %s %s", $self->samtools_path, $sort_params, $self->file_name,  $real_output_file);
    my $sam_sort_rv = Genome::Sys->shellcmd(cmd=>$sam_sort_cmd, input_files=>[$self->file_name], output_files=>[$real_output_file.".bam"], skip_if_output_is_present=>0);
    if ($sam_sort_rv != 1) {
        $self->error_message("Bam sort error.  Return value $sam_sort_rv");
        return;
    }

    return 1;
}


1;
