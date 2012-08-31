package Genome::Model::Tools::Assembly::Ace::Merge;

use strict;
use warnings;

use Genome;
use IO::File;
use Cwd;
use Data::Dumper;

class Genome::Model::Tools::Assembly::Ace::Merge {
    is => 'Genome::Model::Tools::Assembly::Ace',
    has => [
        ace_files => {
            type => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'comma separated string of ace file names',
        },
        ace_list => {
            type => 'Text',
            is_optional => 1,
            doc => 'file of list of ace files to export contigs from',
        },
        directory => {
            type => 'Text',
            is_optional => 1,
            doc => 'directory where ace files are located',
        },
        increment_contig_numbers => {
            type => 'Integer',
            is_optional => 1,
            default_value => 1000000,
            doc => 'Increment merged contig numbers by this number to avoid same numbered contigs in merged acefile. Use zero (0) to not increment. Using this option may result in corrupt acefiles and/or lost contigs. Use at your own risk.',
        },
    ],
};

sub help_brief {
    'Tool to export contig(s) from ace file(s)'
}

sub help_detail {
    return <<"EOS"
gmt assembly ace merge --ace Felis_catus-3.0.pcap.ace
gmt assembly ace merge --ace-list acefiles.txt  --directory /gscmnt/999/assembly/my_assembly
gmt assembly ace merge --acefile-names file.ace.0,file.ace.2,file.ace.3
EOS
}

sub execute {
    my $self = shift;

    $self->directory(cwd()) unless $self->directory;

    my $acefiles; #array ref
    unless (($acefiles) = $self->get_valid_input_acefiles()) {
        $self->error_message("Failed to validate ace input(s)");
        return;
    }

    my $rv = $self->merge_acefiles(
        ace_files => $acefiles,
        increment => $self->increment_contig_numbers,
    );
    unless ( $rv ) { 
        $self->error_message("Failed to merge ace files");
        return;
    }
    return 1;
}

1;
