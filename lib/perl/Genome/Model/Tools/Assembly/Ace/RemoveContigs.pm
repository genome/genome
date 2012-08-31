package Genome::Model::Tools::Assembly::Ace::RemoveContigs;

use strict;
use warnings;

use Genome;
use IO::File;
use Cwd;
use Data::Dumper;

class Genome::Model::Tools::Assembly::Ace::RemoveContigs {
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
        contigs => {
            type => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'Comma separated list of contigs to export',
        },
	contigs_list => {
	    type => 'Text',
	    doc => 'file of list of contig names to export',
	},
	merge => {
	    type => 'Boolean',
	    is_optional => 1,
	    doc => 'merge exported contigs into a single ace file',
	},
	directory => {
	    type => 'Text',
	    is_optional => 1,
	    doc => 'directory where ace files are located',
	},
    ],
};

sub help_brief {
    'Tool to export contig(s) from ace file(s)'
}

sub help_detail {
    return <<"EOS"
gmt assembly ace remove-contigs --ace Felis_catus-3.0.pcap.ace --contigs-list contigs.txt
gmt assembly ace remove-contigs --ace-list acefiles.txt --contigs-list contigs.txt --merge
EOS
}

sub execute {
    my $self = shift;

    $self->directory(cwd()) unless $self->directory;

    $self->status_message("Validating ace files(s)");
    my $acefiles; #array ref
    unless (($acefiles) = $self->get_valid_input_acefiles()) {
	$self->error_message("Failed to validate ace input(s)");
	return;
    }
    $self->status_message("Validating contigs list");
    my $contig_names = {};
    unless (($contig_names) = $self->get_valid_contigs_from_list()) {
	$self->error_message("Failed to validate contigs list");
	return;
    }
    $self->status_message("Filtering contigs in ace file(s)");
    my $new_aces; #array ref
    unless (($new_aces) = $self->filter_ace_files($acefiles, $contig_names, 'remove')) {
	$self->error_message("Failed to parse ace files");
	return;
    }
    
    if ($self->merge) {
	$self->status_message("Merging ace files");
	unless ($self->merge_acefiles($new_aces)) {
	    $self->error_message("Failed to merge ace files");
	    return;
	}
    }

    return 1;
}

1;
