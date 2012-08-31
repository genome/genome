package Genome::Model::Tools::Assembly::Ace::ExportContigs;

use strict;
use warnings;

use Genome;
use IO::File;
use Cwd;
use Data::Dumper;

class Genome::Model::Tools::Assembly::Ace::ExportContigs {
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
            is_optional => 1,
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

sub help_synopsis {
    return <<"EOS"
gmt assembly ace export-contigs --ace-files Felis_catus-3.0.pcap.ace --contigs Contig2,Contig5.8
gmt assembly ace export-contigs --ace-list acefiles.txt --contigs-list contigs.txt --merge
EOS
}

sub help_detail {
    return <<EOS
This tool reads in a text file of contig names and creates
a new ace file from those contigs.  If there are multiple
input acefiles, one ace file will be created from contigs
exported from that ace file.  If a single output acefile
is needed, --merge option will merge all exported contigs
together into a single ace file.
EOS
}

sub execute {
    my $self = shift;

    $self->directory(cwd()) unless $self->directory;

    my $acefiles; #array ref
    $self->status_message("Validating ace files(s)");
    unless (($acefiles) = $self->get_valid_input_acefiles()) {
	$self->error_message("Failed to validate ace input(s)");
	return;
    }

    $self->status_message("Validating contigs list");
    unless ( $self->contigs_list or $self->contigs ) {
        $self->error_message("You must supply a contig name or a list of contigs");
        return;
    }
    my $contig_names = {};
    unless (($contig_names) = $self->get_valid_contigs_from_list()) {
	$self->error_message("Failed to validate contigs list");
	return;
    }

    $self->status_message("Filtering contigs in ace file(s)");
    my $new_aces; #array ref
    unless (($new_aces) = $self->filter_ace_files($acefiles, $contig_names, 'export')) {
	$self->error_message("Failed to parse ace files");
	return;
    }

    if ($self->merge) {
	$self->status_message("Merging ace files");
	unless ($self->merge_acefiles(acefiles => $new_aces)) {
	    $self->error_message("Failed to merge ace files");
	    return;
	}
    }

    return 1;
}

1;
