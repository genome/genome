package Genome::Model::Tools::Velvet::CreateGapFile;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use Sort::Naturally;
use Data::Dumper 'Dumper';


class Genome::Model::Tools::Velvet::CreateGapFile {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Assembly build directory',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to process',
        },
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to create gap.txt file from velvet created contigs.fa file';
}

sub help_detail {
    return <<EOS
gmt velvet create-gap-file --assembly-directory /gscmnt/111/assembly/e_coli_velvet_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    unless ( $self->create_edit_dir ) { # refactor
	$self->error_message("assembly edit_dir does not exist and could not create one");
	return;
    }
    
    unlink $self->gap_sizes_file if -s $self->gap_sizes_file;
    my $fh = Genome::Sys->open_file_for_writing( $self->gap_sizes_file );

    for my $contig ( @{$self->contig_infos} ) {
        # gets contig info in assembled order
        if ( exists $contig->{processed_gap_size} ) {
            $fh->print ($contig->{pcap_name}.' '.$contig->{processed_gap_size}."\n");
        }
    }

    $fh->close;

    return 1;
}

1;
