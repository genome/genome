package Genome::Model::Tools::Velvet::CreateReadsFiles;

use strict;
use warnings;

use Genome;
use IO::File;
use IO::Seekable;
use Bio::SeqIO;
use AMOS::AmosLib;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Velvet::CreateReadsFiles {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
	assembly_directory => {
	    is => 'Text',
	    doc => 'Assembly directory',
	},
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to export reads for',
        },
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign',
            is_optional => 1,
        }
    ],
};

sub help_brief {
    'Tool to create pcap style readinfo.txt and reads.placed files'
}

sub help_detail {
    return <<EOS
gmt velvet create-reads-files --assembly_directory /gscmnt/111/velvet_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    #make edit_dir
    unless ( $self->create_edit_dir ) {
	$self->error_message("Assembly edit_dir does not exist and could not create one");
	return;
    }

    #readinfo.txt
    unlink $self->read_info_file;
    my $ri_fh = Genome::Sys->open_file_for_writing( $self->read_info_file );

    #reads.placed file
    unlink $self->reads_placed_file;
    my $rp_fh = Genome::Sys->open_file_for_writing( $self->reads_placed_file );

    #velvet output Sequences file
    my $seq_fh = Genome::Sys->open_file_for_reading( $self->velvet_sequences_file );

    #velvet output afg file handle
    my $afg_fh = Genome::Sys->open_file_for_reading( $self->velvet_afg_file );

    #stores read names and seek pos in hash or array indexed by read index
    my $seek_positions = $self->load_sequence_seek_positions( $self->velvet_sequences_file );
    unless ($seek_positions) {
	$self->error_message("Failed to load read names and seek pos from Sequences file");
	return;
    }

    my $scaf_info;
    #gets velvet named contigs names .. w/o min length filtering
    unless( $scaf_info = $self->get_scaffold_info_from_afg_file ) {
        $self->error_message("Failed to get scaffold info from contigs.fa file");
        return;
    }

    my $contig_number = 0;
    while (my $record = getRecord($afg_fh)) {
	my ($rec, $fields, $recs) = parseRecord($record);
	#iterating through contigs
	if ($rec eq 'CTG') {
            my $seq = $fields->{seq};
            $seq =~ s/\n//g;
            my $contig_length = length $seq;

            #filter contigs less than min length
            my $velvet_contig_name = $fields->{eid};
            $velvet_contig_name =~ s/\-/\./;

            # contig filtered out by min_length
            next if $scaf_info->{$velvet_contig_name}->{is_removed};

            my $sctg_start = $scaf_info->{$velvet_contig_name}->{contig_start_position};
            my $pcap_name = $scaf_info->{$velvet_contig_name}->{pcap_name};

	    #iterating through reads
	    for my $r (0 .. $#$recs) {
		my ($srec, $sfields, $srecs) = parseRecord($recs->[$r]);
		if ($srec eq 'TLE') {
		    #sfields:
		    #'src' => '19534',  #read id number
		    #'clr' => '0,90',   #read start, stop 0,90 = uncomp 90,0 = comp
		    #'off' => '75'      #read off set .. contig start position

		    # read start, stop pos and orientaion
		    my ($ctg_start, $ctg_stop, $c_or_u) = $self->_read_start_stop_positions($sfields); 

                    #read seek position in velvet Sequences file
                    my $seek_pos = ${$seek_positions}[$sfields->{src}];
                    unless ( defined $seek_pos ) {
                        $self->error_message("Failed to get read name and/or seek position for read id: ".$sfields->{src});
			return;
		    }

                    # read name and length from Sequences file
		    my ($read_length, $read_name) = $self->_read_length_and_name_from_sequences_file( $seek_pos, $seq_fh );

		    # print to readinfo.txt file
		    $ri_fh->print("$read_name $pcap_name $c_or_u $ctg_start $read_length\n");

		    # convert C U to 1 0 for reads.placed file
		    $c_or_u = ($c_or_u eq 'U') ? 0 : 1;

		    # calculate contig start pos in supercontig
                    my $read_sctg_start = $sctg_start + $ctg_start - 1;
                    # print to reads.placed file
                    my ($pcap_supercontig_number) = $pcap_name =~ /Contig(\d+)\./;
		    $rp_fh->print("* $read_name 1 $read_length $c_or_u $pcap_name Supercontig$pcap_supercontig_number $ctg_start $read_sctg_start\n");
		}
	    }
	}
    }

    $seq_fh->close;
    $afg_fh->close;
    $ri_fh->close;
    $rp_fh->close;

    return 1;
}

sub _read_start_stop_positions {
    my ($self, $fields) = @_;

    my ($start, $stop) = split(',', $fields->{clr});
    unless (defined $start and defined $stop) {
	$self->error_message("Failed to get read start, stop positions from record: ".$fields->{clr});
	return;
    }
    #read complementation
    my $c_or_u = ($start > $stop) ? 'C' : 'U';
    #re-direct start, stop to physical contig positions .. regardless of u or c
    ($start, $stop) = $start < $stop ? ($start, $stop) : ($stop, $start);
    $start += $fields->{off} + 1;
    $stop += $fields->{off} + 1;
    
    return $start, $stop, $c_or_u;
}

sub _read_length_and_name_from_sequences_file {
    my ($self, $seek_pos, $fh ) = @_;

    $fh->seek( $seek_pos, 0 );
    my $io = Bio::SeqIO->new( -fh => $fh, -format => 'fasta', -noclose => 1 );
    my $seq = $io->next_seq;
    unless ( $seq ) {
        $self->error_message("Failed to get seq object at seek position $seek_pos");
        return;
    }

    return length $seq->seq, $seq->primary_id;
}

1;
