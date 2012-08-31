package Genome::Model::Tools::Velvet::CreateContigsFiles;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';
use AMOS::AmosLib;

use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;

class Genome::Model::Tools::Velvet::CreateContigsFiles {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Main assembly directory .. above edit_dir',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to process',
        }
    ],
};

sub help_brief {
    'Tool to create contigs.bases and contigs.qual files from velvet afg file';
}

sub help_detail {
    return <<EOS
gmt velvet create-contigs-files --assembly_directory /gscmnt/111/velvet_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    #check to make sure velvet afg file is there
    unless (-s $self->velvet_afg_file ) {
	$self->error_message("Can't find velvet afg file: ".$self->velvet_afg_file );
	return;
    }

    #make edit_dir
    unless ( $self->create_edit_dir ) {
        $self->error_message("assembly edit dir does not exist and could not create one");
        return;
    }
    #get contig/supercontig lengths for contigs
    my $scaf_info;
    unless( $scaf_info = $self->get_scaffold_info_from_afg_file ) {
        $self->error_message( "Failed to get scaffold info from velvet afg file" );
        return;
    }

    #write out contigs.bases file
    my $f_io = Bio::SeqIO->new(-format => 'fasta', -file => ">".$self->contigs_bases_file);
    #write out contigs.quals file
    my $q_io = Bio::SeqIO->new(-format => 'qual', -file => ">".$self->contigs_quals_file);
    #fh to read in afg file
    my $afg_fh = Genome::Sys->open_file_for_reading( $self->velvet_afg_file );

    while (my $record = getRecord($afg_fh)) {
	my ($rec, $fields, $recs) = parseRecord($record);
	if ($rec eq 'CTG') { #contigs
            #get fasta
            my $seq = $fields->{seq};
            $seq =~ s/\n//g; #contig seq is written in multiple lines .. remove end of line
            #skip contig less than min length

            my $contig_name = $fields->{eid};
            $contig_name =~ s/\-/\./;

            # contig removed by min_length filter
            next if exists $scaf_info->{$contig_name}->{is_removed};

            # convert to pcap name
            my $contig_id = $scaf_info->{$contig_name}->{pcap_name};
            #write fasta
	    my $seq_obj = Bio::Seq->new(-display_id => $contig_id, -seq => $seq);
	    $f_io->write_seq($seq_obj);
	    #write qual
	    my $qual = $fields->{qlt};
	    $qual =~  s/\n//g;
	    my @quals;
	    for my $i (0..length($qual)-1) {
                unless (substr($seq, $i, 1) eq '*') {
                    # use default qual of 20
                    push @quals, 20;#ord(substr($qual, $i, 1)) - ord('0');
                }
            }
	    my $qual_obj = Bio::Seq::Quality->new(-display_id => $contig_id, -seq => $seq, -qual => \@quals);
	    $q_io->write_seq($qual_obj);
	}
    }

    $afg_fh->close;
    return 1;
}

1;
