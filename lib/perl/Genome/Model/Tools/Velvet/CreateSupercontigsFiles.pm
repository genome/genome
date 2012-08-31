package Genome::Model::Tools::Velvet::CreateSupercontigsFiles;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';
use Bio::Seq;
use Bio::SeqIO;
use AMOS::AmosLib;

class Genome::Model::Tools::Velvet::CreateSupercontigsFiles {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Assembly directory',
        },
        min_contig_length => {
            is => 'Text',
            doc => 'Min contig length to process',
        },
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign',
            is_optional => 1,
        }
    ],
};

sub help_brief {
    'Tool to create pcap style supercontigs.agp and supercontigs.fasta files from velvet contigs.fa file'
}

sub help_detail {
    return <<EOS
gmt velvet create-supercontigs-files --assembly-directory /gscmnt/111/velvet_asm --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    #create edit_dir
    unless ( $self->create_edit_dir ) {
	$self->error_message("Assembly edit_dir does not exist and could not create one");
	return;
    }

    #filehandle to print supercontigs.agp file
    unlink $self->supercontigs_agp_file;
    my $agp_fh = Genome::Sys->open_file_for_writing($self->supercontigs_agp_file);

    #IO to output supercontigs.fasta
    my $fa_out = Bio::SeqIO->new(-format => 'fasta', -file => ">".$self->supercontigs_fasta_file);

    #contig/supercontig lengths and gap sizes
    my $scaf_info;
    unless( $scaf_info = $self->get_scaffold_info_from_afg_file ) {
        $self->error_message( "Failed to get scaffolding info from afg file" );
        return;
    }

    my $supercontig_fasta;
    my @ctg_names; 
    #read in contigs from afg file
    my $afg_fh = Genome::Sys->open_file_for_reading($self->velvet_afg_file);
    while (my $record = getRecord($afg_fh)) {
	my ($rec, $fields, $recs) = parseRecord($record);
	if ($rec eq 'CTG') {

            my $contig_seq = $fields->{seq};
            $contig_seq =~ s/\n//g;

            my $contig_name = $fields->{eid};
            $contig_name =~ s/\-/\./;

            # skip scaffolds less than min length
            next if $scaf_info->{$contig_name}->{is_removed};
            my ( $sctg_num, $ctg_num ) = split ('-', $fields->{eid});

            # for contigs >= min length
            if ( length $contig_seq >= $self->min_contig_length ) {
                # append base string to scaffold
                $supercontig_fasta .= $contig_seq;
                my $gap_size = ( exists $scaf_info->{$contig_name}->{processed_gap_size} ) ?
                    $scaf_info->{$contig_name}->{processed_gap_size} :
                    0;
                $supercontig_fasta .= ( 'N' x $gap_size ) if $gap_size;
                # track original contig numbers used
                push @ctg_names, $ctg_num;
            }
            
            #append gap/seq and go to next contig if exists
            #next if exists $scaf_info->{$contig_name}->{next_contig_name};
            next if not $scaf_info->{$contig_name}->{is_last_scaffold_contig};

            #write seq
            my $pcap_name = 'Contig'. ( $sctg_num - 1 );
            my $seq = Bio::Seq->new( -seq => $supercontig_fasta, -id => $pcap_name );
            $fa_out->write_seq( $seq );

            #write agp
            unless( $self->_write_agp($seq, $agp_fh, \@ctg_names) ){
                $self->error_message( "Failed to write agp for $pcap_name" );
                return;
            }
            $supercontig_fasta = ''; #reset for next scaffold
        }
    }
    $agp_fh->close;

    return 1;
}

sub _write_agp {
    my ($self, $seq, $fh, $ctg_names) = @_;

    #writing fasta in prev step already changed seq->id to pcap name
    my $pcap_sctg = $seq->id;

    #put a space between bases and gaps for splitting
    my $string = $self->_separate_bases_and_gaps($seq->seq);
    unless ($string) {
	$self->error_message("Failed to separate bases from gaps for sequence: ".$seq->seq);
	return;
    }

    my $start = 1;
    my $stop = 0;

    my $fragment_order = 0;
    foreach (split(/\s+/, $string)) {
	$fragment_order++;
	if ($_ =~ /^[ACTG]/) { #bases
	    my $length = length $_;
            my $contig_number = shift @$ctg_names;
	    my $contig_name = $pcap_sctg.'.'.++$contig_number;
	    $stop = $start - 1 + $length;
	    #printing: Contig1 1       380     1       W       Contig1.1       1       380     +
	    $fh->print( "$pcap_sctg\t$start\t$stop\t$fragment_order\tW\t$contig_name\t1\t$length\t+\n");
	    $start += $length;
	}
	elsif ($_ =~ /^N/) { #gaps
	    my $length = length $_;
	    $stop = $start - 1 + $length;
	    #printing: Contig1 381     453     2       N       73      scaffold        yes     paried-ends 
	    $fh->print("$pcap_sctg\t$start\t$stop\t$fragment_order\tU\t$length\tscaffold\tyes\tpaired-ends\n");
	    $start += $length;
	}
	else { #probably not necessary
	    $self->error_message("Found none base or gap string in sequence: ".$_);
	    return;
	}
    }

    return 1;
}

sub _separate_bases_and_gaps {
    my ($self, $seq) = @_;
    #TODO - better way?
    $seq =~ s/AN/A N/g;
    $seq =~ s/CN/C N/g;
    $seq =~ s/GN/G N/g;
    $seq =~ s/TN/T N/g;
    $seq =~ s/NA/N A/g;
    $seq =~ s/NC/N C/g;
    $seq =~ s/NG/N G/g;
    $seq =~ s/NT/N T/g;
    return $seq;
}

1;
