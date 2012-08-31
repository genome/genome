package Genome::Model::Tools::Newbler::CreateUnplacedReadsFiles;

use strict;
use warnings;

use Genome;
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Newbler::CreateUnplacedReadsFiles {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Newbler assembly directory',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to consider',
        },
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign when newbler does not assign one',
            is_optional => 1,
            default_value => 10,
        }
    ],
};

sub help_brief {
    'Tool to create pcap style reads.unplaced and reads.unplaced.fasta file for newbler assemblies';
}

sub help_detail {
    return <<"EOS"
gmt newbler create-unplaced-reads-files --assembly-directory /gscmnt/111/newbler_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    my $unplaced_reads;
    unless( $unplaced_reads = $self->_unplaced_reads ) {
        $self->error_message( "Failed to create reads.unplaced file" );
        return;
    }

    unless( $self->_unplaced_reads_fasta( $unplaced_reads ) ) {
        $self->error_message( "Failed to create reads.unplaced.fasta file" );
        return;
    }

    return 1;
}

sub _unplaced_reads_fasta {
    my ( $self, $unplaced_reads ) = @_;

    my @input_fastqs = $self->input_fastq_files;
    my $fasta_writer = Bio::SeqIO->new( -format => 'fasta', -file => '>'.$self->reads_unplaced_fasta_file );
    for my $file ( @input_fastqs ) {
        my $reader = Genome::Model::Tools::Sx::FastqReader->create( file => $file );
        while ( my $seq = $reader->read ) {
            my $read_name = $seq->{id};
            if( exists $unplaced_reads->{$read_name} ) {
                my $fasta = Bio::Seq->new( -seq => $seq->{seq}, -id => $seq->{id} );
                $fasta_writer->write_seq( $fasta );
            }
        }
    }

    return 1;
}

sub _unplaced_reads {
    my $self = shift;

    #filter contigs by min length
    my $valid_scaffolds;
    unless( $valid_scaffolds = $self->get_scaffolding_info ) {
        $self->error_message("Failed to get valid scaffolds");
        return;
    }

    #store unplaced reads for later look up .. could be memory hog if there
    #are lots of reads but using array takes really long time
    my %unplaced_reads; 

    my $fh = Genome::Sys->open_file_for_reading( $self->newb_read_status_file );
    unlink $self->reads_unplaced_file;
    my $ru_fh = Genome::Sys->open_file_for_writing( $self->reads_unplaced_file );
    $fh->getline; #skip header
    SEQ: while ( my $line = $fh->getline ) {
        my @tmp = split( /\s+/, $line );
        #$tmp[0] = read name
        #$tmp[1] = reads status Assembled,Singleton,Repeat ..etc
        #$tmp[2] = contig name
        if ( $tmp[1] eq 'Assembled' or $tmp[1] eq 'PartiallyAssembled' ) {
            if ( not exists $valid_scaffolds->{$tmp[2]} ) {
                $ru_fh->print( "* $tmp[0] Filtered\n" );
                $unplaced_reads{$tmp[0]} = 1;
           }
           next SEQ;
        }
        $unplaced_reads{$tmp[0]} = 1;
        $ru_fh->print( " * $tmp[0] $tmp[1]\n" );
    }
    $ru_fh->close;
    $fh->close;

    return \%unplaced_reads;
}

1;
