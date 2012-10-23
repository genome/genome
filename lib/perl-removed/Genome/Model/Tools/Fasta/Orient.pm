package Genome::Model::Tools::Fasta::Orient;

use strict;
use warnings;

use Genome;

require Bio::SearchIO;
require Bio::Seq;
require Bio::SeqIO;
require Bio::Seq::PrimaryQual;
use Data::Dumper;
require File::Temp;
require Genome::Model::Tools::WuBlast::Blastn;
require Genome::Model::Tools::WuBlast::Xdformat::Create;
require IO::File;

my @SENSE_TYPES = (qw/ sense anti_sense /);

class Genome::Model::Tools::Fasta::Orient {
    is  => 'Genome::Model::Tools::Fasta',
    has_optional => [	 
    map(
        {
            $_.'_fasta_file' => {
                is => 'String',
                is_input => 1,
                doc => ucfirst( join('-', split(/_/, $_)) ).' FASTA file',
            }
        } @SENSE_TYPES
    ),
    ],
};

sub confirmed_fasta_file {
    return $_[0]->fasta_file_with_new_suffix('confirmed');
}

sub confirmed_qual_file {
    return $_[0]->qual_file_with_new_suffix('confirmed');
}

sub unconfirmed_fasta_file {
    return $_[0]->fasta_file_with_new_suffix('unconfirmed');
}

sub unconfirmed_qual_file {
    return $_[0]->qual_file_with_new_suffix('unconfirmed');
}

sub help_brief {
    return 'Orients FASTA (and Quality) files by blastn given sense and anti-sense sequences';
}

sub help_detail { 
    return <<EOS 
    Orients FASTA (and Quality) files by given sense and anti-sense sequences.  Produces 4 files.  A 'confirmed' set of files (fasta and qual) for the sequences which their orientaion is confirmed, and an 'unconfirmed' set of files (fasta and qual) for which the orientation was not confirmed.'
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->sense_fasta_file or $self->anti_sense_fasta_file ) {
        $self->error_message("Need sense or anti-sense fasta files to query");
        return;
    }

    for my $sense_type ( @SENSE_TYPES ) {
        my $fasta_method = $sense_type.'_fasta_file';
        if ( $self->$fasta_method and !-e $self->$fasta_method ) {
            $self->error_message( sprintf("$sense_type FASTA file (%s) does not exist.", $self->$fasta_method) );
            return;
        }
    }

    return $self;
}

sub execute {
    my $self = shift;

    # Temp dir and files
    my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    my $database = sprintf('%s/blast_db', $tmp_dir);
    my $query_file = sprintf('%s/query_file.fasta', $tmp_dir);

    # Create blast db
    my $xdformat = Genome::Model::Tools::WuBlast::Xdformat::Create->create(
        database => $database,
        overwrite_db => 1,
        fasta_files => [ $self->fasta_file ],
    )
        or return;
    $xdformat->execute
        or return;

    # Blast & parse
    my %orientation_confirmation;
    for my $sense_type ( @SENSE_TYPES ) {
        my $fasta_method = $sense_type.'_fasta_file';

        next unless defined $self->$fasta_method;
        
        my $blastn = Genome::Model::Tools::WuBlast::Blastn->create(
            database => $database,
            query_file => $self->$fasta_method, 
            M => 1, # these params optimized for sort sequences.
            N => -3,
            Q => 3,
            R => 1,
            V => 100000,
            B => 100000,
        )
            or return;
        $blastn->execute
            or return;

        my $search_io = Bio::SearchIO->new(
            '-file' => $blastn->output_file,
            '-format' => 'blast',
        );

        while ( my $result = $search_io->next_result ) {
            while( my $hit = $result->next_hit ) {
                while ( my $hsp = $hit->next_hsp ){
                    my $query = $hsp->query;
                    my $orientation_confirmation = 0;
                    if ( $sense_type eq 'sense' ) {
                        # If this hit is on the - strand, it needs to complemeted
                        $orientation_confirmation = 1 if $query->strand == -1;
                    }
                    else { #anti sense
                        # If this hit is on the + strand, it needs to complemeted
                        $orientation_confirmation = 1 if $query->strand == 1;
                    }
                    #my $subject_id = $hsp->subject->seq_id;
                    # TODO Verify?
                    # if ( exists $orientation_confirmation{$subject_id}
                    #        and $orientation_confirmation{$subject_id} != $orientation_confirmation ) {
                    #    $self->error_message();
                    #}
                    $orientation_confirmation{ $hsp->subject->seq_id } = $orientation_confirmation;
                    #last; # TODO which one to last?
                }
            }
        }

        unlink $blastn->output_file if -e $blastn->output_file;
    }

    # Write FASTA and Qual
    $self->_write_fasta_file(\%orientation_confirmation);
    $self->_write_qual_file(\%orientation_confirmation) if $self->have_qual_file;

    $self->status_message( 
        sprintf(
            "<== Confirmed Files ==>\nFasta %s\nQual: %s\n<== Uncomfirmed Files ==>\nFasta%s\nQual:%s\n", 
            $self->confirmed_fasta_file,
            $self->confirmed_qual_file,
            $self->unconfirmed_fasta_file,
            $self->unconfirmed_qual_file,
        )
    );

    return 1;
}

sub _write_fasta_file {
    my ($self, $orientation_confirmation) = @_;

    # Open fasta reading 
    my $bioseq_in = $self->get_fasta_reader( $self->fasta_file )
        or return;
    # Open confirmed fasta for writing 
    my $confirmed_fasta = $self->confirmed_fasta_file;
    unlink $confirmed_fasta if -e $confirmed_fasta;
    my $confirmed_bioseq_out = $self->get_fasta_writer($confirmed_fasta)
        or return;
    # Open unconfirmed fasta for writing 
    my $unconfirmed_fasta = $self->unconfirmed_fasta_file;
    unlink $unconfirmed_fasta if -e $unconfirmed_fasta;
    my $unconfirmed_bioseq_out = $self->get_fasta_writer($unconfirmed_fasta)
        or return;

    while ( my $bioseq = $bioseq_in->next_seq ) { 
        if ( exists $orientation_confirmation->{ $bioseq->id } ) {
            if ( $orientation_confirmation->{ $bioseq->id } ) {
                $bioseq->alphabet('dna') if $bioseq->alphabet eq 'protein';
                $bioseq = $bioseq->revcom;
            }
            $confirmed_bioseq_out->write_seq($bioseq);
        }
        else {
            $unconfirmed_bioseq_out->write_seq($bioseq);
        }
    }

    return 1;
}

sub _write_qual_file {
    my ($self, $orientation_confirmation) = @_;

    # Open qual reading 
    my $bioseq_in = $self->get_qual_reader( $self->qual_file )
        or return;
    # Open confirmed qual for writing 
    my $confirmed_qual = $self->confirmed_qual_file;
    unlink $confirmed_qual if -e $confirmed_qual;
    my $confirmed_bioseq_out = $self->get_qual_writer($confirmed_qual)
        or return;
    # Open unconfirmed qual for writing 
    my $unconfirmed_qual = $self->unconfirmed_qual_file;
    unlink $unconfirmed_qual if -e $unconfirmed_qual;
    my $unconfirmed_bioseq_out = $self->get_qual_writer($unconfirmed_qual)
        or return;

    while ( my $bioseq = $bioseq_in->next_seq ) { 
        if ( exists $orientation_confirmation->{ $bioseq->id } ) {
            if ( $orientation_confirmation->{ $bioseq->id } ) {
                $bioseq = $bioseq->revcom;
            }
            $confirmed_bioseq_out->write_seq($bioseq);
        }
        else {
            $unconfirmed_bioseq_out->write_seq($bioseq);
        }
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
