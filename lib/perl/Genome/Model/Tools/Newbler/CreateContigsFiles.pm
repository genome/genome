package Genome::Model::Tools::Newbler::CreateContigsFiles;

use strict;
use warnings;

use Genome;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Seq::Quality;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Newbler::CreateContigsFiles {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Newbler assembly directory',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to export',
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
    'Tool to create pcap style contigs.bases and contigs.quals file for newbler assemblies';
}

sub help_detail {
    return <<"EOS"
gmt newbler create-contigs-files --assembly-directory /gscmnt/111/newbler_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    #make edit_dir in assembly dir
    unless ( -d $self->consed_edit_dir ) {
        $self->create_consed_dir;
    }

    #filter contigs by min length
    my $scaffolds = $self->get_scaffolding_info;
    if ( not $scaffolds ) {
        $self->error_message( "Failed to get scaffolding info" );
        return;
    }

    #read in
    my $f_i = Bio::SeqIO->new( -format => 'fasta', -file => $self->all_contigs_fasta_file );
    my $q_i = Bio::SeqIO->new( -format => 'qual', -file => $self->all_contigs_qual_file );
    #print out
    my $f_o = Bio::SeqIO->new( -format => 'fasta', -file => '>'.$self->contigs_bases_file );
    my $q_o = Bio::SeqIO->new( -format => 'qual', -file => '>'.$self->contigs_quals_file );

    my $supercontig_number = 0;
    SEQ: while ( my $seq = $f_i->next_seq ) { #fasta
        while ( my $qual = $q_i->next_seq ) { #qual
            #make sure fasta and qual are in same order
            if ( not $seq->primary_id eq $qual->primary_id ) {
                $self->error_message( "Fasta and qual files are out of order: got from fasta, ".$seq->primary_id.", from quality got, ".$qual->primary_id );
                return;
            }

            next SEQ unless my $new_name = $scaffolds->{ $seq->primary_id }->{pcap_name};

            #write fasta
            my %fasta_params = (
                -seq => $seq->seq,
                -id  => $new_name,
            );
            $fasta_params{-desc} = $seq->desc if $seq->desc;
            my $new_seq = Bio::Seq->new( %fasta_params );
            $f_o->write_seq( $new_seq );

            #write qual
            my %qual_params = (
                -seq  => $seq->seq,
                -qual => $qual->qual,
                -id   => $new_name,
            );
            $qual_params{-desc} = $qual->desc if $qual->desc;
            my $new_qual = Bio::Seq::Quality->new( %qual_params );
            $q_o->write_seq( $new_qual );

            next SEQ;
        }
    }

    return 1;
}

1;
