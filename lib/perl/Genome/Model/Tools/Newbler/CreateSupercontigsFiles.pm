package Genome::Model::Tools::Newbler::CreateSupercontigsFiles;

use strict;
use warnings;

use Genome;
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Newbler::CreateSupercontigsFiles {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Newbler assembly directory'
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to include in post assembly',
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
    'Tool to make pcap style supercontigs.fasta and supercontigs.agp files';
}

sub help_detail {
    return <<"EOS"
gmt newbler create-supercontigs-files --assembly-directory /gscmnt/111/newbler_ecoli_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    #make assembly dir consed/edit_dir
    unless ( -d $self->consed_edit_dir ) {
        $self->create_consed_dir;
    }

    #filter contigs my min length
    my $scaffolds = $self->get_scaffolding_info;
    # returns .. while filtering out min_contig contigs
    # and re-adjust gap size when contigs are removed
    #'contig00009 => {
    #      'supercontig' => 'scaffold00001',
    #      'contig_length' => '24537',
    #      'pcap_name' => 'Contig0.9',
    #      'contig_name' => 'contig00009',
    #      'gap_length' => '894'
    #    };
    unless( $scaffolds ) {
        $self->error_message( "Failed to parse newbler scaffolds file" );
        return;
    }

    #hash of array of contig names for each scaffold
    my %scaffold_contigs;

    for my $contig ( sort keys %$scaffolds ) {
        my $pcap_contig_name = $scaffolds->{ $contig }->{pcap_name};
        my ($pcap_supercontig_name) = $pcap_contig_name =~ /Contig(\d+)\.\d+/;
        push @{$scaffold_contigs{$pcap_supercontig_name}}, $contig;
    }

    #fasta reader/writer
    my $f_i = Bio::SeqIO->new( -format => 'fasta', -file => $self->all_contigs_fasta_file );
    my $f_o = Bio::SeqIO->new( -format => 'fasta', -file => '>'.$self->supercontigs_fasta_file );

    #supercontiga agp file
    unlink $self->supercontigs_agp_file;
    my $agp_out = Genome::Sys->open_file_for_writing( $self->supercontigs_agp_file );
    for my $scaffold ( sort {$a <=> $b} keys %scaffold_contigs ) {
        my $c = 0;
        my $supercontig_fasta;
        my $agp_position = 0;
        my $contig_start_pos = 0;
        my $contig_end_pos = 0;
        foreach my $contig ( @{$scaffold_contigs{$scaffold}} ) {
            #print "$contig expected\n";
            my $seq = $f_i->next_seq;
            #skip until we get the contig we want
            if ( not $seq->primary_id eq $contig ) {
                until ( $seq->primary_id eq $contig ) {
                    $seq = $f_i->next_seq;
                }
            }
            unless ( $seq ) {
                $self->error_message( "Did not get fasta seq object .. newbler 454AllContigs.fna may be out of sync with 454Scaffolds.txt file or file may have changed after assembling" );
                return;
            }
            #print $seq->primary_id." got\n";
            $c++;
            #append seq fasta to super contig fasta
            $supercontig_fasta .= uc $seq->seq;
            #contig name/seq length for agp contig line
            my $seq_length = length $seq->seq;
            my $pcap_contig_name = $scaffolds->{ $contig }->{pcap_name};
            #agp start/stop positions
            $contig_start_pos = $contig_end_pos + 1;
            $contig_end_pos = $contig_start_pos + $seq_length - 1;
            #print contig line to agp file
            $agp_out->print(
                'Contig'."$scaffold\t$contig_start_pos\t$contig_end_pos\t".++$agp_position."\tW\t$pcap_contig_name\t1\t$seq_length\t+\n"
            );
            unless ( $c == scalar @{$scaffold_contigs{$scaffold}} ) {
                #append Xs of gap lengths to fasta unless it's the last
                #contig in scaffold which is assigned default gap size
                my $gap_length = $scaffolds->{$contig}->{gap_length};

                #append Xs to supercontig fasta for gaps
                $supercontig_fasta .= ('X' x $gap_length );

                #agp start/stop positions for agp fragment line
                $contig_start_pos = $contig_end_pos + 1;
                $contig_end_pos = $contig_start_pos + $gap_length - 1;

                #print fragment line to agp file
                $agp_out->print(
                    'Contig'."$scaffold\t$contig_start_pos\t$contig_end_pos\t".++$agp_position."\tN\t$gap_length\tscaffold\tyes\tpaired-ends\n"
                );
            }
        }
        my $seq_obj = Bio::Seq->new( -seq => $supercontig_fasta, id => 'Contig'.$scaffold );
        $f_o->write_seq( $seq_obj );
    }

    $agp_out->close;

    return 1;
}

1;
