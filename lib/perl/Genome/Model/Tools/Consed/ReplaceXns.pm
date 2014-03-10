package Genome::Model::Tools::Consed::ReplaceXns;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Consed::ReplaceXns {
    is => 'Command',
    has => [
        ace_in => {
            is  => 'Text',
            doc => 'Input ace file',
        },
        ace_out => {
            is => 'Text',
            doc => 'Output ace file name, not not specified is give name <ACE_IN>.xns_replaced',
            is_mutable => 1,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to replace consensus Ns and Xs with a A C G or T base best represented by reads reads at that consensus position'
}

sub help_synopsis {
    return <<EOS
        gmt consed replace-xns --ace-in /gscmnt/100/assembly/edit_dir/454Contigs.ace
        gmt consed replace-xns --ace-in /gscmnt/100/assembly/edit_dir/454Contigs.ace -ace-out /gscmnt/100/assembly/edit_dir/454Contigs.ace.xns_replaced
EOS
}

sub execute {
    my $self = shift;

    #check ace in
    if ( not -s $self->ace_in ) {
        $self->error_message("Failed to find input ace or file is zero size: ".$self->ace_in);
        return;
    }

    #set ace out
    if ( not $self->ace_out ) {
        $self->ace_out( $self->ace_in.'.xns_replaced' );
    }

    #ace reader
    my $reader = Genome::Model::Tools::Consed::AceReader->create(
        file => $self->ace_in,
    );
    if ( not $reader ) {
        $self->error_message("Failed to create ace reader for file: ".$self->ace_in);
    }

    #ace writer
    my $writer = Genome::Model::Tools::Consed::AceWriter->create(
        file => $self->ace_out,
    );
    if ( not $writer ) {
        $self->error_message("Failed to creater ace writer for file: ".$self->ace_out);
        return;
    }
    
    #iterate through contigs
    while ( my $contig = $reader->next_contig ) {
        my @padded_consensus = split( '', $contig->{consensus} );

        #find XNs in consensus
        if ( my @xns_positions = $self->xn_positions_from_consensus(\@padded_consensus) ) {
            $self->debug_message('Processing '.$contig->{name});
            for my $xns_pos_base ( @xns_positions ) {
                my ( $xns_pos, $xns_base ) = split ( ' ', $xns_pos_base );
                $self->debug_message("Found $xns_base at consensus position $xns_pos");

                #pick new consensus from bases called by reads at cons position
                my $new_base = $self->new_base_for_xns($contig, $xns_pos);
                $self->debug_message("Got $new_base to replace $xns_base");
                
                #update padded consensus with new base
                my $new_padded_cons = $self->update_padded_consensus(
                    \@padded_consensus,
                    $xns_pos,
                    $new_base
                );
                $contig->{consensus} = $new_padded_cons;

                #convert padded pos to unpadded pos
                my $unpadded_position = $self->padded_to_unpadded_position(
                    \@padded_consensus,
                    $xns_pos
                );

                #update unpadded consensus
                my $new_unpadded_cons = $self->update_unpadded_consensus(
                    $contig->{unpadded_consensus},
                    $unpadded_position,
                    $new_base
                ); 
                $contig->{unpadded_consensus} = $new_unpadded_cons;

                #update consensus quality to zero
                @{$contig->{base_qualities}}[$unpadded_position - 1] = 0;
            }
        }
        
        #write contig to new ace out
        $writer->add_contig( contig => $contig );
    }

    #transfer tags if any
    if ( $reader->contig_tags ) {
        $self->debug_message("Transfering contig tags");
        $writer->add_contig_tags( $reader->contig_tags );
    }

    $self->debug_message("Writing new ace file .. may take a few minutes");

    $writer->close;

    return 1;
}

sub update_unpadded_consensus {
    my ( $self, $consensus, $pos, $base ) = @_;

    #make sure base at unpadded position is xn
    my $existing_base = substr( $consensus, $pos - 1, 1 );
    if ( not $existing_base =~ /^[xn]$/i ) {
        $self->error_message("Expected X or N at unpadded position, $pos, but got $existing_base");
        return;
    }

    my @unpadded_consensus = split( '', $consensus );
    $unpadded_consensus[$pos - 1] = $base;

    return join ('', map {$_} @unpadded_consensus);
}

sub update_padded_consensus {
    my ( $self, $consensus, $pos, $base ) = @_;

    @$consensus[$pos - 1] = $base;

    return join ('', map {$_} @$consensus );
}

sub new_base_for_xns {
    my ( $self, $contig, $xns_pos ) = @_;

    my $reads = $contig->{reads};
    my %base_counts;
    for my $read ( keys %$reads ) {
        #skip reads not aligned to xns position
        next unless $reads->{$read}->{start} < $xns_pos and $reads->{$read}->{stop} > $xns_pos;
                    
        my $read_base_position = $xns_pos - $reads->{$read}->{start};
        my $read_base = substr( $reads->{$read}->{sequence}, $read_base_position, 1 );

        #track # of times a base is called by a read to pick the most freq represented one
        next if not $read_base =~ /^[acgt]$/i;
        $base_counts{$read_base}++;
    }

    #pick bases most represented
    my @candidate_bases;
    my $prev_count = 0;

    foreach my $base( sort {$base_counts{$b} <=> $base_counts{$a}} keys %base_counts ) {
        next if $base_counts{$base} < $prev_count;
        push @candidate_bases, $base;
        $prev_count = $base_counts{$base};
    }
    
    #randomly select a base from >= 1 possible bases .. selecting from all 4 if no read shows a base
    @candidate_bases = qw/ a c g t / if not @candidate_bases;

    return $candidate_bases[ int rand (scalar @candidate_bases) ];
}

sub xn_positions_from_consensus {
    my ( $self, $padded_consensus ) = @_;

    my $base_position = 0;
    my @xns_positions;

    for my $base ( @$padded_consensus ) {
        $base_position++;
        if ( $base =~ /^[xn]$/i ) {
            push @xns_positions, $base_position.' '.$base;
        }
    }

    return @xns_positions;
}

sub padded_to_unpadded_position {
    my ( $self, $padded_bases, $padded_position ) = @_;

    my $pads_count = 0;

    for ( my $i = 0; $i < $padded_position; $i ++ ) {
        $pads_count++ if @$padded_bases[$i] eq '*';
    }

    return $padded_position - $pads_count;
}

1
