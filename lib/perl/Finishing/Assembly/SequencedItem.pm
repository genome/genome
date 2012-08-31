package Finishing::Assembly::SequencedItem;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Item';

use Data::Dumper;
use Finishing::Assembly::Sequence;

my %seq :name(_sequence:p) :isa('object');

#- HACK FOR OLD INTERFACE -#
sub sequence
{
    my $self = shift;
    $self->warn_msg("Accessing the sequence directly is deprecated");
    return $self;
}

#- SEQUENCE -#
sub START
{
    my $self = shift;

    my $base_string = $self->proxy->get_method('base_string')->();
    my $qualities = $self->proxy->get_method('qualities')->();
    return $self->_sequence
    (
        Finishing::Assembly::Sequence->new
        (
            base_string => $base_string || undef,
            qualities => $qualities || undef,
        )
    );

    return 1;
}

sub base_string
{
    my ($self, $new_bases) = @_;

    if ( $new_bases ) # update the data source
    {
        $self->_set_base_string($new_bases);
    }

    return $self->_sequence->base_string($new_bases); # update the sequence
}

sub get_base_at_position
{
    my ($self, $position) = @_;

    $self->fatal_msg("You must supply a base position to get that base") unless
	$position =~ /^\d+$/;

    my @bases = split ('', $self->_sequence->base_string);

    return $bases[$position];
}


sub _set_base_string #: PRIVATE
{
    my ($self, $new_bases) = @_;

    my $method = $self->proxy->get_method('base_string', $new_bases);
    #$self->fatal_msg("Can't get base_string method from proxy") unless $method;

    return $method->();
}

sub qualities
{
    my ($self, $new_quals) = @_;

    if ( $new_quals ) # update the data src
    {
        $self->_set_qualities($new_quals);
    }

    return $self->_sequence->qualities($new_quals) # update the sequence
}

sub _set_qualities #: PRIVATE
{
    my ($self, $new_quals) = @_;

    my $method = $self->proxy->get_method('qualities', $new_quals);
    #$self->fatal_msg("Can't get qualities method from proxy") unless $method;
    
    return $method->();
}

sub length
{
    my ($self, $length) = @_;

    my $method = $self->proxy->get_method('length', $length);

    $self->undef_attribute('_seq') if $length;

    return $method->() if $method;

    return $self->_sequence->length;
}

#- COMPLEMENT -#
sub complement : CUMULATIVE
{
    my $self = shift;

    $self->_sequence->complement;
    $self->base_string( $self->_sequence->base_string );
    $self->qualities( $self->_sequence->qualities );
    $self->complemented( not $self->complemented );

    return 1;
}

#- LENGTH -#
sub base_count
{
    return shift->_sequence->length;
}

sub padded_base_count
{
    return shift->_sequence->length;
}

sub padded_length
{
    return shift->_sequence->length;
}

sub unpadded_length
{
    return shift->_sequence->unpadded_length;
}

sub unpadded_base_count
{
    my $self = shift;

    return $self->_sequence->unpadded_length;
}

#- BIOSEQ -#
sub to_bioseq
{
    my ($self, %p) = @_;
    # p has keys: start, stop, type, no_xs

    my ($bases_method, $quals_method, $length_method) = ( exists $p{type} and $p{type} eq 'padded' )
    ? ('padded_base_string', 'qualities_at_padded_position', 'padded_length')
    : ('unpadded_base_string', 'qualities_at_unpadded_position', 'unpadded_length');

    my $start = $p{start} || 1;
    my $stop = $p{stop} || $self->$length_method;
    
    my $bases = uc substr
    (
        $self->$bases_method, 
        $start - 1,
        $stop - $start + 1
    );
    $bases =~ s/x/N/ig if $p{no_xs};

    return Bio::Seq::Quality->new
    (
        '-id' =>  $self->name,
        '-desc' => "$start to $stop",
        '-seq' => $bases,
        '-qual' => join(" ", @{ $self->$quals_method($start, $stop - $start + 1) }),
        '-alphabet' => 'dna'
    );
}

#- BASES -#
sub padded_base_string
{
    return base_string(@_);
}

sub padded_base_string_no_xs
{
    my $self = shift;

    return $self->_sequence->padded_base_string_no_xs;
}

sub unpadded_base_string
{
    my ($self, $bases) = @_;

    my $unpad_bases = $self->_sequence->unpadded_base_string($bases);

    $self->_set_base_string( $self->_sequence->base_string ) if $bases;

    return $unpad_bases;
}

sub unpadded_base_string_no_xs
{
    my $self = shift;

    return $self->_sequence->unpadded_base_string_no_xs;
}

sub padded_base_string_xn_positions
{
    my $self = shift;
    return $self->_sequence->padded_base_string_xn_positions;
}

#- QUALITY -#
sub unpadded_qualities
{
    return qualities(@_);
}

sub unpadded_base_qualities
{
    return qualities(@_);
}

sub unpadded_base_quality
{
    return qualities(@_);
}

sub qualities_at_unpadded_position
{
    my $self = shift;

    return $self->_sequence->qualities_at_unpadded_position(@_);
}

sub padded_qualities
{
    my ($self, $new_quals) = @_;

    my $quals = $self->_sequence->padded_qualities($new_quals);

    $self->_set_qualities( $self->_sequence->qualities ) if $new_quals;

    return $quals;
}

sub padded_quality
{
    return padded_qualities(@_);
}

sub padded_base_qualities
{
    return padded_qualities(@_);
}

sub padded_base_quality
{
    return padded_qualities(@_);
}

sub qualities_at_padded_position
{
    my $self = shift;

    return $self->_sequence->qualities_at_padded_position(@_);
}

#- POSITIONS -#
sub unpad_position_to_pad_position
{
    my $self = shift;

    return $self->_sequence->unpad_position_to_pad_position(@_);
}

sub pad_position_to_unpad_position
{
    my $self = shift;

    return $self->_sequence->pad_position_to_unpad_position(@_);
}

1;

#$HeadURL$
#$Id$
