package Finishing::Assembly::Sequence;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my %bases :name(base_string:o) :isa(string);
my %quals :name(qualities:o) :ds(aryref) :empty_ok(1) :default([]);

my %length :name(_length:p) :isa('int');
my %unpad_bases :name(_unpadded_base_string:p) :isa(string);
my %unpad_length :name(_unpadded_length:p) :isa('int');
my %pad_quals :name(_padded_qualities:p) :ds(aryref) :empty_ok(1);
my %pad_pos :name(_pad_pos:p) :ds(aryref) :empty_ok(1);
my %unpad_pos :name(_unpad_pos:p) :ds(aryref) :empty_ok(1);

sub complement
{
    my $self = shift;
    
	my $new_bases = reverse $self->base_string;
	$new_bases =~ tr/actgACTG/tgacTGAC/;
	$self->base_string($new_bases);
    $self->qualities([ reverse( @{ $self->qualities } ) ]);

    $self->_reset_derived_attributes;# affected attributes
    
    return 1;
}

sub _reset_derived_attributes #: PRIVATE
{
    my $self = shift;

    $self->undef_attribute('_unpadded_base_string');
    $self->undef_attribute('_padded_qualities');
    $self->undef_attribute('_pad_pos');
    $self->undef_attribute('_unpad_pos');

    return 1;
}

#- BASES -#
sub padded_base_string
{
    my $self = shift;

    return $self->base_string(@_);
}

sub padded_base_string_no_xs
{
    my $self = shift;

    my $bases = $self->base_string;
    $bases =~ s/x/n/ig;

    return $bases;
}

sub bases_at_position
{
    my ($self, $position, $length) = @_;
    
    $position--;
    $length = 1 unless $length;

    return substr($self->base_string, $position, $length);
}

sub unpadded_base_string
{
    my ($self, $new_string) = @_;

    if ( $new_string )
    {
        return $self->_sync_base_string($new_string);
    }

    return $self->_unpadded_base_string if $self->_unpadded_base_string;

    my $bases = $self->base_string;
    $bases =~ s/\*//g;

    return $self->_unpadded_base_string($bases);
}

sub unpadded_base_string_no_xs
{
    my $self = shift;

    my $bases = $self->unpadded_base_string;
    $bases =~ s/x/n/ig;

    return $bases;
}

sub bases_at_unpadded_position
{
    my ($self, $position, $length) = @_;
    
    $position--;
    $length = 1 unless $length;

    return substr($self->unpadded_base_string, $position, $length);
}

sub padded_base_string_xn_positions
{
    my $self = shift;

    my @tmp = split ('', $self->padded_base_string);
    my @xn_pos;
    my $count = 0;

    foreach my $base (@tmp)
    {
	$count++;
	if ($base =~ /^[xn]$/i)
	{
	    push @xn_pos, $count;
	}
    }
    return \@xn_pos;
}

#- LENGTH -#
sub length
{
    my $self = shift;

    #my $length = $self->_length;  In some cases, this is not synced with base_string
    #return $length if defined $length;

    return $self->_length( CORE::length($self->base_string) );
}

sub unpadded_length
{
    my $self = shift;

    my $unpad_length = $self->_unpadded_length; 
    
    return $unpad_length if $unpad_length;
    
    return $self->_unpadded_length( CORE::length($self->unpadded_base_string) );
}

#- QUALITY -#
sub unpadded_qualities
{
    my $self = shift;

    return $self->qualities(@_);
}

sub qualities_at_unpadded_position
{
	my ($self, $start, $length) = @_;

    $start--;
    my $stop = $start + ( $length || 1 ) - 1;
    my $quals = $self->qualities;
    my @quals;
    push @quals, $quals->[$_] for $start..$stop; 

    return \@quals;
}

sub padded_qualities
{
    my ($self, $new_quals) = @_;

    if ( $new_quals )
    {
        return $self->_sync_qualities($new_quals);
    }
    
    return [] unless @{ $self->qualities };
    
    $self->_set_conversion_arrays;

    return $self->_padded_qualities;
}

sub qualities_at_padded_position
{
	my ($self, $start, $length) = @_;
    
    $start--;
    my $stop = $start + ( $length || 1 ) - 1;
    my $padded_quals = $self->padded_qualities;
    my @quals;
    push @quals, $padded_quals->[$_] for $start..$stop; 

    return \@quals;
}

#- POSITIONS -#
sub unpad_position_to_pad_position
{
    my ($self, $pos) = @_;

    $self->_set_conversion_arrays;

    return $self->_unpad_pos->[$pos];
}

sub pad_position_to_unpad_position
{
	my ($self, $pos) = @_;

    $self->_set_conversion_arrays;

    return $self->_pad_pos->[$pos];
}

#- CONVERSION ARRAYS -#
sub _set_conversion_arrays
{
    my $self = shift;

    return 1 if $self->_pad_pos 
        and $self->_unpad_pos
        and $self->_padded_qualities;

    my $quals = $self->qualities;
    my (@pad_pos, @unpad_pos, @pad_quals);
    my $pad = 0;
    my $unpad = 0;
    foreach my $base ( split(//, $self->padded_base_string) )
    {
        if ( $base eq '*' )
        {
            $pad_quals[$pad] = int( (($quals->[$unpad - 1] + $quals->[$unpad]) / 2 ) + 0.5 );
            # prev_qual + next_qual / 2 
            # add .5 and use int to round
        }
        else
        {
            $pad_quals[$pad] = $quals->[ $unpad++ ];
        }
        $pad++;
        $pad_pos[$pad] = $unpad;
        $unpad_pos[$unpad] = $pad;
    }

    $self->_padded_qualities(\@pad_quals);
    $self->_pad_pos(\@pad_pos);
    $self->_unpad_pos(\@unpad_pos);

    return 1;
}

sub _sync_base_string
{
    my ($self, $unpadded_string) = @_;

    if ( $self->base_string )
    {
        $self->fatal_msg("New unpadded base string is a different length than the current stored unpadded base string")
        unless CORE::length($unpadded_string) eq $self->unpadded_length;
    }
    else
    {
        $self->base_string($unpadded_string);
        return $self->_unpadded_base_string($unpadded_string);
    }

    my @unpadded_bases = split(//, $unpadded_string);

    my $new_padded_bases;
    my $i = 0;
    foreach my $old_base ( split(//, $self->base_string) )
    {
        if ( $old_base eq '*' )
        {
            $new_padded_bases .= '*';
        }
        else
        {
            $new_padded_bases .= $unpadded_bases[$i];
            $i++;
        }
    }

    $self->base_string($new_padded_bases);
    
    return $self->_unpadded_base_string($unpadded_string);
}

sub _sync_qualities
{
    my ($self, $padded_qualities) = @_;

    $self->fatal_msg("No base string to use as de-padding reference when setting padded qualities") unless $self->base_string;

    if (@$padded_qualities) {
        $self->fatal_msg("New padded qualities is a different length than the current padded base string") 
            unless @$padded_qualities eq $self->length;
    }

    my @qualities;
    my $i = 0;
    foreach my $base ( split(//, $self->base_string) )
    {
        push @qualities, $padded_qualities->[$i] unless $base eq '*';
        $i++;
    }

    $self->qualities(\@qualities);
    
    return $self->_padded_qualities($padded_qualities);
}

1;

#$HeadURL$
#$Id$
