package VarScan::FisherTest;

use warnings;
use strict;

=head1 NAME

VarScan::FisherTest

=head1 VERSION

    Version 1.03

=head1 SYNOPSIS

    This module

=cut

our $VERSION = '1.03';
use Math::Pari ();

my $Tolerance = 1;
$Tolerance /= 2 while 1 + $Tolerance/2 > 1;


=head1 FUNCTIONS

=cut

################################################################################################
=head2 calculate_p_value - Given two sets of read counts, determine a p-value
#
################################################################################################
=cut

sub calculate_p_value {
    my ( $a, $b, $c, $d, $ts ) = @_;
    my $test = $a*( $a + $b + $c + $d ) - ( $a + $b )*( $a + $c );
    
    return 1 if $test < 0 and $ts;
    # below here, $test < 0 implies !$ts;
    
    my $p_val;

    if ( $test < 0 )
    {
        if ( $d < $a )
        {
            $p_val = _fishers_exact( $d, $c, $b, $a, $ts, 1 );
        }
        else
        {
            $p_val = _fishers_exact( $a, $b, $c, $d, $ts, 1 );
        }
    }
    else
    {
        if ( $b < $c )
        {
            $p_val = _fishers_exact( $b, $a, $d, $c, $ts, 0 );
        }
        else
        {
            $p_val = _fishers_exact( $c, $d, $a, $b, $ts, 0 );
        }
    }

    return $p_val;
}

sub _fishers_exact {
    my ( $a, $b, $c, $d, $ts, $complement ) = @_;
    die "Bad args\n" if $ts && $complement;
    
    my ( $aa, $bb, $cc, $dd ) = ( $a, $b, $c, $d );
    my $first = my $delta = _single_term( $aa, $bb, $cc, $dd );
    my $p_val = 0;

    {
        $p_val += $delta;
        last if $aa < 1;
        $delta *= ( ( $aa-- * $dd-- )/( ++$bb * ++$cc ) );
        redo;
    }

    if ( $ts )
    {
        my $m = $b < $c ? $b : $c;
        ($aa, $bb, $cc, $dd ) = ( $a + $m, $b - $m, $c - $m, $d + $m );
        $delta = _single_term( $aa, $bb, $cc, $dd );
        my $bound = -$Tolerance;
        while ( $bound <= ( $first - $delta )/$first && $aa > $a )
        {        
            $p_val += $delta;
            $delta *= ( ( $aa-- * $dd-- )/( ++$bb * ++$cc ) );
        }
    }
    elsif ( $complement )
    {
        $p_val = 1 - $p_val + $first;
    }

    return $p_val;
}

sub _single_term
{
    my ( $a, $b, $c, $d ) = @_;
    my ( $r1, $r2 ) = ($a + $b, $c + $d);
    my ( $c1, $c2 ) = ($a + $c, $b + $d);
    my $N = $r1 + $r2;
    
    return  exp( _ln_fact( $r1 ) + _ln_fact( $r2 ) +
                 _ln_fact( $c1 ) + _ln_fact( $c2 ) -
                 _ln_fact( $N ) -
                 ( _ln_fact( $a ) + _ln_fact( $b ) +
                   _ln_fact( $c ) + _ln_fact( $d ) ) );
}


{
    my $two_pi;
    my $pi_over_3;
    my $half;
    
    BEGIN {
        $two_pi    = Math::Pari::PARI( 2 * atan2 0, -1 );
        $pi_over_3 = Math::Pari::PARI( atan2( 0, -1 )/3 );
        $half      = Math::Pari::PARI( 0.5 );
    }

    sub _ln_fact {
        my $n = Math::Pari::PARI( shift() );
        die "Bad args to _ln_fact: $n" if $n < 0;
        my $ln_fact;
        eval {
            $ln_fact = log Math::Pari::factorial( $n );
        };

        if ( $@ ) {
            die $@ unless $@ =~ /\QPARI:   ***   exponent overflow/;
            # Gosper's approximation; from
            # http://mathworld.wolfram.com/StirlingsApproximation.html
            $ln_fact = $half * log( $two_pi*$n + $pi_over_3 )
            +  $n * log( $n )
            -  $n;
        }
        return $ln_fact;
    }
}


################################################################################################
=head2 format_p_value - Format the p-value to a reasonable length
#
################################################################################################
=cut

sub format_p_value
{
    my $p_val = shift(@_);
    
    ## Reformat p-value to consistent scientific notation ##
    
    $p_val = sprintf("%E", $p_val) if($p_val < 0.0001);
    
    if($p_val =~ 'E')
    {
        (my $numeric, my $exponent) = split(/E/, $p_val);
        $numeric = sprintf("%.1f", $numeric);
        $p_val = $numeric . 'E' . $exponent;
    }
    else
    {
            $p_val = sprintf("%.5f", $p_val);
    }
    
    return($p_val);
}


=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Center at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use.

=cut

1; # End of VarScan::FisherTest
