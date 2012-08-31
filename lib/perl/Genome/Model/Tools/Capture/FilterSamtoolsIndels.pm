
package Genome::Model::Tools::Capture::FilterSamtoolsIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FilterSamtoolsIndels - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/23/2009 by D.K.
#	MODIFIED:	10/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::FilterSamtoolsIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of variants in indel format", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive filtered indels", is_optional => 0, is_input => 1, is_output => 1 },
		min_coverage     => { is => 'Text', doc => "Minimum coverage to allow indel [4]", is_optional => 1 },
		min_reads2     => { is => 'Text', doc => "Minimum read support to allow indel [2]", is_optional => 1 },
		min_var_freq     => { is => 'Text', doc => "Minimum allele frequency to allow indel [0.10]", is_optional => 1 },
		min_indel_score     => { is => 'Text', doc => "Minimum allele frequency to allow indel [100]", is_optional => 1 },		
		max_var_freq     => { is => 'Text', doc => "Maximum allele frequency to allow indel [1.0]", is_optional => 1 },
		verbose		=> { is => 'Text', doc => "Print filter-passing indels as they come", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filters glfSomatic indel calls"                 
}

sub help_synopsis {
    return <<EOS
This command filters indels from somaticSniper 
EXAMPLE:	gmt analysis somatic-pipeline merge-snvs-with-annotation --variants-file [file] --annotation-file [file] --output-file [file]
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $variants_file = $self->variants_file;
	my $output_file = $self->output_file;
	my $min_coverage = 4;
	$min_coverage = $self->min_coverage if($self->min_coverage);

	my $min_reads2 = 2;
	$min_reads2 = $self->min_reads2 if($self->min_reads2);

	my $min_indel_score = 100;
	$min_indel_score = $self->min_indel_score if($self->min_indel_score);

	my $min_var_freq = 0.10;
	$min_var_freq = $self->min_var_freq if($self->min_var_freq);

	my $max_var_freq = 1.0;
	$max_var_freq = $self->max_var_freq if(defined($self->max_var_freq));

	my %stats = ();
	$stats{'num_indels'} = $stats{'num_pass_filter'} = 0;
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	## Parse the variants file ##
	
	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my @lineContents = split(/\t/, $line);
		my $indel_score = $lineContents[5];
		my $reads1 = $lineContents[11];
		my $reads2 = $lineContents[10];
		(my $allele1, my $allele2) = split(/\//, $lineContents[3]) if($lineContents[3] && $lineContents[3] =~ '/');

		## Adjust for "formatted" indels ##

		if($lineContents[2] =~ /[0-9]/)
		{
			$allele1 = $lineContents[3];
			$allele2 = $lineContents[4];
			$indel_score = $lineContents[6];
			$reads1 = $lineContents[12];
			$reads2 = $lineContents[11];			
		}

		if($reads1 || $reads2)
		{
			my $coverage = $reads1 + $reads2;
			my $var_freq = $reads2 / ($reads1 + $reads2);
			
			$stats{'num_indels'}++;
			
			if(($allele1 && $allele1 =~ 'N') || ($allele2 && $allele2 =~ 'N'))
			{
				## Filter out N-allele indels ##	
			}
			elsif($coverage >= $min_coverage && $var_freq >= $min_var_freq && $reads2 >= $min_reads2 && $indel_score >= $min_indel_score)
			{
				if($self->verbose)
				{
					$var_freq = sprintf("%.2f", $var_freq * 100) . '%';
					print "$line\t$var_freq\n"
		#			print "$chrom\t$chr_start\t$chr_stop\t$indel_type-$var_allele\t$normal_reads1 $normal_reads2\t$tumor_reads1 $tumor_reads2 $tumor_freq\% $p_value $somatic_score\n";
				}
			
				print OUTFILE "$line\t$var_freq\n";
				
				$stats{'num_pass_filter'}++;
				
			}
		}


	}

	close($input);
		
	close(OUTFILE);
	
	print "$stats{'num_indels'} indels\n";
	print "$stats{'num_pass_filter'} passed filter\n";
}



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


1;

