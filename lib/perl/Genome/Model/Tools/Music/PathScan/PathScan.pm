package Genome::Model::Tools::Music::PathScan::PathScan;

# DEBUG
# print "USING LOCAL MCW VERSION OF REGULAR +/-\n";
# DEBUG

#__STANDARD PERL PACKAGES
use strict;
use Carp;

#__CONSTANT OF PI -- NEEDED IN RAMANUJAN APPROX FOR POISSON PROBABILITY MASSES
#  (SEE "PATH-SCAN TEST" NOTES PP 29-31)
use constant PI => 4*atan2 1, 1;
use constant LOG_PI_OVER_2 => log (PI) / 2;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

PathScan - the Path-Scan significance test for mutations in
groups of putative cancer genes

=head1 SYNOPSIS

	use PathScan;

	my $pmobj = PathScan->new ($list_of_gene_lengths);
	my $pval = $pmobj->path_scan ($actual_hits, $background_mutation_rate);

=head1 DESCRIPTION

This package calculates the so-called path-scan statistic
P-value for sets of putative cancer genes under the null
hypothesis that somatic mutations found in data are the result
of a random process characterized by the background mutation
rate.
This test is applied to, for example, a biologically-relevant group of genes,
say all the genes in a particular pathway, for which somatic mutation data are
available.
A low p-value would imply that the null hypothesis should be
rejected.
In other words, the result suggests that the mutation configuration
in this pathway is probably not the result of a strictly random
process.

=head2 Nature of the Path-Scan Test

This statistic considers individual genes in a "binary" fashion,
i.e. a gene is either mutated (has one or more mutations) or it is not
mutated.
I<The number of mutations in a mutated gene is not
considered.>
This is the "path-scan" aspect of the
test.

Why is such information
discarded?
The somatic background mutation rate is typically
very small compared to the size of the average
gene.
Consequently, the expected number of mutations
in any given gene is very low, much less than one, in
fact.
Under the null hypothesis, most genes will have no
mutations.
Genes with one (or just a few) may be interesting, but when many
genes in a biologically-relevant group (say a pathway) have one (or
just a few) mutations, that could be a sign of some underlying I<non-random>
process.
In other words, this test is useful in cases where many
genes in a group might each contribute a small component
(i.e. a small fitness advantage) in the context of the disease
process.
What this test is not concerned with (and will I<not> detect)
is the case where a single, specific gene has a non-random
association and it reflects this fact via a large number of
mutations.
Other single-gene tests should presumably flag such
cases.
The path-scan test should, therefore, be thought
of as just one tool within a larger statistical
"toolbox".

=head2 Assumptions in the Test

The main assumption is that a single background
mutation rate applies to the set of genes of
interest.
That is, the rate does not vary among genes, among
chromosomes (if more than one hosts genes of interest),
etc.

=head1 BUGS AND OPPORTUNITIES FOR EXTENSION

Coefficients are recalculated for every individual test, but it
would be good for these to persist between tests, adding more as
necessary (i.e. if a subsequent test involves more genes than the current
one).

=head1 AUTHOR

Michael C. Wendl

S<mwendl@wustl.edu>

Copyright (C) 2007, 2008 Washington University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

=head1 METHODS

The available methods are as follows.

=cut

################################################################################
##                                                                            ##
##                      P R O G R A M M E R   N O T E S                       ##
##                                                                            ##
################################################################################
#
#  The obj schematic resembles:
#
#  $obj = {
#
#  #__LIST OF GENE LENGTHS (RAW DATA)
#
#     one list-within-a-list implies doing the exact solution, while 2 or more
#     lists-within-a-list imply doing the approximate "binning" approach
#
#     genes => [
#                 [3500, 4234, 5609, 4550, 10763, 9879,...],
#                 [33500, 44234, 35609, 34550, 110763, 49879,...],
#                       :
#     ],
#
#
#  #__CONTEXT: 0=EXACT OR >0 FOR APPROXIMATE
#
#     this indicator is 0 in the exact context, but in the approximate context
#     it gives the number of lists that have been passed
#
#     context_approx => 1,
#
#
#  #__EXPECTATION
#
#     this is the expected number of mutated genes and is only used for the
#     asymptotic (Poisson) model --- it will computed whenever the gene sizes
#     are passed as a single list (i.e. as when we do the exact solution), but
#     it is only used for the Poisson solution
#
#     expected_genes => 0.097,
#
#
#  #__TOTAL COUNT OF GENES (FOR CONVENIENCE)
#
#     this is the total count of all genes, whether in a single list (exact)
#     context, or multiple lists (approximate "binning" context)
#
#     num_genes => 35,
#
#
#  #__AVERAGE GENE LENGTHS AND SIZES OF LISTS
#
#     these will only be defined when passing multiple lists of genes, i.e.
#     for the approximate "binning" context - each element in 'avg_lengths' and
#     'list_sizes' give, respectively, the average length and the number of
#     genes in the corresponding list of gene raw data above
#
#     avg_lengths => [5445.454, 48984.39, ...],
#     list_sizes => [21, 34, 32, 54, 19, ...],
#
#
#  #__TOTAL LENGTH OF ALL GENES CONCATENATED AND SIZE OF MOST POPULOUS LIST
#
#     this is the total length of all genes, whether in a single list (exact)
#     context, or multiple lists (approximate "binning" context), and the number
#     of elements in the most populous list, respectively -- if there is just
#     one list, then the latter will be equivalent to the value in 'num_genes'
#
#     total_length => 75983,
#     most_populous_list = 35,
#
#
#  #__LIST OF CORRESPONDING "MODIFIED" BERNOULLI "PROBABILITY OF FAILURE" VALUES
#
#     in the exact-solution context, each element is the bernoulli probability
#     for the corresponding gene length in the single list-within-a-list (above)
#     but in the approximate-"binning"-solution-context, each value is the
#     bernoulli probability for the average gene length in the corresponding
#     list of genes, e.g. the zero-th value is for the zero-th list, the first
#     value for the first list, etc. - the code will interpret the context
#     automatically
#
#     mpvals => [0.003596, 0.004243, 0.005625, 0.004560, 0.01082, 0.009928,...],
#
#
#  #__BACKGROUND MUTATION RATE I.E. PROBABILITY OF A MUTATION IN A GIVEN NUCLEO
#     mutation_prob => 0.000001,
#
#
#  #__BINOMIAL COEFFICIENTS
#
#     this will only be defined when passing multiple lists of genes, i.e.
#     for the approximate "binning" solution - we only define the symmetric
#     half of pascal's triangle
#
#     binom_coeffs => [[1], [1], [1,2], [1,3], [1,4,6], [1,5,10], ....],
#
#  };

################################################################################
##                                                                            ##
##                         P U B L I C   M E T H O D S                        ##
##                                                                            ##
################################################################################

#  ===
#  NEW   create a new path-scan object
#  ===   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 new

This is the usual object constructor, which
optionally takes all the gene-length data as
input.
If you want to use the exact probability solution or the
asymptotic approximate solution, pass all lengths in a single list
reference

	my $pmobj = PathScan->new ([3434, 54565, 6445, ...]);

but if you want to use the convolution approximation method, divide the list of
gene sizes into the desired number of bins and pass each of these as a
reference

	my $pmobj = PathScan->new ([3434, 54565], [6445, ...]);

In other words, the way you pass these arguments at partially determines
the context in which you will obtain your P-value for this set of
genes.
The latter choice is typically betetr, as it
gives good accuracy and good computational
efficiency.
Conversely, the exact solution is identically correct, but can be difficult to
compute.
The asymptotic approximation is always computationally
efficient, but not necessarily accurate for small test
sets.

=cut

sub new {
   my $class = shift;
   my (@list_of_gene_length_lists) = @_;

#__OBJECT TEMPLATE
   my $self = {};

#__PROCESS GENE LENGTH DATA IF SPECIFIED
   &store_genes ($self, @list_of_gene_length_lists)
      if scalar @list_of_gene_length_lists;

#__BLESS INTO CLASS AND RETURN OBJECT
   bless $self, $class;
   return $self;
}

#  ==========
#  PATH SCAN   path-scan (tail) probability value
#  ==========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 path_scan

This function calculates the path-scan statistic in one of the
appropriate contexts (exact or convolution approximation, as described
above).
It takes the actual number of "hits" you've observed in
the data, i.e. the number of genes that have a mutated
status.

	my $pval = $pmobj->path_scan (7);

If you have not yet done the pre-processing with respect to the
background mutation rate (see below), then pre-processing can be
executed implicitly by passing the rate as the second
argument.

	my $pval = $pmobj->path_scan (7, 0.000001);

=cut

sub path_scan {
   my $obj = shift;
   my ($actual_hits, $mutation_prob) = @_;

#__SOME BASIC CHECKS
   croak "argument list: '$actual_hits' must be an integer"
      unless $actual_hits =~ /^-?\d+$/;
   croak "argument list: '$actual_hits' cant be negative"
      unless $actual_hits >= 0;
   croak "argument list: '$actual_hits' cant exceed number of " .
      "genes '$obj->{'num_genes'}'" if $actual_hits > $obj->{'num_genes'};

#__PREPROCESS IF NECESSARY
   if ($mutation_prob) {
      $obj->preprocess ($mutation_prob);
   }

#__CALCULATE BINOMIAL COEFFICIENTS IF NECESSARY
   if ($obj->{'context_approx'}) {
      $obj->binom_coeffs unless defined $obj->{'binom_coeffs'};
   }

#__CALCULATE P-VAL WITH THE MINIMUM OF EFFORT
   my $half_of_num_genes = $obj->{'num_genes'} / 2;
   my $path_scan_pval = 0;

#__COMPUTE ACTUAL TAIL IF NUM HITS IS IN THIS RANGE
   if ($actual_hits > $half_of_num_genes) {

   #__COMPUTE THE ACTUAL UPPER TAIL
      for (my $k = $actual_hits; $k <= $obj->{'num_genes'}; $k++) {
         my $pval;
         if ($obj->{'context_approx'}) {
            $pval = $obj->p_value_binomial_approx ($k);
         } else {
            $pval = $obj->p_value_exact ($k);
         }
         $path_scan_pval += $pval;
      }

   #__RETURN RESULT
      return $path_scan_pval;

#__COMPUTE 1 - COMPLIMENTARY TAIL NUM HITS IS IN THIS RANGE
   } else {

   #__COMPUTE COMPLIMENTARY (LOWER) TAIL
      for (my $k = 0; $k < $actual_hits; $k++) {
         my $pval;
         if ($obj->{'context_approx'}) {
            $pval = $obj->p_value_binomial_approx ($k);
         } else {
            $pval = $obj->p_value_exact ($k);
         }
         $path_scan_pval += $pval;
      }

   #__RETURN 1 - THIS VAL AS THE RESULT
      return 1 - $path_scan_pval;
   }
}

#  =============
#  CDF_TRUNCATED   truncated cummulative distribution function
#  =============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  NOTE ON INDEX OF p_value_binomial_approx FUNCTION
#
#  we have the identity for the probability masses of
#
#     P(0) + P(1) + P(2) + ... + P(M) = 1
#
#  where M is the total number of genes, so that the tailed p-values in the
#  CDF list, ordered from most extreme to least extreme are
#  
#  P(K >= 0) = 1
#  P(K >= 1) = 1 - P(0)               = P(K >= 0) - P(0)
#  P(K >= 2) = 1 - P(0) - P(1)        = P(K >= 1) - P(1)
#  P(K >= 3) = 1 - P(0) - P(1) - P(2) = P(K >= 2) - P(2)
#       :      :    :      :      :        :         :
#  P(K >= k) =                        = P(K >= k-1) - P(k-1)
#
#  So, the index of p_value_binomial_approx is '$k-1' and each succeeding
#  element is UNSHIFTED onto the list to get the desired ordering from most
#  extreme to least extreme.
#
#  Method "cdf_asymptot" has the same indexing for its probability calling
#  function.

=head2 cdf_truncated

This function returns the cummulative distribution in the
context of the convolution approximation I<truncated> such that
it contains only enough information to process the given number of
hits.

	my $pvals_list = $pmobj->cdf_truncated ($hits);

The list is ordered from most extreme to least extreme probability
tail values, i.e. the last value in the list is always
unity.
However, tailed p-values more extreme than that associated with the argument
are not, in fact, calculated, but rather are replaced with the flag
-1.
This saves processing time and also reduces the chances of numerical
overflow for large pathways, as the full CDF must ultimately raise
an "mval" (>1) to a power equal to the number of genes in the
pathway.
The method assumes you have already done the pre-processing with respect to the
background mutation rate.

=cut

sub cdf_truncated {
   my $obj = shift;
   my ($hits) = @_;

#__MAKE SURE HITS IS AN INTEGER BETWEEN 0 AND THE NUMBER OF GENES
   croak "hit number must be integer > 0" unless $hits =~ /^\d+$/ && $hits >= 0;
   croak "hit number must be less than total number of genes"
      unless $hits <= $obj->{'num_genes'};

#__INCREMENT BY ONE FOR ANY METHODS THAT MAY ALSO NEED THE NEXT MORE EXTREME
#  TAILED P-VALUE, E.G. THE LANCASTER CORRECTION IN THE POPULATION CALCULATION,
#  UNLESS WE'RE ALREADY AT THE TOTAL NUMBER OF GENES, E.G. FOR A 1-GENE PATHWAY
#  LIKE HSA04112
   $hits++ unless $hits == $obj->{'num_genes'};

#__MAKE SURE PRE-PROCESSING HAS ALREADY BEEN DONE
   croak "cannot call 'cdf_truncated' without first pre-processing"
      unless defined $obj->{'mutation_prob'};

#__CALCULATE BINOMIAL COEFFICIENTS IF NECESSARY
   if ($obj->{'context_approx'}) {
      $obj->binom_coeffs unless defined $obj->{'binom_coeffs'};
   }

#__CALCULATE BINOMIAL COEFFICIENTS IF NECESSARY
#  if ($obj->{'context_approx'}) {
#     $obj->binom_coeffs unless defined $obj->{'binom_coeffs'};
#  }

#__CDF INITIALIZED WITH UNITY --- MORE EXTREME VALS WILL BE UNSHIFTED IN FRONT
   my $path_scan_pval = 1;
   my $cdf = [$path_scan_pval];

#__COMPUTE CDF VALS PROGRESSIVELY MORE EXTREME PUSHING EACH TO FRONT OF LIST
   for (my $k = 1; $k <= $obj->{'num_genes'}; $k++) {

   #__COMPUTE ACTUAL TAILED P-VALUE IF WE'RE WITHIN TRUNCATION RANGE
      if ($k <= $hits) {
         $path_scan_pval -= $obj->p_value_binomial_approx ($k-1);
         unshift @{$cdf}, $path_scan_pval;

   #__OTHERWISE JUST INSERT A FLAG TO FUNCTION AS A PLACEHOLDER FOR OTHER
   #  METHODS THAT EXPECT THE *FORM* OF THE LIST TO BE A FULL CDF
      } else {
         unshift @{$cdf}, -1;
      }
   }

#__RETURN TRUNCATED CDF
   return $cdf;
}

#  ===
#  CDF   cummulative distribution function
#  ===   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 cdf

This function returns the cummulative distribution in one of the
appropriate contexts (exact or convolution approximation, as described
above).
There are no arguments,

	my $pvals_list = $pmobj->cdf;

unless you have not yet done the pre-processing with respect to the
background mutation rate (see below), in which case the pre-processing can be
executed implicitly by passing the rate as the sole
argument.

	my $pvals_list = $pmobj->cdf (0.000001);

The list is ordered from most extreme to least extreme probability
tail values, i.e. the last value in the list is always
unity.

=cut

sub cdf {
   my $obj = shift;
   my ($mutation_prob) = @_;
   my $cdf = [];

#__PREPROCESS IF NECESSARY
   if ($mutation_prob) {
      $obj->preprocess ($mutation_prob);
   }

#__CALCULATE BINOMIAL COEFFICIENTS IF NECESSARY
   if ($obj->{'context_approx'}) {
      $obj->binom_coeffs unless defined $obj->{'binom_coeffs'};
   }

#__CALCULATE P-VAL WITH THE MINIMUM OF EFFORT
   my $path_scan_pval = 0;

#__COMPUTE CDF STARTING WITH MOST THE EXTREME STATE WORKING TOWARD LEAST EXTREME
   for (my $k = $obj->{'num_genes'}; $k >= 0; $k--) {
      my $pval;
      if ($obj->{'context_approx'}) {
         $pval = $obj->p_value_binomial_approx ($k);
      } else {
         $pval = $obj->p_value_exact ($k);
      }
      $path_scan_pval += $pval;
      push @{$cdf}, $path_scan_pval;
   }
   return $cdf;
}

#  ============
#  CDF_ASYMPTOT   cummulative distribution function using asymptotic analysis
#  ============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 cdf_asymptot

This function returns the cummulative distribution based on asymptotic
analysis.
There are no arguments, i.e.

	my $pvals_list = $pmobj->cdf_asymptot;

unless you have not yet done the pre-processing with respect to the
background mutation rate (see below), in which case the pre-processing can be
executed implicitly by passing the rate as the sole
argument.

	my $pvals_list = $pmobj->cdf_asymptot (0.000001);

The list is ordered from most extreme to least extreme probability
tail values, i.e. the last value in the list is always
unity.

Note that asymptotic analysis gives a function (the Poisson) having infinite
support.
The infinite tail probability for all values past the most extreme
physical case are all bundled into that most extreme
p-value.

=cut

sub cdf_asymptot {
   my $obj = shift;
   my ($mutation_prob) = @_;

#__CDF INITIALIZED WITH UNITY --- MORE EXTREME VALS WILL BE UNSHIFTED IN FRONT
   my $path_scan_pval = 1;
   my $cdf = [$path_scan_pval];

#__PREPROCESS IF NECESSARY
   if ($mutation_prob) {
      $obj->preprocess ($mutation_prob);
   }

#__COMPUTE CDF VALS PROGRESSIVELY MORE EXTREME PUSHING EACH TO FRONT OF LIST
#  SEE PROGRAMMING NOTES OF CDF_TRUNCATED METHOD THAT EXPLAIN INDEXING OF THE
#  PROBABILITY CALL
   for (my $k = 1; $k <= $obj->{'num_genes'}; $k++) {
      $path_scan_pval -= $obj->p_value_asymptot_approx ($k-1);
      unshift @{$cdf}, $path_scan_pval;
   }
   return $cdf;
}

#  ===================
#  PATH SCAN ASYMPTOT  asymptotic path-scan probability value (CDF)
#  ===================  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 path_scan_asymptot

This function calculates the path-scan statistic in the asymptotic (Poisson)
context.
It takes the actual number of "hits" you've observed in
the data, i.e. the number of genes that have a mutated
status.

	my $pval = $pmobj->path_scan_asymptot (7);

If you have not yet done the pre-processing with respect to the
background mutation rate (see below), then pre-processing can be
executed implicitly by passing the rate as the second
argument.

	my $pval = $pmobj->path_scan_asymptot (7, 0.000001);

You must set up the object, somewhat paradoxically,
I<as if you will be doing the calculation in the exact
context>.
(This is a consequence of how data are stored internally within the
object.)

=cut

#  The Poisson distribution has infinite support, so to be rigorous, we
#  should sum probability masses to a very large number. Strictly speaking,
#  we get the same (rigorous) result by simply computing the complement and
#  subtracting that from one.

sub path_scan_asymptot {
   my $obj = shift;
   my ($actual_hits, $mutation_prob) = @_;

#__SOME BASIC CHECKS
   croak "argument list: '$actual_hits' must be an integer"
      unless $actual_hits =~ /^-?\d+$/;
   croak "argument list: '$actual_hits' cant be negative"
      unless $actual_hits >= 0;
   croak "argument list: '$actual_hits' cant exceed number of " .
      "genes '$obj->{'num_genes'}'" if $actual_hits > $obj->{'num_genes'};

#__PREPROCESS IF NECESSARY
   if ($mutation_prob) {
      $obj->preprocess ($mutation_prob);
   }

#__POISSON HAS INFINITE SUPPORT SO COMPUTE COMPLIMENTARY (LOWER) TAIL AS PVAL
   my $path_scan_pval = 0;
   for (my $k = 0; $k < $actual_hits; $k++) {
      $path_scan_pval += $obj->p_value_asymptot_approx ($k);
   }

#__RETURN 1 - THIS VAL AS THE RESULT
   return 1 - $path_scan_pval;
}

################################################################################
##                                                                            ##
##                   S E M I - P R I V A T E   M E T H O D S                  ##
##                                                                            ##
##  methods that are not ordinarily called externally but can be if needed    ##
##  because they are cast according to the object-oriented interface          ##
##                                                                            ##
################################################################################

=head2 additional methods

The basic functionality of this package is encompassed in the methods described
above.
However, some lower-level functions can also sometimes be
useful.

=cut

#  =============
#  P VALUE EXACT   returns probability of exactly k genes mutated by chance
#  =============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  number of required recursion levels is based on number of hits = k

=head3 p_value_exact

This function returns the exact value of the
probability I<mass> for a specific number of
hits.

	$pval_exact = $pmobj->p_value_exact (7);

You must make sure to call this only if you've
configured the object in the exact context (see
above).

=cut

sub p_value_exact {
   my $obj = shift;
   my ($k) = @_;

#__SOME BASIC CHECKS
   croak "argument list: '$k' must be an integer" unless $k =~ /^-?\d+$/;
   croak "argument list: '$k' cant be negative" unless $k >= 0;
   croak "argument list: '$k' cant exceed number of genes '$obj->{'num_genes'}'"
      if $k > $obj->{'num_genes'};

#__COMPUTE THE LEADING TERM
   my $leading_term = exp (- $obj->{'mutation_prob'} * $obj->{'total_length'});

#__THIS IS ALSO THE EXACT P-VALUE IN THE SPECIFIC CASE OF K=0
   return $leading_term if $k == 0;

#__CALCULATE EXACT P-VALUE FOR K>0
   my $total_subpart_of_pval = $obj->_recursive_exact_pval_calculation ($k, 0, 1, 1);
   return $leading_term * $total_subpart_of_pval;
}

#  =======================
#  P VALUE BINOMIAL APPROX   returns approx prob of k genes mutated by chance
#  =======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  NOTE: we have not implemented the approximate method for the case of just
#  a single bin, although a simple solution exists (path-scan test notes pp 16)
#  because the assumption is that, if the user just passes a single list, they
#  want the exact solution. In other words, we assume that the approximate
#  solution is only desired when there are at least 2 lists.
#
#  number of required recursion levels is based on the number of lists
#  that the data have been binned into

=head3 p_value_binomial_approx

This function returns the convolution approximated
value (i.e. using the binomial binning approximation)
of the probability I<mass> for a specific number of
hits.

	$pval_exact = $pmobj->p_value_binomial_approx (7);

You must make sure to call this only if you've
configured the object in the approximate binomial context (see
above).
Also, you must explicitly calculate the
necessary binomial coefficients beforehand
(see C<binom_coeffs>).

=cut

sub p_value_binomial_approx {
   my $obj = shift;
   my ($k) = @_;

#__SOME BASIC CHECKS
   croak "argument list: '$k' must be an integer" unless $k =~ /^-?\d+$/;
   croak "argument list: '$k' cant be negative" unless $k >= 0;
   croak "argument list: '$k' cant exceed number of genes '$obj->{'num_genes'}'"
      if $k > $obj->{'num_genes'};
   croak "need to first calculate binomial coefficients"
      unless defined $obj->{'binom_coeffs'};

#__COMPUTE THE LEADING TERM
   my $leading_term = exp (- $obj->{'mutation_prob'} * $obj->{'total_length'});

#__THIS IS ALSO THE EXACT P-VALUE IN THE SPECIFIC CASE OF K=0
   return $leading_term if $k == 0;

#__COMPUTE RECURSIVE PORTION OF THE SOLUTION
   my $num_lists = $obj->{'context_approx'};
   my $total_subpart_of_pval =
      $obj->_recursive_binom_pval_calculation ($num_lists - 1, 0, $k);

#  THIS CHUNK DOES NOT WORK IN GENERAL BECAUSE A CDF IS NOT GUARENTEED TO START
#  FROM K=0
#
#
#__CALCULATE EXACT P-VALUE FOR K>0
#  my $prob_mass = $leading_term * $obj->{'binomial_multiplier'} *
#                  $total_subpart_of_pval;
#__TAKE CARE OF THE POWER TERM FOR NEXT ITERATION (STEP-BY-STEP BUILD-UP)
#  -1 MODIFIER RECONCILES MATH LIST NUMBERING WITH PERL LIST NUMBERING
#  $obj->{'binomial_multiplier'} *= $obj->{'mpvals'}->[$num_lists-1];
#__RETURN P
#  return $prob_mass;

#__RETURN (-1 MODIFIER RECONCILES MATH LIST NUMBERING WITH PERL LIST NUMBERING)
   return $leading_term * $obj->{'mpvals'}->[$num_lists-1]**$k *
          $total_subpart_of_pval;
}

#  =======================
#  P VALUE ASYMPTOT APPROX   returns poisson approx prob of k genes mutated
#  =======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head3 p_value_asymptot_approx

This function returns the asymptotic approximated
value (i.e. using the Poisson limit approximation)
of the probability I<mass> for a specific number of
hits.

	$pval_exact = $pmobj->p_value_asymptot_approx (7);

Somewhat paradoxically, you must make sure to call this only
if you've configured the object in the exact context (see
above).

=cut

sub p_value_asymptot_approx {
   my $obj = shift;
   my ($k) = @_;

#__SOME BASIC CHECKS
   croak "argument list: '$k' must be an integer" unless $k =~ /^-?\d+$/;
   croak "argument list: '$k' cant be negative" unless $k >= 0;
   croak "argument list: '$k' cant exceed number of genes '$obj->{'num_genes'}'"
      if $k > $obj->{'num_genes'};

#__SOME SHORTHAND
   my $avg = $obj->{'expected_genes'};

#__ASYMPTOTIC PROBABILITY MASS
   my $pmass;
   if ($k) {

   #__PROBABILITY MASS
      my $ramanuj = log($k * (1 + 4*$k*(1 + 2*$k))) / 6;
      my $arg_h0 = $k * log ($avg / $k) - $avg + $k;
      $arg_h0 -= $ramanuj;
      $arg_h0 -= LOG_PI_OVER_2;
      $pmass = exp ($arg_h0);
   } else {
      $pmass = exp (-$avg);
   }

#__RETURN PROBABILITY MASS
   return $pmass;
}

#  ===========
#  STORE GENES   run some basic consistency checks on raw data then store them
#  ===========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  CONFIGURED FOR BOTH EXACT AND APPROXIMATE CONTEXTS
#
#  this has to take the whole data struct (list-of-lists) because it can be
#  called externally and we don't expect the user to have to do the proper
#  looping

=head3 store_genes

Stores the raw gene length
data.
Use this if you did not pass these data to
C<new> before you call any calculation
methods.
Works in the same way as C<new>, described
above.
Specifically, the context is partially determined by whether you
pass a single list (exact context or asymptotic approximation)

	$pmobj->store_genes ([3434, 54565, 6445, ...]);

or more than one list (convolution approximate context)

	$pmobj->store_genes ([3434, 54565], [6445, ...]);

=cut

sub store_genes {
   my $obj = shift;
   my (@list_of_gene_length_lists) = @_;

#__INIT
   my ($total_length, $most_populous_list) = (0, 0);
   $obj->{'num_genes'} = 0;

#__DETERMINE WHETHER USER WANTS TO RUN IN EXACT CONTEXT OR APPROXIMATE CONTEXT
#  IF APPROXIMATE THE RECORD THE NUMBER OF LISTS
   $obj->{'context_approx'} = 0;
   if (scalar @list_of_gene_length_lists > 1) {
      $obj->{'context_approx'} = scalar @list_of_gene_length_lists;
   }

#__PROCESS THE DATA
   foreach my $gene_list (@list_of_gene_length_lists) {

   #__MAKE SURE THIS IS A LIST
      croak "argument must be list reference" unless ref $gene_list eq "ARRAY";

   #__PROCESS EACH GENE IN THIS LIST
      my $list_concat_length = 0;
      foreach my $length (@{$gene_list}) {

      #__MAKE SURE THIS LENGTH IS A POSITIVE INTEGER
         croak "argument list: '$length' must be an integer"
            unless $length =~ /^-?\d+$/;
         croak "argument list: '$length' must be positive" unless $length > 0;

      #__TALLY THIS LENGTH TO THE TOTAL LENGTH OF THIS LIST
         $list_concat_length += $length;
      }

   #__SOME ADDITIONAL INFORMATION NEEDED IF IN APPROXIMATE CONTEXT
      if ($obj->{'context_approx'}) {

      #__RECORD THE AVERAGE LENGTH IN THIS LIST AND THE LIST SIZE
         my $num_elements = scalar @{$gene_list};
         push @{$obj->{'avg_lengths'}}, $list_concat_length/$num_elements;
         push @{$obj->{'list_sizes'}}, $num_elements;

      #__DISCERN WHICH SUB-LIST HAS THE MOST MEMBERS
         $most_populous_list = $num_elements
            if $num_elements > $most_populous_list;
      }

   #__TALLY THE NUMBER OF GENES IN THIS LIST TO THE GRAND TOTAL
      $obj->{'num_genes'} += scalar @{$gene_list};

   #__SAVE THIS GENE LIST TO THE OBJECT
      push @{$obj->{'genes'}}, $gene_list;

   #__TALLY TO OVERALL LENGTH OF ALL THE GENES
      $total_length += $list_concat_length;
   }

#__RECORD GRAND TOTAL LENGTH OF ALL GENES CONCATENATED TOGETHER
   $obj->{'total_length'} = $total_length;

#__RECORD SIZE OF MOST POPULOUS LIST
   $obj->{'most_populous_list'} = $most_populous_list;

#
#  NOTES
#
#    Here, we should probably now make sure to erase
#       $obj->{'mutation_prob'} and $obj->{'mpvals'} = $mpvals
#    if they exist, because their contents would no longer be valid (i.e. as
#    based on an "old" set of gene lengths)
#
}

#  ============
#  BINOM COEFFS   calculate binomial coefficients
#  ============   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head3 binom_coeffs

Calculates the binomial coefficients needed
in the binomial (convolution) approximate
solution.

	$pmobj->binom_coeffs;

The internal data structure is essentially the
symmetric half of the appropriately-sized Pascal
triangle.
Considerable memory is saved by not storing the full
triangle.

=cut

### This is called automatically, if necessary, before the probability
### calculation, so there is not typically a need to call it
### manually.

sub binom_coeffs {
   my $obj = shift;
   croak "need to know most populous list first"
      unless defined $obj->{'most_populous_list'};
   carp "already have binomial coefficients" if defined $obj->{'binom_coeffs'};

#__SET-UP FOR COEFFICIENTS
   my $bin_coeffs = [1];
   $obj->{'binom_coeffs'} = [];
   push (@{$obj->{'binom_coeffs'}}, $bin_coeffs);

#__CALCULATE COEFFICIENTS UP TO THAT REQUIRED BY THE MOST POPULOUS LIST
   for (my $i = 1; $i <= $obj->{'most_populous_list'}; $i++) {
      $bin_coeffs = &next_bin_coeff_row ($i, $bin_coeffs);
      push (@{$obj->{'binom_coeffs'}}, $bin_coeffs);
   }
}

#  ==========
#  PREPROCESS   calculate list of modified bernoulli probabilities for each gene
#  ==========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  CONFIGURED FOR BOTH EXACT AND APPROXIMATE CONTEXTS

=head3 preprocess

Calculates the Bernoulli kernel probabilities for the individual genes
or gene bins

	$pmobj_binom->preprocess ($background_mutation_rate);

The data structure can be re-configured to run the test with different
background mutation rates by just re-calling this routine with a different
value

	$pmobj_binom->preprocess ($new_background_mutation_rate);

=cut

sub is_float {
    my $val = shift;
    if ($val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
        return 1;
    }
    else {
        return 0;
    }
}

sub preprocess {
   my $obj = shift;
   my ($mutation_prob) = @_;

#__PRELIMINARY VALIDATION
   croak "need background mutation rate" unless $mutation_prob;
   croak "background mutation '$mutation_prob' rate must be a p-val"
      unless (Genome::Model::Tools::Music::PathScan::PathScan::is_float($mutation_prob) &&
      $mutation_prob > 0 && $mutation_prob < 1);
   croak "preprocessing: no data" unless defined $obj->{'genes'};

#__JUST RETURN SILENTLY IF WE'RE ABOUT TO COMPUTE WHAT HAS ALREADY BEEN DONE
#  carp "data already preprocessed" if defined $obj->{'mpvals'};
#
#  NOTE: IF THE CALLER WANTS TO NOW ANALYZE A NEW BACKGROUND MUTATION RATE
#        TEST FOR EQUIVALENCE TO THE OLD ONE USING "STRING" CONTEXT
   if (defined $obj->{'mpvals'} && $obj->{'mutation_prob'} eq $mutation_prob) {
      return;
   }

#__CALCULATE THE MODIFIED BERNOULLI PROBABILITY VALUE FOR EACH GENE
   my $mpvals = [];

#__CALCULATE IN THIS FASHION IF IN APPROXIMATE CONTEXT
   if ($obj->{'context_approx'}) {
      foreach my $length (@{$obj->{'avg_lengths'}}) {
         my $pval = exp ($mutation_prob * $length) - 1;
         if ($pval > 0) {
            push @{$mpvals}, $pval;
         } else {
            croak "unexpected bernoulli pval '$pval' from average gene length '$length'";
         }
      }

   #__ALSO CALCULATE EXPECTED NUMBER OF GENES BEING MUTATED IN CASE WE NEED TO
   #  REVERT TO ASYMPTOTIC ANALYSIS
      my $expected_mutated_genes = 0;
      foreach my $gene_list (@{$obj->{'genes'}}) {
         foreach my $length (@{$gene_list}) {
            $expected_mutated_genes += 1 - exp (- $mutation_prob * $length);
         }
      }
      $obj->{'expected_genes'} = $expected_mutated_genes;

#__ELSE CALCULATE IN THIS FASHION IF IN EXACT CONTEXT OR IF USING THE
#  ASYMPTOTIC (POISSON) APPROXIMATION
   } else {
      my $expected_mutated_genes = 0;

   #__THERE IS ONLY 1 LIST OF GENES HERE (INDEX 0 IN THE DATA STRUCTURE)
      foreach my $length (@{$obj->{'genes'}->[0]}) {

      #__ACTUAL PVAL FOR THIS GENE BEING MUTATED TALLIED TOWARD EXPECTED VALUE
         $expected_mutated_genes += 1 - exp (- $mutation_prob * $length);

      #__ALGORITHMIC MODEL PVALUE (SEE NOTES)
         my $pval = exp ($mutation_prob * $length) - 1;

      #__STORE ALGORITHMIC MODELED PVALUE
         if ($pval > 0) {
            push @{$mpvals}, $pval;
         } else {
            croak "unexpected bernoulli pval '$pval' from gene length '$length'";
         }
      }

   #__STORE EXPECTED NUMBER OF MUTATED GENES
      $obj->{'expected_genes'} = $expected_mutated_genes;
   }

#__STORE RESULT
   $obj->{'mutation_prob'} = $mutation_prob;
   $obj->{'mpvals'} = $mpvals;
}

################################################################################
##                                                                            ##
##            M E T H O D S   M E A N T   T O   B E   P R I V A T E           ##
##                                                                            ##
##  methods that cannot be called in a contextually meaningful way from an    ##
##  external application using the object-oriented interface                  ##
##                                                                            ##
################################################################################

#  ================================
#  RECURSIVE EXACT PVAL CALCULATION   recursion for handling summation of terms
#  ================================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub _recursive_exact_pval_calculation {
   my $obj = shift;
   my ($k, $prev_recursion_level, $start_index, $coeff_for_this_level) = @_;

#__THE CURRENT LEVEL OF RECURSION
   my $curr_recursion_level = $prev_recursion_level + 1;
   my $local_contribution = 0;

#__LOOP AT THE CURRENT RECURSION LEVEL USING APPROPRIATE START/STOP INDECES
   for (my $i = $start_index;
           $i <= $obj->{'num_genes'} - $k + $curr_recursion_level;
           $i++) {

   #__RECURSE FURTHER IF NECESSARY: NOT
      if ($curr_recursion_level < $k) {

      #__RECURSION: NOTE LAST ARG RESOLVES DIFFERENCE BETWEEN PHYSICAL GENE
      #  NUMBERING (STARTING AT 1) VS NUMBERING IN PERL LIST (STARTING AT 0)
         $local_contribution +=
            $obj->_recursive_exact_pval_calculation (
               $k,
               $curr_recursion_level,
               $i+1, # see mcw notes
               $coeff_for_this_level * $obj->{'mpvals'}->[$i-1] # resolve number
            );

   #__ELSE WE'RE "AT THE BOTTOM" SO TAKE THE NECESSARY PRODUCT: NOTE AGAIN
   #  WE RESOLVE DIFFERENCE BETWEEN PHYSICAL GENE NUMBERING (STARTING AT 1)
   #  VS NUMBERING IN PERL LIST (STARTING AT 0)
      } else {
         $local_contribution +=
                  $coeff_for_this_level * $obj->{'mpvals'}->[$i-1];
      }
   }

#__RETURN RESULT TO THE ANTECEDENT LEVEL
   return $local_contribution;
}

#  ================================
#  RECURSIVE BINOM PVAL CALCULATION   recursion for handling summation of terms
#  ================================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  HERE $obj->{'mpvals'}->[$i-1] IS THE BERNOULLI PROBABILITY FOR THE i-th BIN
#  I.E. THE -1 MODIFIER RECONCILES MATH LIST NUMBERING WITH PERL LIST NUMBERING
#
#  note that the summation structure in the solution means that binomial
#  coefficients that are out-of-bounds can nevertheless be called upon, so
#  we must explicitly test when this is the case and skip these (see code below)

sub _recursive_binom_pval_calculation {
   my $obj = shift;
   my ($recurs_limit, $prev_recursion_level, $k_stop_index) = @_;

#__THE CURRENT LEVEL OF RECURSION
   my $curr_recursion_level = $prev_recursion_level + 1;
   my $local_contribution = 0;

#__SUBSCRIPT
   my $current_list_number = $recurs_limit + 1 - $curr_recursion_level;

#__EFFICIENCY -- WE DO NOT NEED TO COMPUTE THE POWER TERM WITHIN THE LOOP
#  THE POWER VALUE CORRESPONDS TO THE LOOP ITERATOR, SO WE CAN INCREMENT
#  THESE TOGETHER -- SEE SOLUTION NOTES
   my $fraction_coeff = 1;
   my $multiplier = $obj->{'mpvals'}->[$current_list_number-1] /
                    $obj->{'mpvals'}->[$current_list_number];
 
#__LOOP AT THE CURRENT RECURSION LEVEL USING APPROPRIATE START/STOP INDECES
   for (my $k = 0;
           $k <= $k_stop_index;
           $k++) {

   #__RECURSE FURTHER IF NECESSARY: NOT
      if ($curr_recursion_level < $recurs_limit) {

      #__RECURSION: CURRENT ITERATION IN THIS LOOP IS STOP INDEX FOR NEXT
         if (
            $obj->{'list_sizes'}->[$current_list_number] >= $k_stop_index - $k
         ) {
            $local_contribution += $fraction_coeff
               * &bincoeff ($obj, $obj->{'list_sizes'}->[$current_list_number], $k_stop_index - $k)
               * $obj->_recursive_binom_pval_calculation (
                  $recurs_limit,
                  $curr_recursion_level,
                  $k
               );
         }

   #__ELSE WE'RE "AT THE BOTTOM" SO TAKE THE NECESSARY PRODUCT: NOTE AGAIN
   #  WE RESOLVE DIFFERENCE BETWEEN PHYSICAL GENE NUMBERING (STARTING AT 1)
   #  VS NUMBERING IN PERL LIST (STARTING AT 0)

      } else {
         if (
            $obj->{'list_sizes'}->[$current_list_number-1] >= $k &&
            $obj->{'list_sizes'}->[$current_list_number] >= $k_stop_index - $k
         ) {
            $local_contribution += $fraction_coeff
            * &bincoeff ($obj, $obj->{'list_sizes'}->[$current_list_number-1], $k)
            * &bincoeff ($obj, $obj->{'list_sizes'}->[$current_list_number], $k_stop_index - $k);
         }
      }

   #__TAKE CARE OF THE POWER TERM FOR NEXT ITERATION (STEP-BY-STEP BUILD-UP)
      $fraction_coeff *= $multiplier;
   }

#__RETURN RESULT TO THE ANTECEDENT LEVEL
   return $local_contribution;
}

#  ==================
#  NEXT BIN COEFF ROW   compute 1/2 row i of binomial coefficients given row i-1
#  ==================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub next_bin_coeff_row {
   my ($i, $im1_row) = @_;

#__FIRST ELEMENT I,0 IN ALL ROWS DEFINED AS UNITY
   my $next_row = [1];

#__IF I > 1 COMPUTE REST OF NEXT (I) ROW USING PASCAL TRIANGLE ON PREV (I-1) ROW
   if ($i > 1) {

   #__COMPUTE STOPPING POINT BASED ON SYMMETRY
      my $i_mirror = $i / 2;
      my $i_end = int $i_mirror;

   #__FILL IN INTERMEDIATE ELEMENTS I,1 TO I,I_END USING PASCALS TRIANGLE
      for (my $iposition = 1; $iposition <= $i_end; $iposition++) {

      #__START WITH "LEFT SIDE" VAL OF PREVIOUS ROW
         my $element = $im1_row->[$iposition-1];

      #__COMPUTE NEXT ELEMENT USING THE PASCAL TRIANGLE METHOD
         if ($iposition == $i_mirror) {
            $element += $im1_row->[$iposition-1];
         } else {
            $element += $im1_row->[$iposition];
         }

      #__SAVE
         push (@{$next_row}, $element);
      }
   }

#__RETURN LIST OF HALF SYMETRIC NEXT ROW OF BINOMIAL COEFFICIENTS
   return $next_row;
}

# ========
# BINCOEFF   return binomial coefficient using 1/2 symetric stored table
# ========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  RETURNS C_{$i_top, $k_bottom}

sub bincoeff {
   my $obj = shift;
   my ($i_top, $k_bottom) = @_;

#__SYMMETRY CUTOFF IS HALF THE VALUE OF THE I'TH ROW
   my $i_half = $i_top / 2;

#__IF I,K IS WITHIN THE STORED SYMMETRIC HALF OF TRIANGLE THEN SIMPLY RETURN VAL
   if ($k_bottom <= $i_half) {
      return $obj->{'binom_coeffs'}->[$i_top]->[$k_bottom];

#__ELSE COMPUTE SYMMETRIC REFLECTION OF NON-STORED COMPONENT AND RETURN THAT VAL
   } else {
      my $k_reflect = $i_top - $k_bottom;
      return $obj->{'binom_coeffs'}->[$i_top]->[$k_reflect];
   }
}

################################################################################
##                                                                            ##
##             T R A I L I N G   P O D   D O C U M E N T A T I O N            ##
##                                                                            ##
################################################################################

=head1 EXAMPLES

The following examples may be helpful in using this
package.
In each case, assume we have first executed some required preliminary
code.

	#__USE THE PACKAGE
	   use PathScan;

	#__SOME DATA FOR AN "EXACT CONTEXT" CALCULATION
	   my $genes_exact = [
	      4000, 4000, 4000, 4000, 4000,
	      15000, 15000, 15000, 15000, 15000,
	      35000, 35000, 35000, 35000, 35000
	   ];

	#__SOME DATA FOR AN "APPROXIMATE CONTEXT" CALCULATION
	   my @genes_binned = (
	      [4000, 4000, 4000, 4000, 4000],
	      [15000, 15000, 15000, 15000, 15000],
	      [35000, 35000, 35000, 35000, 35000]
	   );

=head2 simple path-scan test

Here, we compare the values returned by both the exact
and approximate algorithms over the whole domain of
possible hits for a case where the answers should be
identical.

	#__SET BACKGROUND MUTATION RATE
	   my $rho = 0.00002;

	#__CONFIGURE OBJECTS IN "EXACT" AND "APPROXIMATE" CONTEXTS
	   my $pmobj_exact = PathScan->new ($genes_exact);
	   $pmobj_exact->preprocess ($rho);

	   my $pmobj_binom = PathScan->new (@genes_binned);
	   $pmobj_binom->preprocess ($rho);

	#__CALCULATE AND TALLY THE MAXIMUM DIFFERENCE
	   my $maxdiff = 0;
	   for (my $i = 0; $i <= scalar @{$genes_exact}; $i++) {
	      my $pm_pval_exact = $pmobj_exact->path_scan($i);
	      my $pm_pval_binom = $pmobj_binom->path_scan($i);
	      my $diff = abs ($pm_pval_exact - $pm_pval_binom);
	      $maxdiff = $diff if $diff > $maxdiff;
	      print "$i hits: $pm_pval_exact   $pm_pval_binom   $diff\n";
	   }
	   print "MAXIMUM DIFFERENCE IS $maxdiff\n";

=head2 testing at different background rates

This example shows how to run the test for a fixed number of hits, say 7 in this
case, for various different background mutation rates.

	#__CONFIGURE OBJECT
	   my $pmobj_binom = PathScan->new (@genes_binned);

	#__CALCULATE
	   for (my $rho = 0.00001; $rho <= 0.0001; $rho += 0.00001) {
	      my $pm_pval_binom = $pmobj_binom->path_scan(7, $rho);
	      print "7 hits at background $rho : P = $pm_pval_binom\n";
	   }

Note that we did not run C<preprocess> explicitly,
but rather let the C<path_scan> method call it implicitly
for each new value of the background mutation
rate.

=head2 computing asymptotic approximate solution

The asymptotic (Poisson) approximate probabiltiy value is straightforward to
compute.

	#__SET BACKGROUND MUTATION RATE
	   my $rho = 0.00002;

	#__CONFIGURE OBJECT
	   my $pmobj_poisson = PathScan->new ($genes_exact);
	   $pmobj_poisson->preprocess ($rho);

	#__P-VALUE FOR 7 OBSERVED MUTATED GENES
	   my $pm_pval_poisson = $pmobj_exact->path_scan_asymptot (7);

=head2 accessing individual probability masses

The probability masses for specific numbers of hits can also be
calculated.


	#__SET BACKGROUND MUTATION RATE
	   my $rho = 0.00002;

	#__CONFIGURE OBJECT
	   my $pmobj_exact = PathScan->new ($genes_exact);
	   $pmobj_exact->preprocess ($rho);

	#__CALCULATE MASSES
	   my $total_prob = 0;
	   for (my $i = 0; $i <= scalar @{$genes_exact}; $i++) {
	      my $pval_exact = $pmobj_exact->p_value_exact ($i);
	      $total_prob += $pval_exact;
	      print "$i hits : probability mass =  $pval_exact\n";
	   }
	   print "total probability = $total_prob\n";

=cut

################################################################################
##                                                                            ##
##                                -  E N D  -                                 ##
##                                                                            ##
################################################################################

1;
