package Genome::Model::Tools::Music::PathScan::CombinePvals;

#__STANDARD PERL PACKAGES
use strict;
use Carp;
use Statistics::Distributions;

#__CONSTANT OF PI -- NEEDED IN RAMANUJAN APPROX FOR POISSON PROBABILITY MASSES
#  (SEE "PATHSCAN TEST" NOTES PP 29-31)
#  use constant PI => 4*atan2 1, 1;
#  use constant LOG_PI_OVER_2 => log (PI) / 2;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

CombinePvals - combining probabilities from independent tests of
significance into a single aggregate figure

=head1 SYNOPSIS

	use CombinePvals;

	my $obj = CombinePvals->new ($reference_to_list_of_pvals);

	my $pval = $obj->method_name;

	my $pval = $obj->method_name (@arguments);

=head1 DESCRIPTION

There are a variety of circumstances under which one might have a number
of different kinds of tests and/or separate instances of the same kind
of test for one particular null hypothesis, where each of these tests returns a
p-value.
The problem is how to properly condense this list of
probabilities into a single value so as to be able to make
a statistical inference, e.g. whether to reject the null
hypothesis.
This problem was examined heavily starting about the 1930s, during
which time numerous mathematical contintencies were treated, e.g. dependence
vs. independence of tests, optimality, inter-test weighting, computational
efficiency, continuous vs. discrete tests and combinations thereof,
etc.
There is quite a large mathematical literature on
this topic (see L</"REFERENCES"> below) and any one
particular situation might incur some of the above
subtleties.
This package concentrates on some of the more straightforward
scenarios, furnishing various methods for combining
p-vals.
The main consideration will usually be the trade-off between
the exactness of the p-value (according to strict frequentist
modeling) and the computational efficiency, or even its actual
feasibility.
Tests should be chosen with this factor in
mind.

Note also that this scenario of combining p-values (many
tests of a single hypothesis) is fundamentally different
from that where a given hypothesis is tested multiple
times.
The latter instance usually calls for some method of multiple testing
correction.

=head1 REFERENCES

Here is an abbreviated list of the substantive works on the topic of combining
probabilities.

=over

=item *

Birnbaum, A. (1954)
I<Combining Independent Tests of Significance>,
Journal of the American Statistical Association B<49>(267), 559-574.

=item *

David, F. N. and Johnson, N. L. (1950)
I<The Probability Integral Transformation When the Variable is Discontinuous>,
Biometrika B<37>(1/2), 42-49.

=item *

Fisher, R. A. (1958)
I<Statistical Methods for Research Workers>, 13-th Ed. Revised,
Hafner Publishing Co., New York.

=item *

Lancaster, H. O. (1949)
I<The Combination of Probabilities Arising from Data in Discrete Distributions>,
Biometrika B<36>(3/4), 370-382.

=item *

Littell, R. C. and Folks, J. L. (1971)
I<Asymptotic Optimality of Fisher's Method of Combining Independent Tests>,
Journal of the American Statistical Association B<66>(336), 802-806.

=item *

Pearson, E. S. (1938)
I<The Probability Integral Transformation for Testing Goodness of Fit and
Combining Independent Tests of Significance>,
Biometrika B<30>(12), 134-148.

=item *

Pearson, E. S. (1950)
I<On Questions Raised by the Combination of Tests Based on Discontonuous
Distributions>,
Biometrika B<37>(3/4), 383-398.

=item *

Pearson, K. (1933)
I<On a Method of Determining Whether a Sample Of Size N Supposed to
Have Been Drawn From a Parent Population Having a Known Probability
Integral Has Probably Been Drawn at Random>
Biometrika B<25>(3/4), 379-410.

=item *

Van Valen, L. (1964)
I<Combining the Probabilities from Significance Tests>,
Nature B<201>(4919), 642.

=item *

Wallis, W. A. (1942)
I<Compounding Probabilities from Independent Significance Tests>,
Econometrica B<10>(3/4), 229-248.

=item *

Zelen, M. and Joel, L. S. (1959)
I<The Weighted Compounding of Two Independent Significance Tests>,
Annals of Mathematical Statistics B<30>(4), 885-895.

=back

=head1 AUTHOR

Michael C. Wendl

S<mwendl@wustl.edu>

Copyright (C) 2009 Washington University

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

=head1 GENERAL REMARKS ON METHODS

The available methods are listed
below.
Each of computational techniques assumes that tests, as
well as their associated p-values, are independent of one
another and none considers any form of differential
weighting.

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
#  #__PROBABILITY VALUES FOR THE INDIVIDUAL TESTS
#     pvals => [0.103, 0.078, 0.03, 0.2,...],
#
#  #__PRODUCT OF THE INDIVIDUAL PROBABILITY VALUES
#     big_q => 0.103 * 0.078 * 0.03 * 0.2 * ...,
#
#  #__THE ACTUAL NUMBER OF TESTS (SAVED FOR CONVENIENCE)
#     num_tests = integer,
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

=head1 CONSTRUCTOR METHODS

These methods return an object in the CombinePvals
class.

=cut

################################
#  BEGIN: CONSTRUCTOR METHODS  #
################################

#  ===
#  NEW   create a new object
#  ===   ~~~~~~~~~~~~~~~~~~~

=head2 new

This is the usual object constructor, which takes a mandatory,
but otherwise un-ordered (reference to a) list of the p-values
obtained by a set of independent
tests.

	my $obj = CombinePvals->new ([0.103, 0.078, 0.03, 0.2,...]);

The method checks to make sure that all elements are actual
p-values, i.e. they are real numbers and they have values bounded by 0 and
1.

=cut

sub new {
   my $class = shift;
   my ($pvals) = @_;

#__OBJECT TEMPLATE
   my $self = {};

#__PROCESS PVALS IF THEY'RE SPECIFIED
   if (defined $pvals && $pvals) {

   #__MAKE SURE THIS IS A LIST
      croak "argument must be list reference" unless ref $pvals eq "ARRAY";

   #__SAVE LIST
      $self->{'pvals'} = $pvals;

   #__PROCESS THE INPUT
      my ($big_q, $num_tests) = (1, 0);
      foreach my $pval (@{$pvals}) {

      #__MAKE SURE THIS IS A PVAL
         croak "'$pval' is not a p-val" unless &is_a_pval ($pval);

      #__COUNT NUMBER OF TESTS
         $num_tests++;

      #__SAVE THE PRODUCT --- THIS IS THE DERIVED TEST STATISTIC
         $big_q *= $pval;
      }
      $self->{'big_q'} = $big_q;
      $self->{'num_tests'} = $num_tests;

#__OTHERWISE CROAK
   } else {
      croak "must specify a list of pvals as an argument";
   }

#__BLESS INTO CLASS AND RETURN OBJECT
   bless $self, $class;
   return $self;
}

##############################
#  END: CONSTRUCTOR METHODS  #
##############################

=head1 EXACT ENUMERATIVE PROCEDURES FOR STRICTLY DISCRETE DISTRIBUTIONS

When all the individual p-vals are derived from tests based on discrete
distributions, the "standard" continuum methods cannot be used in the strictest
sense.
Both Wallis (1942) and Lancaster (1949) discuss the option of
full enumeration, which will only be feasible when there are a
limited number of p-values and their range is not too
large.
Feasibility experiments are suggested, depending
upon the type of hardware and size of
calculation.

=cut

# Again, these methods do some rudimentary checking, but the calling
# program is responsible for making sure all elements are actual
# p-values, i.e. real numbers, have values bounded by 0 and 1,
# etc.
# They are also responsible for making sure all p-values are
# listed in decreasing order of extremity, as illustrated
# below.

####################################################################
#  BEGIN: EXACT ENUMERATIVE PROCEDURES FOR DISCRETE DISTRIBUTIONS  #
####################################################################

#  ====================
#  EXACT ENUM ARBITRARY
#  ====================
#
#  exact enumerative solution for a set of p-values obtained from an
#  arbitrary set of not-necessarily-the-same *discrete* distributions

=head2 exact_enum_arbitrary

This routine is designed for combining p-values
from completely arbitrary discrete probability
distributions.
It takes a list-of-lists data structure, each list being the probability
tails I<ordered from most extreme to least extreme> (i.e. as a probability
cummulative density function) associated with each individual
test.
However, the ordering of the lists themselves is not
important.
For instance, Wallis (1942) gives the example of two binomials, a one-tailed
test having tail values of 0.0625, 0.3125, 0.6875, 0.9375, and 1, and a
two-tailed test having tail values 0.125, 0.625, and
1.
We would then call this method using

	my $pval = $obj->exact_enum_arbitrary (
	   [0.0625, 0.3125, 0.6875, 0.9375, 1],
	   [0.125, 0.625, 1]
	);

The internal computational method is relatively
straightforard and described in detail by
Wallis (1942).
Note that this method does "all-by-all" multiplication,
so it is the least efficient, although entirely
exact.

=cut

sub exact_enum_arbitrary {
   my $obj = shift;
   my (@pvals_lists) = @_;
   my $pval = 0;

#__NUMBER OF LISTS SHOULD BE SAME AS NUMBER OF PVALS PASSED TO CONSTRUCTOR
   my $num_lists = scalar @pvals_lists;
   croak "number of lists passed ($num_lists) not equal to number of tests" .
         "in 'new' constructor ($obj->{'num_tests'})"
         unless $num_lists == $obj->{'num_tests'};

#__CHECK LIST INPUT
   my $list_num = 0;
   foreach my $list (@pvals_lists) {
      $list_num++;
      my $previous_pval = 0;
      foreach my $test_pval (@{$list}) {

      #__MAKE SURE THIS IS A PVAL
         croak "'$test_pval' in distribution $list_num is not a p-val"
            unless &is_a_pval ($test_pval);

      #__MAKE SURE THIS PVAL IS LARGER THAN PREVIOUS ONE: I.E. THIS IS A C.D.F.
         croak "distribution $list_num is not in ascending order (not a CDF)"
            unless $test_pval > $previous_pval;
 
      #__RESET
         $previous_pval = $test_pval;
      }
   }

#__COMBINE INDIVIDUAL-TEST-P-VALS INTO A SINGLE P-VAL
   $pval = $obj->_recursive_exact_enum_arbitrary ([@pvals_lists], 0, 1, 1);

#__RETURN P-VAL
   return $pval;
}

#  ====================
#  EXACT ENUM IDENTICAL
#  ====================
#
#  exact enumerative solution for a set of p-values obtained from one,
#  or rather, a set of identical *discrete* distributions

=head2 exact_enum_identical

This routine is designed for combining a set of p-values that
all come from a single probability
distribution.

	NOT IMPLEMENTED YET

=cut

##################################################################
#  END: EXACT ENUMERATIVE PROCEDURES FOR DISCRETE DISTRIBUTIONS  #
##################################################################

=head1 TRANSFORMS FOR CONTINUOUS DISTRIBUTIONS

The mathematical literature furnishes several straightforward
options for combining p-vals if all of the distributions
underlying all of the individual tests are
continuous.

=cut

####################################################
#  BEGIN: TRANSFORMS FOR CONTINUOUS DISTRIBUTIONS  #
####################################################

#  ===========================
#  FISHER CHI-SQUARE TRANSFORM
#  ===========================
#
#  Fisher's solution using the chi-square transform, valid strictly for
#  continuum distributions, but can be used approximately for discrete
#  distributions. Accuracy increases with the support of the distributions.

=head2 fisher_chisq_transform

This routine implements R.A. Fisher's (1958, originally 1932) chi-square
transform method for combining p-vals from continuous
distributions, which is essentially a CPU-efficient approximation of
K. Pearson's log-based result (see e.g. Wallis (1942) pp
232).
Note that the underlying distributions are
not actually relevant, so no arguments are
passed.

	my $pval = $obj->fisher_chisq_transform;

This is certainly the fastest and easiest method for combining p-vals,
but its accuracy for discrete distributions will not usually be very
good.
For such cases, an exact or a corrected method are better
choices.

=cut

sub fisher_chisq_transform {
   my $obj = shift;

#__GO THROUGH LIST OF INDIVIDUAL-TEST-P-VALS ACCUMULATING FISHER'S LOG
#  TRANSFORM TO CHI-SQUARE STATISTIC
#  my $chisq = 0;
#  foreach my $test_pval (@{$obj->{'pvals'}}) {
####  $chisq += - 2 * log ($test_pval);
#     $chisq -= 2 * log ($test_pval);
#  }

#__WALLIS (1942) MAKES THIS CLEVER SIMPLIFICATION
   my $chisq = -2 * log ($obj->{'big_q'});

#__TRANSFORM: TWICE THE DEGREES OF FREEDOM
   my $dof = 2 * $obj->{'num_tests'};

#__NOW GET P-VAL FROM A CHI-SQUARE TEST
   my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);

#__RETURN P-VAL
   return $pval;
}

##################################################
#  END: TRANSFORMS FOR CONTINUOUS DISTRIBUTIONS  #
##################################################

=head1 CORRECTION PROCEDURES FOR DISCRETE DISTRIBUTIONS: LANCASTER'S MODELS

Enumerative procedures quickly become infeasible if the
number of tests and/or the support of each test grow
large.
A number of procedures have been described for
correcting the methodologies designed for continuum
testing, mostly in the context of applying so-called continuity
corrections.
Essentially, these seek to "spread" dicrete data out into a pseudo-continuous
configuration as appropriate as possible, and then apply standard
transforms.
Accuracy varies and should be suitably established in each
case.

The methods in this section are due to H.O. Lancaster (1949), who discussed
two corrections based upon the idea of describing how a chi-square transformed
statistic varies between the points of a discrete
distribution.
Unfortunately, these methods require one to pass some extra information
to the routines, i.e. not only the CDF (the p-val of each test), but the
CDF value associated with the next-most-extreme
statistic.
These two pieces of information are the basis of
interpolating.
For example, if an underlying distribution has the possible tail values of
0.0625, 0.3125, 0.6875, 0.9375, 1 and the test itself has a value of
0.6875, then you would pass I<both> 0.3125 I<and> 0.6875 to the
routine.
I<In all cases, the lower value, i.e. the more
extreme one, precedes higher value in the argument
list.>
While there generally will be some extra inconvenience in obtaining
this information, the accuracy is much improved over Fisher's
method.

=cut

#  PROGRAMMING NOTE ON LANCASTER'S METHODS
#
#  Each method differs substantively by only a few lines of code, so there
#  are a lot of extra lines here that are required to offer the user 3
#  individually-named methods. This should be fixed when time permits, for
#  example, perhaps pass the name of the correction method as an argument too.

################################################################################
#  BEGIN: CORRECTION PROCEDURES FOR DISCRETE DISTRIBUTIONS: LANCASTER'S MODELS #
################################################################################

#  ==========================================================
#  LANCASTER'S MEAN-CONTINUITY-CORRECTED CHI-SQUARE TRANSFORM
#  ==========================================================

=head2 lancaster_mean_corrected_transform

This method is based on the mean value of the chi-squared transformed
statistic.

	my $pval = $obj->lancaster_mean_corrected_transform (@cdf_pairs);

Its accuracy is good, but the method is not strictly defined if one
of the tests has either the most extreme or second-to-most-extreme
statistic.

=cut

sub lancaster_mean_corrected_transform {
   my $obj = shift;
   my (@fxm1_and_fx_pvals) = @_;

#__NUMBER OF CDF PAIRS SHOULD BE SAME AS NUMBER OF PVALS PASSED TO CONSTRUCTOR
   my $num_lists = scalar @fxm1_and_fx_pvals;
   croak "number of pairs passed ($num_lists) not equal to number of tests" .
         "in 'new' constructor ($obj->{'num_tests'})"
         unless $num_lists == $obj->{'num_tests'};

#__ACCUMULATE LANCASTER'S MEAN CHI-SQUARED STATISTIC
   my ($chisq, $list_num) = (0, 0);
   foreach my $fxm1_and_fx_pair (@fxm1_and_fx_pvals) {
      $list_num++;
      my ($fxm1, $fx) = @{$fxm1_and_fx_pair};

   #__MAKE SURE BOTH ARE PVALS
      croak "'$fxm1' in cdf pair $list_num is not a p-val"
         unless &is_a_pval ($fxm1);
      croak "'$fx' in cdf pair $list_num is not a p-val"
         unless &is_a_pval ($fx);

   #__MAKE SURE THEY'RE ORDERED AS EXPECTED
      croak "cdf pair $list_num is not in ascending order"
         unless $fx > $fxm1;

   #__MEAN CORECTION
      $chisq += 2 * (1 - ($fx * log($fx) - $fxm1 * log($fxm1))/($fx - $fxm1));
   }

#__TRANSFORM: TWICE THE DEGREES OF FREEDOM
   my $dof = 2 * $obj->{'num_tests'};

#__NOW GET P-VAL FROM A CHI-SQUARE TEST
   my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);

#__RETURN P-VAL
   return $pval;
}

#  ============================================================
#  LANCASTER'S MEDIAN-CONTINUITY-CORRECTED CHI-SQUARE TRANSFORM
#  ============================================================

=head2 lancaster_median_corrected_transform

This method is based on the median value of the chi-squared transformed
statistic.

	my $pval = $obj->lancaster_median_corrected_transform (@cdf_pairs);

Its accuracy may sometimes be not quite as good as when using the
average, but the method is strictly defined for I<all> values of the
statistic.

=cut

sub lancaster_median_corrected_transform {
   my $obj = shift;
   my (@fxm1_and_fx_pvals) = @_;

#__NUMBER OF CDF PAIRS SHOULD BE SAME AS NUMBER OF PVALS PASSED TO CONSTRUCTOR
   my $num_lists = scalar @fxm1_and_fx_pvals;
   croak "number of pairs passed ($num_lists) not equal to number of tests" .
         "in 'new' constructor ($obj->{'num_tests'})"
         unless $num_lists == $obj->{'num_tests'};

#__ACCUMULATE LANCASTER'S MEAN CHI-SQUARED STATISTIC
   my ($chisq, $list_num) = (0, 0);
   foreach my $fxm1_and_fx_pair (@fxm1_and_fx_pvals) {
      $list_num++;
      my ($fxm1, $fx) = @{$fxm1_and_fx_pair};

   #__MAKE SURE BOTH ARE PVALS
      croak "'$fxm1' in cdf pair $list_num is not a p-val"
         unless &is_a_pval ($fxm1);
      croak "'$fx' in cdf pair $list_num is not a p-val"
         unless &is_a_pval ($fx);

   #__MAKE SURE THEY'RE ORDERED AS EXPECTED
      croak "cdf pair $list_num is not in ascending order"
         unless $fx > $fxm1;

   #__MEDIAN CORRECTION
      if ($fxm1) {
         $chisq -= 2 * log (($fx + $fxm1)/2);
      } else {
         $chisq += 2 * (1 - log ($fx));
      }
   }

#__TRANSFORM: TWICE THE DEGREES OF FREEDOM
   my $dof = 2 * $obj->{'num_tests'};

#__NOW GET P-VAL FROM A CHI-SQUARE TEST
   my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);

#__RETURN P-VAL
   return $pval;
}

#  ===========================================================
#  LANCASTER'S MIXED-CONTINUITY-CORRECTED CHI-SQUARE TRANSFORM
#  ===========================================================

=head2 lancaster_mixed_corrected_transform

This method is a mixture of both the mean and median
methods.
Specifically, mean correction is used wherever it
is well-defined, otherwise median correction is
used.

	my $pval = $obj->lancaster_mixed_corrected_transform (@cdf_pairs);

This will be a good way to handle certain
cases.

=cut

sub lancaster_mixed_corrected_transform {
   my $obj = shift;
   my (@fxm1_and_fx_pvals) = @_;

#__NUMBER OF CDF PAIRS SHOULD BE SAME AS NUMBER OF PVALS PASSED TO CONSTRUCTOR
   my $num_lists = scalar @fxm1_and_fx_pvals;
   croak "number of pairs passed ($num_lists) not equal to number of tests" .
         "in 'new' constructor ($obj->{'num_tests'})"
         unless $num_lists == $obj->{'num_tests'};

#__ACCUMULATE LANCASTER'S MEAN CHI-SQUARED STATISTIC
   my ($chisq, $list_num) = (0, 0);
   foreach my $fxm1_and_fx_pair (@fxm1_and_fx_pvals) {
      $list_num++;
      my ($fxm1, $fx) = @{$fxm1_and_fx_pair};

   #__MAKE SURE BOTH ARE PVALS
      croak "'$fxm1' in cdf pair $list_num is not a p-val"
         unless &is_a_pval ($fxm1);
      croak "'$fx' in cdf pair $list_num is not a p-val"
         unless &is_a_pval ($fx);

   #__MAKE SURE THEY'RE ORDERED AS EXPECTED
   #
   #  NOTE: WE ALLOW FOR EQUIVALENCE OF ADJACENT VALUES (WITHIN FLOATING-POINT
   #        PRECISION) FOR THOSE CASES WHERE THE CDF IS LOCALLY EXTREMELY
   #        ASYMPTOTIC

      croak "cdf pair $list_num is not in ascending order"
         if $fx < $fxm1;
   #     unless $fx > $fxm1;

   #__NO CORRECTION NEEDED IF VALS IDENTICAL WITHIN FLOATING-POINT PRECISION
      if ($fx == $fxm1) {
         $chisq -= 2 * log ($fx); # Fisher

   #__ELSE APPLY LANCASTER'S CONTINUITY CORRECTION (MIXED DEPENDING UPON VALS)
      } else {

      #__USE LANCASTER'S MEAN IF POSSIBLE
         if ($fxm1) {
            $chisq += 2 * (1 - ($fx * log($fx) - $fxm1 * log($fxm1))/($fx - $fxm1));
      #__OTHERWISE USE LANCASTER'S MEDIAN
         } else {
            $chisq += 2 * (1 - log ($fx));
         }
      }
   }

#__TRANSFORM: TWICE THE DEGREES OF FREEDOM
   my $dof = 2 * $obj->{'num_tests'};

#__NOW GET P-VAL FROM A CHI-SQUARE TEST
   my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);

#__RETURN P-VAL
   return $pval;
}

###############################################################################
#  END: CORRECTION PROCEDURES FOR DISCRETE DISTRIBUTIONS: LANCASTER'S MODELS  #
###############################################################################

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

#  ======================
#  EXACT_ENUM_ARBITRARY_2   2-distrib precursor of exact_enum_arbitrary
#  ======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head3 exact_enum_arbitrary_2

Hard-wired precursor of I<exact_enum_arbitrary> for 2
distributions.
Does no pre-checking, but may be useful for
comparing to the output of the general
program.

=cut

sub exact_enum_arbitrary_2 {
   my $obj = shift;
   my (@pvals_lists) = @_;
   my $pval = 0;

#__NUMBER OF LISTS SHOULD BE SAME AS NUMBER OF PVALS PASSED TO CONSTRUCTOR
   my $num_lists = scalar @pvals_lists;
   croak "number of lists passed ($num_lists) not equal to number of tests" .
         "in 'new' constructor ($obj->{'num_tests'})"
         unless $num_lists == $obj->{'num_tests'};
   croak "you must have exactly 2 tests" unless $obj->{'num_tests'} == 2;

#__TWO-LIST SPECIAL CASE
   my $list1 = $pvals_lists[0];
   my $list2 = $pvals_lists[1];

#__TRAVERSE LIST 1
   for (my $i = 0; $i <= $#{$list1}; $i++) {

   #__TAIL VALUE AND THE PROBABILITY OF THIS TAIL VALUE
      my $ptail1 = $list1->[$i];
      my $probability_of_ptail1;
      if ($i > 0) {
         $probability_of_ptail1 = $list1->[$i] - $list1->[$i-1];
      } else {
         $probability_of_ptail1 = $list1->[$i];
      }

   #__TRAVERSE LIST 2
      for (my $j = 0; $j <= $#{$list2}; $j++) {

      #__TAIL VALUE AND THE PROBABILITY OF THIS TAIL VALUE
         my $ptail2 = $list2->[$j];
         my $probability_of_ptail2;
         if ($j > 0) {
            $probability_of_ptail2 = $list2->[$j] - $list2->[$j-1];
         } else {
            $probability_of_ptail2 = $list2->[$j];
         }

      #__PRODUCT OF TAIL VALUES AND THE PROBABILITY OF THIS PRODUCT
         my $product_tail_pval = $ptail1 * $ptail2;
         my $probability_of_product_tail_pval = $probability_of_ptail1 * $probability_of_ptail2;

      #__TALLY TO RESULTANT COMPOUND P-VAL IF SIGNIFICANT
         $pval += $probability_of_product_tail_pval
            if $product_tail_pval <= $obj->{'big_q'};
      }
   }

#__RETURN P-VAL
   return $pval;
}

#  ======================
#  EXACT_ENUM_ARBITRARY_3   3-distrib precursor of exact_enum_arbitrary
#  ======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head3 exact_enum_arbitrary_3

Hard-wired precursor of I<exact_enum_arbitrary> for 3
distributions.
Does no pre-checking, but may be useful for
comparing to the output of the general
program.

=cut

sub exact_enum_arbitrary_3 {
   my $obj = shift;
   my (@pvals_lists) = @_;
   my $pval = 0;

#__NUMBER OF LISTS SHOULD BE SAME AS NUMBER OF PVALS PASSED TO CONSTRUCTOR
   my $num_lists = scalar @pvals_lists;
   croak "number of lists passed ($num_lists) not equal to number of tests" .
         "in 'new' constructor ($obj->{'num_tests'})"
         unless $num_lists == $obj->{'num_tests'};
   croak "you must have exactly 3 tests" unless $obj->{'num_tests'} == 3;

#__THREE-LIST SPECIAL CASE
   my $list1 = $pvals_lists[0];
   my $list2 = $pvals_lists[1];
   my $list3 = $pvals_lists[2];

#__TRAVERSE LIST 1
   for (my $i = 0; $i <= $#{$list1}; $i++) {

   #__TAIL VALUE AND THE PROBABILITY OF THIS TAIL VALUE
      my $ptail1 = $list1->[$i];
      my $probability_of_ptail1;
      if ($i > 0) {
         $probability_of_ptail1 = $list1->[$i] - $list1->[$i-1];
      } else {
         $probability_of_ptail1 = $list1->[$i];
      }

   #__TRAVERSE LIST 2
      for (my $j = 0; $j <= $#{$list2}; $j++) {

      #__TAIL VALUE AND THE PROBABILITY OF THIS TAIL VALUE
         my $ptail2 = $list2->[$j];
         my $probability_of_ptail2;
         if ($j > 0) {
            $probability_of_ptail2 = $list2->[$j] - $list2->[$j-1];
         } else {
            $probability_of_ptail2 = $list2->[$j];
         }

      #__TRAVERSE LIST 3
         for (my $k = 0; $k <= $#{$list3}; $k++) {

         #__TAIL VALUE AND THE PROBABILITY OF THIS TAIL VALUE
            my $ptail3 = $list3->[$k];
            my $probability_of_ptail3;
            if ($k > 0) {
               $probability_of_ptail3 = $list3->[$k] - $list3->[$k-1];
            } else {
               $probability_of_ptail3 = $list3->[$k];
            }

         #__PRODUCT OF TAIL VALUES AND THE PROBABILITY OF THIS PRODUCT
            my $product_tail_pval = $ptail1 * $ptail2 * $ptail3;
            my $probability_of_product_tail_pval = $probability_of_ptail1 * $probability_of_ptail2 * $probability_of_ptail3;

         #__TALLY TO RESULTANT COMPOUND P-VAL IF SIGNIFICANT
            $pval += $probability_of_product_tail_pval
               if $product_tail_pval <= $obj->{'big_q'};
         }
      }
   }

#__RETURN P-VAL
   return $pval;
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

################################################################################
##                                                                            ##
##            M E T H O D S   M E A N T   T O   B E   P R I V A T E           ##
##                                                                            ##
##  methods that cannot be called in a contextually meaningful way from an    ##
##  external application using the object-oriented interface                  ##
##                                                                            ##
################################################################################

#  ==========================================================================
#  ROUTINE FOR DETERMINING WHETHER A VARIABLE REPRESENTS A LEGITIMATE P-VALUE
#  ==========================================================================

sub is_a_pval {
   my ($val) = @_;
# print "VAL IS '$val'\n";

#__MUST BE A FLOAT (REGEXP: PERL COOKBOOK CHAP 2.1) & MUST BE BOUNDED BY 0 AND 1
   if ($val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/
      && $val >= 0 && $val <= 1) {
      return 1;
##################
#    if ($val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
#       print "VAL '$val' IS REAL\n";
#       if ($val >= 0) {
#          print "  VAL '$val' >= 0\n";
#          if ($val <= 1) {
#             print "    VAL '$val' <= 1\n";
#             return 1;
#          } else {
#             print "    VAL '$val' NOT <= 1\n";
# # chop $val;
# # if ($val <= 1) {
# #   print "      VAL '$val' NOW <= 1\n";
# # } else {
# #   print "      VAL '$val' STILL NOT <= 1\n";
# # }
#             return 0;
#          }
#       } else {
#          print "  VAL '$val' NOT >= 0\n";
#          return 0;
#       }
##################

#__ELSE IT IS NOT A PVAL
   } else {
#     print "VAL '$val' IS NOT REAL\n";
      return 0;
   }
}

#  ================================================================
#  RECURSIVE EXACT PVAL CALCULATION FOR ARITRARY TEST DISTRIBUTIONS
#  ================================================================

sub _recursive_exact_enum_arbitrary {
   my $obj = shift;
   my ($list_of_pvals_lists, $prev_recursion_level,
       $product_tail_pval, $probability_of_product_tail_pval) = @_;

#__THE CURRENT LIST
   my $current_list = $prev_recursion_level;
   my $list = $list_of_pvals_lists->[$current_list];

#__THE CURRENT LEVEL OF RECURSION
   my $curr_recursion_level = $prev_recursion_level + 1;
   my $local_pval = 0;

#__LOOP AT THE CURRENT RECURSION LEVEL
   for (my $i = 0; $i <= $#{$list}; $i++) {

   #__TAIL VALUE AND THE PROBABILITY OF THIS TAIL VALUE
      my $ptail = $list->[$i];
      my $probability_of_ptail;
      if ($i > 0) {
         $probability_of_ptail = $list->[$i] - $list->[$i-1];
      } else {
         $probability_of_ptail = $list->[$i];
      }

   #__RECURSE FURTHER IF NECESSARY
      if ($curr_recursion_level < $obj->{'num_tests'}) {

      #__RECURSE AND ACCUMULATE P-VAL CONTRIBUTIONS
         $local_pval +=
            $obj->_recursive_exact_enum_arbitrary (
               $list_of_pvals_lists,
               $curr_recursion_level,
               $product_tail_pval * $ptail,
               $probability_of_product_tail_pval * $probability_of_ptail
            );

   #__ELSE WE'RE "AT THE BOTTOM" SO TAKE THE NECESSARY PRODUCT
      } else {

      #__FINAL PRODUCTS
         my $local_product_tail_pval = $product_tail_pval * $ptail;
         my $local_probability_of_product_tail_pval =
	    $probability_of_product_tail_pval * $probability_of_ptail;

      #__TALLY IF CONDITION IS SATISFIED
         if ($local_product_tail_pval <= $obj->{'big_q'}) {
            $local_pval += $local_probability_of_product_tail_pval;

      #__ELSE SKIP REST OF THIS DISTRIBUTION CUZ SUCCEEDING VALS ARE ALL LARGER
         } else { 
            last; 
         }
      }
   }

#__RETURN RESULT TO THE ANTECEDENT LEVEL
   return $local_pval;
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

################################################################################
##                                                                            ##
##                                -  E N D  -                                 ##
##                                                                            ##
################################################################################

1;
