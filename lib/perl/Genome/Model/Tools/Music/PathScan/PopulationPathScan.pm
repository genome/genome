package Genome::Model::Tools::Music::PathScan::PopulationPathScan;

#__STANDARD PERL PACKAGES
use strict;
use warnings;

use Carp;
use Genome::Model::Tools::Music::PathScan::CombinePvals;
use Genome::Model::Tools::Music::PathScan::PathScan;
# DEBUG
# print "USING LOCAL MCW VERSION OF POPULATION PATHSCAN\n";
# DEBUG

# DEBUG -- PLEASE REMOVE
# use lib "/home/mwendl/work/perl_modules";
# use PostData;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

PopulationPathScan - apply PathScan test to populations rather
than just single individuals

=head1 SYNOPSIS

	use PopulationPathScan;

	my $obj = PopulationPathScan->new ($ref_to_list_of_gene_lengths);

	$obj->assign ($number_of_compartments);

	$obj->preprocess ($background_mutation_rate);

	$pval = $obj->population_pval_approx ($ref_to_list_of_hits_per_sample);
	$pval = $obj->population_pval_exact ($ref_to_list_of_hits_per_sample);

=head1 DESCRIPTION

The C<PathScan> package is implemented strictly
as a test of a set of genes, e.g. a pathway, for a I<single>
individual.
Specifically, knowing the gene lengths in the pathway, the number
of genes that have at least one mutation, and the estimated background
mutation rate, one can test the null hypothesis that these observed
mutations are well-explained simply by the mechanism of random background
mutation.
However, it will often be the case that data for a pathway will be available
for many individuals, meaning that we now have many tests of the given (single)
hypothesis.
(This should not be confused with the scenario of multiple hypothesis
testing.)
The set of values contains much more information than a single value,
suggesting that significance must be judged on the basis of the collective
result.
For example, while no single p-value by itself may exceed the chosen statistical
threshold, the overall set of probabilities may still give the impression of
significance.
Properly combining such numbers is a necessary, but not entirely trivial
task.
This package basically serves as a high-level interface to first
perform individual tests using the methods of C<PathScan>,
and then to properly combine the resulting p-values using the methods of
C<CombinePvals>.

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

=head1 METHODS

The available methods are listed
below.

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
#  #__GENE LENGTHS IN THE POPULATION PATHSCAN TEST
#     gene_lengths => [474, 1038, 285, ...],
#
#  #__THE ACTUAL NUMBER OF GENES IN TEST (SAVED FOR CONVENIENCE)
#     num_genes = 15,
#
#  #__ARGUMENT LIST FOR PATHSCAN COMPUTATION (PathScan) STRUCTURE
#     IS DETERMINED BY THE WAY THE "assign" METHOD IS CALLED
#     path_scan_arg_list = [],
#     path_scan_arg_list = [ [], [], [] ],
#
#  #__ASSIGN LEVEL (ESSENTIALLY THE ARGUMENT OF THE 'ASSIGN' METHOD)
#     assign_level = 1,
#
#  #__CUMMULATIVE DISTRIBUTION FOR THIS SET OF GENES ORDERED MOST EXTREME
#     TO LEAST EXTREME --- COULD BE EITHER THE "COMPLETE" CDF, I.E. THE ENTIRE
#     DISTRIBUTION
#     cdf = [0.003, 0.0234, 0.1001, 0.23, 0.4, 0.8, 0.9, 0.94, 0.97, 0.99, 1],
#
#  #__OR COULD BE A TRUNCATED LIST WITH JUST ENOUGH VALUES TO DO A CALCULATION
#     I.E. WHERE THE MORE EXTREME TAILED PROBABILITY VALUES ARE OMITTED, BEING
#     REPLACED BY A SIMPLE PLACEHOLDER FLAG -1
#     cdf = [-1, -1, -1, -1, -1, -1, -1, 0.94, 0.97, 0.99, 1],
#
#  #__MAXIMUM NUMBER OF MUTATED GENES TAKEN OVER ALL SAMPLES -- SEE PREPROCESS
#     max_hits = 5,
#  };

################################################################################
##                                                                            ##
##                         P U B L I C   M E T H O D S                        ##
##                                                                            ##
################################################################################

#  ===
#  NEW   create a new population path-scan object
#  ===   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 new

The object constructor takes a mandatory, but otherwise un-ordered reference
to a list of gene lengths comprising the biological group (e.g. a pathway)
whose mutation significance is to be analyzed using the PathScan
paradigm.

	my $obj = PopulationPathScan->new ([474, 1038, 285, ...]);

The method checks to make sure that all elements
are legitimate lengths, i.e. integers exceeding
3.

=cut

sub new {
   my $class = shift;
   my ($gene_lengths) = @_;

#__OBJECT TEMPLATE
   my $self = {};

#__PROCESS GENE LENGTHS IF THEY'RE SPECIFIED
   if (defined $gene_lengths && $gene_lengths) {

   #__MAKE SURE THIS IS A LIST
      croak "argument must be list reference"
         unless ref $gene_lengths eq "ARRAY";

   #__SAVE LIST
      $self->{'gene_lengths'} = $gene_lengths;
      $self->{'num_genes'} = scalar @{$gene_lengths};

   #__VALIDATE THE INPUT
      foreach my $gene_length (@{$gene_lengths}) {

      #__MAKE SURE THIS IS A LEGITIMATE LENGTH
         croak "'$gene_length' is not a gene length" unless
            $gene_length =~ /^\d+$/ && $gene_length >= 3;

      }

#__OTHERWISE CROAK
   } else {
      croak "must specify a list of gene lengths as an argument";
   }

#__BLESS INTO CLASS AND RETURN OBJECT
   bless $self, $class;
   return $self;
}

#  ======
#  ASSIGN   assign the manner in which genes will be internally organized
#  ======   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 assign

This method assigns the manner in which genes will be internally
organized for passing to the PathScan calculation
component.
The main consideration here is how the list may be compartmentalized for greater
computational efficiency, though at some loss of accuracy, for the PathScan
calculation.
If the gene list is long, exact calculation is generally
infeasible.
The method takes a single argument representing the number of compartments
(or sub-lists) the lengths will be divided into, e.g. 1 represents a
single list, i.e. exact computation, 2 indicates two lists, 3 three lists,
etc.

	$obj->assign (3);

The values are then organized internally such that the smallest
genes are grouped together, then the slightly larger ones, and so
forth.
Generally, 3 or 4 lists give reasonable balance
between accuracy and computation (Wendl et al., in
progress).

=cut

#  THIS HAS NOT BEEN IMPLEMENTED YET
#
#  The method can also be called without an argument
#  
#	$obj->assign;
#  
#  in which case the gene lengths will put into a number of
#  compartments such that each one has a maximum of 10
#  values.

sub assign {
   my $obj = shift;
   my ($assign_level) = @_;

#__ORDER THE LIST OF VALIDATED GENE LENGTHS ACCORDING TO INCREASING SIZE
   @{$obj->{'gene_lengths'}} = sort _numerical_ @{$obj->{'gene_lengths'}};
   sub _numerical_ {$a <=> $b}

#__ASSIGN TO A SPECIFIC NUMBER OF COMPARTMENTS IF SPECIFIED
   if (defined $assign_level && $assign_level) {
      $obj->{'assign_level'} = $assign_level;

   #__QUICK-PROCESSING IF THERE'S NO COMPARTMENTALIZATION
      if ($assign_level == 1) {
         $obj->{'path_scan_arg_list'} = $obj->{'gene_lengths'};
         return;
      }

   #__REMAINDER AFTER DIVIDING GENE LIST INTO AN INTEGER-NUMBER OF COMPARTMENTS
      my $remain = $obj->{'num_genes'} % $assign_level;

   #__LENGTH OF ALL COMPARTMENTS (EXCEPT LAST ONE IF THERE'S A REMAINDER)
      my $list_length = ($obj->{'num_genes'} - $remain) / $assign_level;

   #__BUILD-UP THE COMPARTMENTALIZED ARGUMENT LIST
      my ($list_number, $gene_number, $compartment) = (1, 1, []);
      foreach my $gene_length (@{$obj->{'gene_lengths'}}) {
         push @{$compartment}, $gene_length;
         $gene_number++;
         if ($gene_number > $list_length) {
            $list_number++;
            push @{$obj->{'path_scan_arg_list'}}, $compartment;
            ($gene_number, $compartment) = (1, []);
            $list_length += $remain if $list_number == $assign_level;
         }
      }

#__ELSE ASSIGN SUCH THAT EACH COMPARTMENT HAS A MAXIMUM SIZE
   } else {
      croak "illegal assignment level";
   }
}

#  ==========
#  PREPROCESS   set-up PathScan calculation and compute CDF
#  ==========   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head2 preprocess

This method pre-processes the population-level calculation, specifically,
it sets up and executes the PathScan module to obtain the CDF associated
with the given gene set and background mutation
rate.
It takes the latter as an
argument.

	$obj->preprocess (0.0000027);

Executing this method will take various amounts of CPU time,
depending upon the level of accuracy and the number of genes in the
calculation.

The method optionally takes the list of the number of mutated genes
in the group for each sample as a second argument, if this information is
known at this point

	$obj->preprocess (0.0000027, [4, 5, 7, 3, 0, ...]);

and it is usually better to use this form because the
internals will compute only a truncated CDF that is just
sufficient to process this list, rather than computing the full
CDF.
Not only is speed improved, but this helps avoid overflow errors for large
pathways.

=cut

sub preprocess {
   my $obj = shift;
   my ($mutation_prob, $list_of_hits) = @_;
   my $max_hits = 0;

#__PRELIMINARY VALIDATION OF ARGUMENT
   croak "need background mutation rate" unless $mutation_prob;
   croak "background mutation '$mutation_prob' rate must be a p-val"
      unless is_a_pval ($mutation_prob);

#__INVOKE NEW PATHSCAN OBJECT USING PRE-COMPUTED ARGUMENT LIST FOR EXACT SOLN
   my $pm_obj;
   if ($obj->{'assign_level'} == 1) {
      $pm_obj = Genome::Model::Tools::Music::PathScan::PathScan->new ($obj->{'path_scan_arg_list'});

#__OR FOR APPROXIMATE SOLUTION
   } elsif ($obj->{'assign_level'} > 1) {
      $pm_obj = Genome::Model::Tools::Music::PathScan::PathScan->new (@{$obj->{'path_scan_arg_list'}});

   #__ALSO FIND THE MAX NUMBER OF MUTATED GENES AMONG ALL SAMPLES IF GIVEN HITS
      if (defined $list_of_hits && $list_of_hits) {

      #__MAKE SURE THIS IS A LIST
         croak "argument must be list reference"
            unless ref $list_of_hits eq "ARRAY";

      #__HARD-SET MAX HITS TO 1 IN CASE 0 SAMPLES HAVE HITS & TRIGGERS TRUNC CDF
         $max_hits = 1;

      #__FIND MAXIMUM NUMBER OF HITS
         foreach my $hits (@{$list_of_hits}) {

         #__MAKE SURE THIS IS A LEGITIMATE HIT NUMBER
            croak "'$hits' is not a hit number" unless
               $hits =~ /^\d+$/ && $hits >= 0 && $hits <= $obj->{'num_genes'};

         #__RECORD MAXIMUM
            $max_hits = $hits if $hits > $max_hits;
         }

      #__SAVE MAX HITS TO OBJECT
         $obj->{'max_hits'} = $max_hits;
      }

#__ELSE WE CANT PROCESS
   } else {
      croak "I dont understand the 'assign' level you used previously";
   }

#__STANDARD PREPROCESSING FOR PathScan OBJECT
   $pm_obj->preprocess ($mutation_prob);

#__COMPUTE AND STORE CDF -- EITHER FULL OR TRUNCATED
   if ($max_hits) {
      $obj->{'cdf'} = $pm_obj->cdf_truncated ($max_hits);
   } else {
      $obj->{'cdf'} = $pm_obj->cdf;
   }

#  $obj->{'cdf'} = $pm_obj->cdf_asymptot;

#__WEIRD HEURISTIC: MAKE SURE LAST VALUE IN CDF LIST IS ALWAYS IDENTICALLY UNITY
#  
#  We have seen sometimes that the last value appears to be unity, i.e. it
#  prints as such, but the 'is_a_pval' rejects it either on the real number
#  regexp, or the <= 1 condition. Here is an actual croak:
#
#  VAL '1' IS NOT REAL
#  '1' in distribution 1 is not a p-val at Statistics/PopulationPathScan.pm line 395
#
#  Please track down this problem when you have a chance, but this practical
#  fix seems to work acceptably for the moment.

   $obj->{'cdf'}->[$#{$obj->{'cdf'}}] = 1;
}

#  =====================
#  POPULATION_PVAL_EXACT   tail prob for the population using exact enumeration
#  =====================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  PROGRAMMING NOTES:
#
#  1. remember that the cdf list returned from PathScan->cdf
#     is ordered from most etreme (lowest p-value, highest number of hits) to
#     least extreme (highest p-value = 1, lowest number of hits = 0). Therefore,
#     the correct p-value corresponding to the actual number of hits cannot
#     be naively looked-up in the list according to order, but must rather be
#     looked up according to the *reverse order*. For example, for the usual
#     binomial (0.5 + 0.5)^4 (see e.g. Wallis (1942) pp 244), we have
#  
#     cdf      = [0.0625,    0.3125,         0.6875,     0.9375,     1]
#     position = [0,         1,              2,          3,          4]
#     meaning  = [all 4 hit, at least 3 hit, at least 2, at least 1, at least 0]
#
#     therefore, the actual "hit" pvalue is in position
#
#        $obj->{'num_genes'} - $hits
#
#  2. CombinePvals does not yet have a method that exploits
#     scenarios, such as this one, where each individual p-val comes from the
#     _same_ distribution. Currently, we must call "exact_enum_arbitrary",
#     which does a full enumeration. Change methods here if the CombinePvals
#     class ever gets such a method.

=head2 population_pval_exact

This method performs the population-level calculation using exact
enumeration.
It takes the list of the number of mutated genes in the
group for each sample, e.g. each patient's whole genome
sequence, for example

	patient 1: 4 genes in the pathway are mutated
	patient 2: 5 genes in the pathway are mutated
	patient 3: 7 genes in the pathway are mutated
	patient 4: 3 genes in the pathway are mutated
	patient 5: 0 genes in the pathway are mutated
	  :     :  :   :    :  :     :     :     :

which is invoked as

	$pval = $obj->population_pval_exact ([4, 5, 7, 3, 0, ...]);

Most scenarios will not actually be able to make use of this method
because enumeration of all possible cases is rarely computationally
feasible.
This method will mostly be useful for examining small test
cases.

=cut

sub population_pval_exact {
   my $obj = shift;
   my ($list_of_hits) = @_;

#__PROCESS HITS IF THEY'RE SPECIFIED
   if (defined $list_of_hits && $list_of_hits) {

   #__MAKE SURE THIS IS A LIST
      croak "argument must be list reference"
         unless ref $list_of_hits eq "ARRAY";

   #__WE NEED 2 LISTS FOR EXACT METHOD
      my ($default_arg_list, $cdf_list) = ([], []);

   #__VALIDATE AND PROCESS THE INPUT INTO ARGUMENT LISTS
      foreach my $hits (@{$list_of_hits}) {

      #__MAKE SURE THIS IS A LEGITIMATE HIT NUMBER
         croak "'$hits' is not a hit number" unless
            $hits =~ /^\d+$/ && $hits >= 0 && $hits <= $obj->{'num_genes'};

      #__TAIL PVAL FOR THIS HIT NUMBER (SEE PROGRAMMING NOTE ABOVE)
         my $pval_x = $obj->{'cdf'}->[$obj->{'num_genes'} - $hits];

      #__STORE IN DEFAULT CombinePvals ARG LIST
         push @{$default_arg_list}, $pval_x;

      #__CDFS GO IN SPECIAL ARG LIST FOR EXACT ENUMERATION
         push @{$cdf_list}, $obj->{'cdf'};
      }

   #__INVOKE NEW COMBINE_PVALS OBJECT USING PRE-COMPUTED ARGUMENT LIST
      my $combine_obj = Genome::Model::Tools::Music::PathScan::CombinePvals->new ($default_arg_list);

   #__COMPUTE OVERALL "GROUP" P-VALUE BASED ON EXACT ENUMERATION
###### DEBUG
#    print "from PopulationPathScan --- args for new\n";
#    &PostData ($default_arg_list);
#    print "from PopulationPathScan --- args for exact_enum_arbitrary\n";
#    &PostData ($cdf_list);
###### DEBUG
      my $pval = $combine_obj->exact_enum_arbitrary (@{$cdf_list});
      return $pval;

#__OTHERWISE CROAK
   } else {
      croak "must specify a list of number of genes mutated for the sample set";
   }
}

#  ======================
#  POPULATION_PVAL_APPROX   tail prob for the population using Lancaster approx
#  ======================   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  PROGRAMMING NOTE:
#
#  remember that the cdf list returned from PathScan->cdf
#  is ordered from most etreme (lowest p-value, highest number of hits) to
#  least extreme (highest p-value = 1, lowest number of hits = 0). Therefore,
#  the correct p-value corresponding to the actual number of hits cannot
#  be naively looked-up in the list according to order, but must rather be
#  looked up according to the *reverse order*. For example, for the usual
#  binomial (0.5 + 0.5)^4 (see e.g. Wallis (1942) pp 244), we have
#  
#  cdf      = [0.0625,    0.3125,         0.6875,     0.9375,     1]
#  position = [0,         1,              2,          3,          4]
#  meaning  = [all 4 hit, at least 3 hit, at least 2, at least 1, at least 0]
#
#  therefore, the actual "hit" pvalue is in position
#
#     $obj->{'num_genes'} - $hits
#
#  and the next-most-extreme (lower p-value) is in position
#
#     $obj->{'num_genes'} - $hits - 1
#
#  for using Lancaster's correction methods

=head2 population_pval_approx

This method performs the population-level
calculation using Lancaster's approximate transform
correction.
It takes, as a mandatory argument, the list of the number of mutated
genes in the group for each sample, e.g. each patient's whole genome
sequence.

	$pval = $obj->population_pval_approx ([4, 5, 7, 3, 0, ...]);

You must pass the list of hits, even if you already
passed this list earlier to the pre-processing
method.
Most cases will use this method because exact combination
of individual probability values is rarely computationally
feasible.
Note that Lancaster's method typically gives much better (more
accurate) results than Fisher's "standard" chi-square
transform.

=over

=item *

Fisher, R. A. (1958)
I<Statistical Methods for Research Workers>, 13-th Ed. Revised,
Hafner Publishing Co., New York.

=item *

Lancaster, H. O. (1949)
I<The Combination of Probabilities Arising from Data in Discrete Distributions>,
Biometrika B<36>(3/4), 370-382.

=back

=cut

sub population_pval_approx {
   my $obj = shift;
   my ($list_of_hits) = @_;

#__PROCESS HITS IF THEY'RE SPECIFIED
   if (defined $list_of_hits && $list_of_hits) {

   #__MAKE SURE THIS IS A LIST
      croak "argument must be list reference"
         unless ref $list_of_hits eq "ARRAY";

   #__WE NEED 2 LISTS FOR LANCASTER'S METHOD
      my ($default_arg_list, $lancaster_list) = ([], []);

   #__VALIDATE AND PROCESS THE INPUT INTO ARGUMENT LISTS
# DEBUG
#print "processing list of hits\n";
# DEBUG
      foreach my $hits (@{$list_of_hits}) {

      #__MAKE SURE THIS IS A LEGITIMATE HIT NUMBER
         croak "'$hits' is not a hit number" unless
            $hits =~ /^\d+$/ && $hits >= 0 && $hits <= $obj->{'num_genes'};

      #__TAIL PVALS FOR THIS HIT NUMBER (SEE PROGRAMMING NOTE ABOVE)
         my $pval_x = $obj->{'cdf'}->[$obj->{'num_genes'} - $hits];
#        my $pval_x_m_1 = $obj->{'cdf'}->[$obj->{'num_genes'} - $hits - 1];
         my $pval_x_m_1;
         my $x_m_1_index = $obj->{'num_genes'} - $hits - 1;
         if ($x_m_1_index >= 0) {
            $pval_x_m_1 = $obj->{'cdf'}->[$x_m_1_index];
         } else {
            $pval_x_m_1 = 0; # dont allow this to inadvertently loop to list end
         }

# DEBUG
# print "  hit number = $hits\n";
# print "     number of genes = $obj->{'num_genes'}\n";
# print "     pval_x = $pval_x\n";
# print "     pval_x_m_1 = $pval_x_m_1\n";
# DEBUG

      #__STORE IN DEFAULT CombinePvals ARG LIST
      #  (THIS IS ACUTALLY JUST A FORMALITY IF USING LANCASTERS METHOD)
         push @{$default_arg_list}, $pval_x;

      #__STORE IN SPECIAL ARG LIST FOR LANCASTERS METHOD: P(X-1) THEN P(X)
         push @{$lancaster_list}, [$pval_x_m_1, $pval_x];
      }

   #__INVOKE NEW COMBINE_PVALS OBJECT USING PRE-COMPUTED ARGUMENT LIST
      my $combine_obj = Genome::Model::Tools::Music::PathScan::CombinePvals->new ($default_arg_list);

   #__COMPUTE OVERALL "GROUP" P-VALUE BASED ON LANCASTERS TRANSFORM CORRECTION
      my $pval = $combine_obj->lancaster_mixed_corrected_transform (@{$lancaster_list});
      return $pval;

#__OTHERWISE CROAK
   } else {
      croak "must specify a list of number of genes mutated for the sample set";
   }
}

################################################################################
##                                                                            ##
##                       P R I V A T E   M E T H O D S                        ##
##                                                                            ##
################################################################################

#  ==========================================================================
#  ROUTINE FOR DETERMINING WHETHER A VARIABLE REPRESENTS A LEGITIMATE P-VALUE
#  ==========================================================================

sub is_a_pval {

    my ($val) = @_;
$DB::single = 1;
#__MUST BE A FLOAT (REGEXP: PERL COOKBOOK CHAP 2.1) & MUST BE BOUNDED BY 0 AND 1
   if (Genome::Model::Tools::Music::PathScan::PathScan::is_float($val)
      && $val >= 0 && $val <= 1) {
      return 1;

#__ELSE IT IS NOT A PVAL
   } else {
      return 0;
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
