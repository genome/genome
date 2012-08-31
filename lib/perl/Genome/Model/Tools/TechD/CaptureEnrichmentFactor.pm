package Genome::Model::Tools::TechD::CaptureEnrichmentFactor;

use strict;
use warnings;

use Carp;
use Genome;

# Todd Wylie  <twylie@genome.wustl.edu>
# Sun Oct 31 22:15:24 CDT 2010

# -----------------------------------------------------------------------------
# C A P T U R E   E N R I C H M E N T   F A C T O R
# -----------------------------------------------------------------------------

# This class essentially wraps up calculations for capture Enrichment Factor
# calculations into a single call. There are 5 mandatory pieces of information
# which must be passed to the class constructor. All of these values are
# currently available from the current capture coverage pipeline
# incarnation. Mandatory values are:
#
#   capture_unique_bp_on_target
#   capture_duplicate_bp_on_target
#   capture_total_bp
#   genome_total_bp
#   target_total_bp
#

# The following public accessors are provided downstream of the "new" method.
#
#   theoretical_max_enrichment_factor()
#   unique_on_target_enrichment_factor()
#   total_on_target_enrichment_factor()
#   capture_unique_bp_on_target()
#   capture_duplicate_bp_on_target()
#   capture_total_bp()
#   genome_total_bp()
#   target_total_bp()
#
# -----------------------------------------------------------------------------

class Genome::Model::Tools::TechD::CaptureEnrichmentFactor {
    is  => ['Command'],
    has => {
        #input
        capture_unique_bp_on_target => {
            is       => 'Number',
            is_input => 1
        },
        capture_duplicate_bp_on_target => {
            is       => 'Number',
            is_input => 1
        },
        capture_total_bp => {
            is       => 'Number',
            is_input => 1
        },
        genome_total_bp => {
            is       => 'Number',
            is_input => 1
        },
        target_total_bp => {
            is  => 'Number',
            is_input => 1
        },

        # calculated
        theoretical_max_enrichment_factor => {
            is        => 'Number',
            calculate => q| $self->_theoretical_max_enrichment_factor(); |
        },
        unique_on_target_enrichment_factor => {
            is        => 'Number',
            calculate => q| $self->_unique_on_target_enrichment_factor(); |
        },
        total_on_target_enrichment_factor => {
            is        => 'Number',
            calculate => q| $self->_total_on_target_enrichment_factor(); |
        },
        capture_total_on_target_bp => {
            is        => 'Number',
            calculate => q| $self->_capture_total_on_target_bp(); |
        },
    },
};


sub execute {
    my $self = shift;

    # print "\nCaptureEnrichmentFactor input:\n";
    # print "capture_unique_bp_on_target: " . $self->capture_unique_bp_on_target . "\n";
    # print "capture_duplicate_bp_on_target: " . $self->capture_duplicate_bp_on_target . "\n";
    # print "capture_total_bp: " . $self->capture_total_bp . "\n";
    # print "target_total_bp: " . $self->target_total_bp . "\n";
    # print "genome_total_bp: " . $self->genome_total_bp . "\n";

    # print "\nCaptureEnrichmentFactor output:\n";
    # print "theoretical_max_enrichment_factor: " . $self->theoretical_max_enrichment_factor . "\n";
    # print "unique_on_target_enrichment_factor: " . $self->unique_on_target_enrichment_factor . "\n";
    # print "total_on_target_enrichment_factor: " . $self->total_on_target_enrichment_factor . "\n";

    return $self;
}


sub _capture_total_on_target_bp {
    my $self = shift;

    # Total on-target capture space includes both unique and duplicate on-target
    # alignment bp.
    my $capture_total_on_target_bp =
        $self->capture_duplicate_bp_on_target +
        $self->capture_unique_bp_on_target;

    return $capture_total_on_target_bp;
}


sub _theoretical_max_enrichment_factor {
    my $self = shift;

    # Theoretical maximum enrichment factor (EF) describes enrichment factor
    # when the numerator (capture bp aligned) is at 100% ideal alignment. This
    # value gives a guideline for comparing calculated EF based on observed
    # data.
    #
    #                           (ideal)
    #      --------------------------------------------------------  =  Theoretical
    #                     (# bp in targets)                             Max. EF
    #                     ----------------- X 100
    #                      (# bp in genome)

    use constant IDEAL => 100;

    my $theoretical_max_enrichment_factor =
        IDEAL /
        (($self->target_total_bp / $self->genome_total_bp) * 100);

    return _round( $theoretical_max_enrichment_factor );
}



sub _unique_on_target_enrichment_factor {
    my $self = shift;

    # NOTE: We are including only unique target values here.
    #
    #                   (# bp on target unique)
    #      -------------------------------------------------- X 100
    #                      total # capture bp
    #      --------------------------------------------------------  =  unique on target
    #                     (# bp in targets)                             EF
    #                     ----------------- X 100
    #                      (# bp in genome)

    my $unique_on_target_enrichment_factor =
        (($self->capture_unique_bp_on_target / $self->capture_total_bp) * 100) /
        (($self->target_total_bp / $self->genome_total_bp) * 100);

    return _round( $unique_on_target_enrichment_factor );
}



sub _total_on_target_enrichment_factor {
    my $self = shift;

    # NOTE: We are including all on target values here, regardless of
    # redundancy.
    #
    #      (# bp on target unique + # bp on target duplicate)
    #      -------------------------------------------------- X 100
    #                      total # capture bp
    #      --------------------------------------------------------  = total on target
    #                     (# bp in targets)                            EF
    #                     ----------------- X 100
    #                      (# bp in genome)

    my $total_on_target_enrichment_factor =
        ((($self->capture_unique_bp_on_target + $self->capture_duplicate_bp_on_target) / $self->capture_total_bp) * 100) /
        (($self->target_total_bp / $self->genome_total_bp) * 100);

    return _round( $total_on_target_enrichment_factor );
}


sub _round { return sprintf( "%.1f", shift ) }

1;  # end of package.


=head1 NAME

CaptureEnrichmentFactor - Class for calculating exome capture Enrichment Factor metrics.


=head1 VERSION

version 0.1


=head1 SYNOPSIS

 #!/usr/bin/env perl

 # Simple test of CaptureEnrichmentFactor class.

 use strict;
 use warnings;
 use CaptureEnrichmentFactor;

 my $myEF = CaptureEnrichmentFactor->new(
                                         capture_unique_bp_on_target    =>  6_795_966_000,
                                         capture_duplicate_bp_on_target =>    834_157_200,
                                         capture_total_bp               => 13_613_105_800,
                                         target_total_bp                =>     45_203_256,
                                         genome_total_bp                =>  2_858_012_809,
                                        );

 my $theoretical_max_enrichment_factor  = $myEF->theoretical_max_enrichment_factor();
 my $unique_on_target_enrichment_factor = $myEF->unique_on_target_enrichment_factor();
 my $total_on_target_enrichment_factor  = $myEF->total_on_target_enrichment_factor();

 # NOTE: Should equal 63.2, 31.6, & 35.4 respectively.

 print join (
             "\t",
             $theoretical_max_enrichment_factor,
             $unique_on_target_enrichment_factor,
             $total_on_target_enrichment_factor,
            ) . "\n";

 __END__


=head1 DESCRIPTION

This class essentially wraps up calculations for capture Enrichment Factor calculations into a single call. There are 5 mandatory pieces of information which must be passed to the class constructor. All of these values are currently available from the current capture coverage pipeline incarnation. Mandatory values are: capture_unique_bp_on_target; capture_duplicate_bp_on_target; capture_total_bp; genome_total_bp; target_total_bp. Three metrics are assessed: 1) theoretical maximum enrichment factor; 2) unique on-target enrichment factor; and 3) total on-target enrichment factor.


=head1 INTERFACE

The following routines are supported:

=head2 B<theoretical_max_enrichment_factor>

Theoretical maximum enrichment factor (EF) describes enrichment factor when the numerator (capture bp aligned) is at 100% ideal alignment. This value gives a guideline for comparing calculated EF based on observed data.

 my $value = $myEF->theoretical_max_enrichment_factor();

 CALCULATION:

                      (ideal)
 --------------------------------------------------------  =  Theoretical
                (# bp in targets)                             Max. EF
                ----------------- X 100
                 (# bp in genome)



=head2 B<unique_on_target_enrichment_factor>

We are including only unique target values for this assessment.

 my $value = $myEF->unique_on_target_enrichment_factor();

 CALCULATION:

              (# bp on target unique)
 -------------------------------------------------- X 100
                 total # capture bp
 --------------------------------------------------------  =  unique on target
                (# bp in targets)                             EF
                ----------------- X 100
                 (# bp in genome)


=head2 B<total_on_target_enrichment_factor>

We are including all target values for this assessment, regardless of redundancy issues.

 my $value = $myEF->total_on_target_enrichment_factor();

CALCULATION:

 (# bp on target unique + # bp on target duplicate)
 -------------------------------------------------- X 100
                 total # capture bp
 --------------------------------------------------------  = total on target
                (# bp in targets)                            EF
                ----------------- X 100
                 (# bp in genome)

=head2 B<capture_unique_bp_on_target>

General value accessor.

  my $value = $myEF->capture_unique_bp_on_target();

=head2 B<capture_duplicate_bp_on_target>

General value accessor.

  my $value = $myEF->capture_duplicate_bp_on_target();

=head2 B<capture_total_bp>

General value accessor.

  my $value = $myEF->capture_total_bp();

=head2 B<genome_total_bp>

General value accessor.

  my $value = $myEF->genome_total_bp();

=head2 B<target_total_bp>

General value accessor.

  my $value = $myEF->target_total_bp();

=head1 CONFIGURATION AND ENVIRONMENT

CaptureEnrichmentFactor requires no independant configuration files or environment variables--aside from the Perl modules dependencies outlined below.


=head1 DEPENDENCIES

 strict
 warnings
 Carp
 constant

=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

None reported. Please report any bugs or feature requests to <twylie@genome.wustl.edu>.


=head1 AUTHOR

Todd Wylie

C<< <twylie@genome.wustl.edu> >>

L<< http://genome.wustl.edu  >>


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010, Todd Wylie C<< <twylie@genome.wustl.edu> >>. All rights reserved.

This software is licensed under the Artistic License 2.0, a copy of which
should have been provided with this distribution.

=head1 DISCLAIMER OF WARRANTY

THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS" AND
WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES. THE IMPLIED WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT ARE
DISCLAIMED TO THE EXTENT PERMITTED BY YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO
COPYRIGHT HOLDER OR CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE
PACKAGE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=head1 NOTE

This software was written using the latest version of GNU Emacs, the
extensible, real-time text editor. Please see
L<http://www.gnu.org/software/emacs> for more information and download
sources.
