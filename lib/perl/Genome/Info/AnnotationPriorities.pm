package Genome::Info::AnnotationPriorities;

use strict;
use warnings;


# NOTE - Genome::Model::Tools::Annotate::TranscriptVariants no longer uses this module.
# It has its own copy of these data structures so it can version control itself

my %variant_priorities =
(
    nonsense                        => 1,
    frame_shift                     => 2,
    frame_shift_del                 => 3,
    frame_shift_ins                 => 4,
    splice_site                     => 5,
    splice_site_del                 => 6,
    splice_site_ins                 => 7,
    in_frame_del                    => 8,
    in_frame_ins                    => 9,
    missense                        => 10,
    nonstop                         => 11,
    silent                          => 12,
    rna                             => 13,
    '5_prime_untranslated_region'   => 14,
    '3_prime_untranslated_region'   => 15,
    splice_region                   => 16,
    splice_region_del               => 17,
    splice_region_ins               => 18,
    intronic                        => 19,
    '5_prime_flanking_region'       => 20,
    '3_prime_flanking_region'       => 21,
    #consensus_error                 => 17,
);

my %transcript_status_priorities = 
(
    reviewed    => 1,
    validated   => 2,
    provisional => 3,
    predicted   => 4,
    model       => 5,
    inferred    => 6,
    known       => 7,
    novel       => 8,
    unknown     => 9,
);

my %transcript_source_priorities = 
(
    genbank => 1,
    ensembl => 2,
);

my %transcript_error_priorities = 
(
    no_errors                               => 1,
    gap_between_substructures               => 2,
    mismatch_between_exon_seq_and_reference => 3,
    bad_bp_length_for_coding_region         => 4,
    overly_large_intron                     => 5,
    rna_with_coding_region                  => 6,
    no_coding_region                        => 7,
    no_stop_codon                           => 8,
    pseudogene                              => 9,
    no_start_codon                          => 10,
);

sub variant_priorities{
    return %variant_priorities;
}

sub transcript_status_priorities{
    return %transcript_status_priorities;
}

sub transcript_source_priorities{
    return %transcript_source_priorities;
}

sub transcript_error_priorities {
    return %transcript_error_priorities;
}

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/annotation_update/Genome/Info/AnnotationPriorities.pm $
#$Id: AnnotationPriorities.pm 55387 2010-02-15 22:16:40Z adukes $

1;
