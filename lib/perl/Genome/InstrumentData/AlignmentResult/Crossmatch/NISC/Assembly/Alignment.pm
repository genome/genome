
### Adapted from code by Nancy Hansen <nhansen@mail.nih.gov>


############################################################
# Alignment.pm: Object to describe and manipulate information
#  	regarding a alignment of two sequences
#
# Author:       Nancy F. Hansen
# Version: $Id: Alignment.pm,v 1.7 2005/04/06 14:25:14 nhansen Exp $
############################################################

package Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::Alignment;
use strict;
use Carp;

use Bio::Seq;

#use Opendb;
#use NISC::Comp::Disc;
#use NISC::Comp::Discrepancy;
#use NISC::Comp::FSDiscrepancy;
#use NISC::Comp::WGSDiscrepancy;

###############################################################################
# Alignment constructor
#
# INPUT: parameters:
#    -align_type - name of program that generated the alignment ('crossmatch', 
#         'blast')
#    -score - a number indicating the strength of the alignment.
#    -substitution_rate - percentage of bases substituted
#    -deletion_rate - percentage of bases deleted
#    -insertion_rate - percentage of bases inserted
#    -comp - U or C, depending on whether the query is aligned with or against
#         the hit
#    -query_file - name of the query file
#    -query - name of the query sequence
#    -query_start - coordinate where the alignment begins on the query
#    -query_end - coordinate where the alignment ends on the query
#    -query_remaining - unaligned sequence at the end of the query
#    -hit_file - name of the hit file
#    -hit - name of the hit sequence
#    -hit_start - coordinate where the alignment begins on the hit
#    -hit_end - coordinate where the alignment ends on the hit
#    -hit_remaining - unaligned sequence at the end of the hit
#    -align_string - string containing align from ace file.
#    -expt_id - tells which dataset we're working with
#
# OUTPUT: a new Alignment object
###############################################################################

sub new {

    my $this = shift;
    my %params = @_;
 
    my $align_type = $params{-align_type};
    my $score = $params{-score};
    my $substitution_rate = $params{-substitution_rate};
    my $deletion_rate = $params{-deletion_rate};
    my $insertion_rate = $params{-insertion_rate};
    my $query_file = $params{-query_file};
    my $query = $params{-query};
    my $query_start = $params{-query_start};
    my $query_end = $params{-query_end};
    my $query_remaining = $params{-query_remaining};
    my $hit_file = $params{-hit_file};
    my $hit = $params{-hit};
    my $hit_start = $params{-hit_start};
    my $hit_end = $params{-hit_end};
    my $hit_remaining = $params{-hit_remaining};
    my $comp = $params{-comp} || '';
    my $query_string = $params{-query_string};
    my $hit_string = $params{-hit_string};
    my $expt_id = $params{-expt_id};

    my $self = { align_type => $align_type, score => $score, 
                 substitution_rate => $substitution_rate,
                 deletion_rate => $deletion_rate, 
                 insertion_rate => $insertion_rate, 
                 query_file => $query_file, 
                 query => $query, query_start => $query_start, 
                 query_end => $query_end, 
                 query_remaining => $query_remaining, 
                 hit_file => $hit_file, 
                 hit => $hit, hit_start => $hit_start,
                 hit_end => $hit_end, hit_remaining => $hit_remaining, 
                 comp => $comp, 
                 query_string => $query_string, hit_string => $hit_string,
                 expt_id => $expt_id };

    my $class = ref($this) || $this;
    bless $self, $class;
    return($self);

} # end new

###############################################################################
# Subroutine to return the type of alignment on which this Alignment
#     object is based.
#
# INPUT: Alignment object 
# OUTPUT: program name (scalar) if defined, else empty string!
###############################################################################
sub align_type {
    my $self = shift;

    return $self->{align_type} || '';
    
} ## end align_type

###############################################################################
# Subroutine to return whether an alignment is complemented or not
#
# INPUT: Alignment object 
# OUTPUT: comp value (scalar 'U' or 'C', or empty string if it wasn't specified)
###############################################################################
sub comp {
    my $self = shift;
    my $comp = $self->{comp} || '';

    return $comp;
    
} ## end comp

###############################################################################
# Subroutine to return the score of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: score value (scalar)
###############################################################################
sub score {
    my $self = shift;

    if (defined (my $new_score = shift))
    {
        $self->{score} = $new_score;
    }

    return $self->{score};
    
} ## end score

###############################################################################
# Subroutine to return a reference to an array of subalignments for this
#     alignment (only used for blastz currently)
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: reference to an array of Alignment objects
###############################################################################
sub sub_alignments {
    my $self = shift;

    if (defined (my $new_sub_alignments = shift))
    {
        $self->{sub_alignments} = $new_sub_alignments;
    }

    return $self->{sub_alignments};
    
} ## end sub_alignments

###############################################################################
# Subroutine to return the score of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: score value (scalar)
###############################################################################
sub score_global {
    my $self = shift;

    if (defined (my $new_score = shift))
    {
        $self->{score} = $new_score;
    }

    return $self->{score};
    
} ## end score_global

###############################################################################
# Subroutine to return the substitution_rate of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: substitution_rate value (scalar)
###############################################################################
sub substitution_rate {
    my $self = shift;

    if (defined (my $new_substitution_rate = shift))
    {
        $self->{substitution_rate} = $new_substitution_rate;
    }

    return $self->{substitution_rate};
    
} ## end substitution_rate

###############################################################################
# Subroutine to return the deletion_rate of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: deletion_rate value (scalar)
###############################################################################
sub deletion_rate {
    my $self = shift;

    if (defined (my $new_deletion_rate = shift))
    {
        $self->{deletion_rate} = $new_deletion_rate;
    }

    return $self->{deletion_rate};
    
} ## end deletion_rate

###############################################################################
# Subroutine to return the insertion_rate of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: insertion_rate value (scalar)
###############################################################################
sub insertion_rate {
    my $self = shift;

    if (defined (my $new_insertion_rate = shift))
    {
        $self->{insertion_rate} = $new_insertion_rate;
    }

    return $self->{insertion_rate};
    
} ## end insertion_rate

###############################################################################
# Subroutine to return the query_file of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_file value (scalar)
###############################################################################
sub query_file {
    my $self = shift;

    if (defined (my $new_query_file = shift))
    {
        $self->{query_file} = $new_query_file;
    }

    return $self->{query_file};
    
} ## end query_file

###############################################################################
# Subroutine to return the query of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query value (scalar)
###############################################################################
sub query {
    my $self = shift;

    if (defined (my $new_query = shift))
    {
        $self->{query} = $new_query;
    }

    return $self->{query};
    
} ## end query

###############################################################################
# Subroutine to return the query_start of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_start value (scalar)
###############################################################################
sub query_start {
    my $self = shift;

    if (defined (my $new_query_start = shift))
    {
        $self->{query_start} = $new_query_start;
    }

    return $self->{query_start};
    
} ## end query_start

###############################################################################
# Subroutine to return the query_end of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_end value (scalar)
###############################################################################
sub query_end {
    my $self = shift;

    if (defined (my $new_query_end = shift))
    {
        $self->{query_end} = $new_query_end;
    }

    return $self->{query_end};
    
} ## end query_end

###############################################################################
# Subroutine to return the query_remaining of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_remaining value (scalar)
###############################################################################
sub query_remaining {
    my $self = shift;

    if (defined (my $new_query_remaining = shift))
    {
        $self->{query_remaining} = $new_query_remaining;
    }

    return $self->{query_remaining};
    
} ## end query_remaining

###############################################################################
# Subroutine to return the query's length in an alignment
#
# INPUT: Alignment object
# OUTPUT: query_align_length value (scalar)
###############################################################################
sub query_align_length {
    my $self = shift;

    return $self->query_end() - $self->query_start();
    
} ## end query_align_length

###############################################################################
# Subroutine to return the query's length
#
# INPUT: Alignment object
# OUTPUT: query_length value (scalar)
###############################################################################
sub query_length {
    my $self = shift;

    return $self->query_end() + $self->query_remaining();
    
} ## end query_length

###############################################################################
# Subroutine to return the query_string of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_string value (scalar)
###############################################################################
sub query_string {
    my $self = shift;

    if (defined (my $new_query_string = shift))
    {
        $self->{query_string} = $new_query_string;
    }

    return $self->{query_string};
    
} ## end query_string

###############################################################################
# Subroutine to return the sequence of the query to the left of the 
#     alignment (when the hit is uncomplemented).
#
# INPUT: Alignment object, -query_seq => Bio::Seq object for query.
# OUTPUT: string representing query sequence to the left.
###############################################################################

sub query_left_side
{
    my $self = shift;
    my %params = @_;
    my $seq_obj = $params{-query_seq}
        or croak "Must pass a Bio::Seq object to query_left_side as -query_seq!\n";

    if ($self->comp() eq 'U') # uncomplemented
    {
        my $start = $self->query_start();
        if ($start - 1) # there's sequence before the alignment:
        {
            return $seq_obj->subseq(1, $start-1);
        }
        else
        {
            return '';
        }
    }
    else # complemented
    {
        my $end = $self->query_end();
        my $length = $seq_obj->length();

        if ($length - $end)
        {
            my $left_side_comp = $seq_obj->subseq($end+1, $length);
            my $comp_seq = Bio::Seq->new( -seq => $left_side_comp );
            $comp_seq->alphabet('dna');
            return $comp_seq->revcom()->seq();
        }
        else
        {
            return '';
        }
    }

} ## end query_left_side

###############################################################################
# Subroutine to return the sequence of the query to the right of the 
#     alignment (when the hit is uncomplemented).
#
# INPUT: Alignment object, -query_seq => Bio::Seq object for query.
# OUTPUT: string representing query sequence to the right.
###############################################################################
sub query_right_side
{
    my $self = shift;
    my %params = @_;
    my $seq_obj = $params{-query_seq}
        or croak "Must pass a Bio::Seq object to query_left_side as -query_seq!\n";

    if ($self->comp() eq 'U') # uncomplemented
    {
        my $end = $self->query_end();
        my $length = $seq_obj->length();

        if ($length - $end) # there's sequence after the alignment:
        {
            return $seq_obj->subseq($end+1, $length);
        }
        else
        {
            return '';
        }
    }
    else # complemented
    {
        my $start = $self->query_start();

        if ($start - 1)
        {
            my $right_side_comp = $seq_obj->subseq(1, $start-1);
            my $comp_seq = Bio::Seq->new( -seq => $right_side_comp );
            $comp_seq->alphabet('dna');
            return $comp_seq->revcom()->seq();
        }
        else
        {
            return '';
        }
    }
} ## end query_right_side

###############################################################################
# Subroutine to return the hit_file of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_file value (scalar)
###############################################################################
sub hit_file {
    my $self = shift;

    if (defined (my $new_hit_file = shift))
    {
        $self->{hit_file} = $new_hit_file;
    }

    return $self->{hit_file};
    
} ## end hit_file

###############################################################################
# Subroutine to return the hit of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit value (scalar)
###############################################################################
sub hit {
    my $self = shift;

    if (defined (my $new_hit = shift))
    {
        $self->{hit} = $new_hit;
    }

    return $self->{hit};
    
} ## end hit

###############################################################################
# Subroutine to return the hit_start of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_start value (scalar)
###############################################################################
sub hit_start {
    my $self = shift;

    if (defined (my $new_hit_start = shift))
    {
        $self->{hit_start} = $new_hit_start;
    }

    return $self->{hit_start};
    
} ## end hit_start

###############################################################################
# Subroutine to return the hit_end of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_end value (scalar)
###############################################################################
sub hit_end {
    my $self = shift;

    if (defined (my $new_hit_end = shift))
    {
        $self->{hit_end} = $new_hit_end;
    }

    return $self->{hit_end};
    
} ## end hit_end

###############################################################################
# Subroutine to return the coordinate of the left end of the hit of an 
#    alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_left_end value (scalar)
###############################################################################
sub hit_left_end {
    my $self = shift;

    if (!defined ($self->{hit_left_end}))
    {
        my $comp = $self->comp();
        if ($comp eq 'U')
        {
            $self->{hit_left_end} = $self->hit_start();
        }
        else
        {
            $self->{hit_left_end} = $self->hit_end();
        }
    }
    return $self->{hit_left_end};
    
} ## end hit_left_end

###############################################################################
# Subroutine to return the coordinate of the right end of the hit of an 
#    alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_right_end value (scalar)
###############################################################################
sub hit_right_end {
    my $self = shift;

    if (!defined ($self->{hit_right_end}))
    {
        my $comp = $self->comp();
        if ($comp eq 'U')
        {
            $self->{hit_right_end} = $self->hit_end();
        }
        else
        {
            $self->{hit_right_end} = $self->hit_start();
        }
    }
    return $self->{hit_right_end};
    
} ## end hit_right_end

###############################################################################
# Subroutine to return the coordinate of the left end of the query of an 
#    alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_left_end value (scalar)
###############################################################################
sub query_left_end {
    my $self = shift;

    if (!defined ($self->{query_left_end}))
    {
        my $comp = $self->comp();
        if ($comp eq 'U')
        {
            $self->{query_left_end} = $self->query_start();
        }
        else
        {
            $self->{query_left_end} = $self->query_end();
        }
    }
    return $self->{query_left_end};
    
} ## end query_left_end

###############################################################################
# Subroutine to return the coordinate of the right end of the query of an 
#    alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: query_right_end value (scalar)
###############################################################################
sub query_right_end {
    my $self = shift;

    if (!defined ($self->{query_right_end}))
    {
        my $comp = $self->comp();
        if ($comp eq 'U')
        {
            $self->{query_right_end} = $self->query_end();
        }
        else
        {
            $self->{query_right_end} = $self->query_start();
        }
    }
} ## end query_right_end

###############################################################################
# Subroutine to return the hit_remaining of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_remaining value (scalar)
###############################################################################
sub hit_remaining {
    my $self = shift;

    if (defined (my $new_hit_remaining = shift))
    {
        $self->{hit_remaining} = $new_hit_remaining;
    }

    return $self->{hit_remaining};
    
} ## end hit_remaining

###############################################################################
# Subroutine to return the hit's length
#
# INPUT: Alignment object
# OUTPUT: hit_length value (scalar)
###############################################################################
sub hit_length {
    my $self = shift;

    return $self->hit_right_end() + $self->hit_remaining();
    
} ## end hit_length

###############################################################################
# Subroutine to return the hit_string of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: hit_string value (scalar)
###############################################################################
sub hit_string {
    my $self = shift;

    if (defined (my $new_hit_string = shift))
    {
        $self->{hit_string} = $new_hit_string;
    }

    return $self->{hit_string};
    
} ## end hit_string

###############################################################################
# Subroutine to return hashes representing polymorphisms found in this 
#    alignment
#
# INPUT: Alignment object
# OUTPUT: reference to a list of hash references, with fields:
#     hit - name of hit (from Alignment object)
#     query - name of query (from Alignment object)
#     hit_left_flank_end - coordinate to left of polymorphism in hit
#     hit_right_flank_start - coordinate to right of polymorphism in hit
#     query_left_flank_end - coordinate to left of polymorphism in query
#     query_right_flank_start - coordinate to right of polymorphism in query
#     hit_allele - sequence between flanks in hit
#     query_allele - sequence between flanks in hit
#     strand - 'U' or 'C' (from comp of Alignment object)
###############################################################################
sub polymorphisms {
    my $self = shift;

    my $hit = $self->hit();
    my $query = $self->query();
    my $hit_start = $self->hit_left_end();
    my $hit_end = $self->hit_right_end();
    my $query_start = $self->query_left_end();
    my $query_end = $self->query_right_end();
    my $hit_pos = $hit_end;
    my $query_pos = $query_end;
    my $hit_string = $self->hit_string();
    my $query_string = $self->query_string();
    my $strand = $self->comp();

    my ($query_allele, $hit_allele) = ('', '');
    my ($query_lfe, $query_rfs, $hit_lfe, $hit_rfs);
    my $in_poly;

    my @polys = ();
    while ($hit_string)
    {
        my $query_base = chop $query_string;
        my $hit_base = chop $hit_string;

        if (($query_base ne $hit_base) && ($query_base !~ /N/i) && ($hit_base !~ /N/i))
        {
            if ($in_poly) # already within polymorphism--continue
            {
                $query_allele = $query_base.$query_allele;
                $hit_allele = $hit_base.$hit_allele;
                #print "At $query/$hit ($strand) $query_pos/$hit_pos: $query_allele/$hit_allele\n";
            }
            else # begin new poly
            {
                $query_allele = $query_base;
                $hit_allele = $hit_base;
                $hit_rfs = ($hit_start > $hit_end) ? 
                               $hit_pos - 1 : $hit_pos + 1;
                $query_rfs = ($query_start > $query_end) ? 
                               $query_pos - 1 : $query_pos + 1;
                $in_poly = 1;
            }
        }
        else
        {
            if ($in_poly) # need to "close" polymorphism
            {
                $hit_allele =~ s:\*::g;
                $query_allele =~ s:\*::g;
                push @polys, {'hit' => $hit,
                              'query' => $query,
                              'hit_left_flank_end' => $hit_pos,
                              'hit_right_flank_start' => $hit_rfs,
                              'query_left_flank_end' => $query_pos,
                              'query_right_flank_start' => $query_rfs,
                              'hit_allele' => $hit_allele,
                              'query_allele' => $query_allele,
                              'strand' => $strand};
                ($hit_allele, $query_allele) = ('','');
                $in_poly = 0;
                $hit_lfe = undef;
                $hit_rfs = undef;
                $query_lfe = undef;
                $query_rfs = undef;
            }
        }

        if (($query_end > $query_start) && ($query_base ne '*'))
        {
            $query_pos--;
        }
        elsif ($query_base ne '*')
        {
            $query_pos++;
        }
        if (($hit_end > $hit_start) && ($hit_base ne '*'))
        {
            $hit_pos--;
        }
        elsif ($hit_base ne '*')
        {
            $hit_pos++;
        }
    }

    return \@polys;

} ## end polymorphisms

###############################################################################
# Subroutine to return the average offset of an alignment (equal to 
#    the average of hit_start - query_start and hit_end - query_end)
#
# INPUT: Alignment object
# OUTPUT: offset value (scalar number)
###############################################################################
sub offset {
    my $self = shift;

    my $q_start = $self->query_start();
    my $q_end = $self->query_end();
    my $h_start = $self->hit_start();
    my $h_end = $self->hit_end();

    my $offset = (($h_start - $q_start) + ($h_end - $q_end))/2;

    return $offset;
    
} ## end offset

###############################################################################
# Subroutine to return the expt_id of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: expt_id value (scalar)
###############################################################################
sub expt_id {
    my $self = shift;

    if (defined (my $new_expt_id = shift))
    {
        $self->{expt_id} = $new_expt_id;
    }

    return $self->{expt_id};
    
} ## end expt_id

###############################################################################
# Subroutine to return the discrepancies of an alignment
#
# INPUT: Alignment object, optional argument sets value
# OUTPUT: discrepancies value (scalar)
###############################################################################
sub discrepancies {
    my $self = shift;

    if (defined (my $new_discrepancies = shift))
    {
        $self->{discrepancies} = $new_discrepancies;
    }

    return $self->{discrepancies};
    
} ## end discrepancies

###############################################################################
# Subroutine to translate the coordinates on the hit by the amount given
#     in the argument.
#
# INPUT: Alignment object, argument tells how many bases to translate
#        n.b.-this routine does not alter the value of hit_remaining!
#
# OUTPUT: Alignment object with new hit coordinates
###############################################################################
sub translate_hit_coords {
    my $self = shift;
    my $bases = shift;

    if (defined ($self->{hit_start}))
    {
        $self->{hit_start} += $bases;
    }

    if (defined ($self->{hit_end}))
    {
        $self->{hit_end} += $bases;
    }

    if (defined ($self->{hit_left_end}))
    {
        $self->{hit_left_end} += $bases;
    }
    
    return $self;

} ## end translate_hit_coords

###############################################################################
# Subroutine to return 1 if two alignments are consistent with plasmid
#     ends as the query--if they hit the same entry, they must be
#     ordered and oriented consistently with the correct distance 
#     (which can be set with -min_distance and -max_distance)
#
# INPUT: Two Alignment objects, optional -min_distance and -max_distance
# OUTPUT: 1 if they could be plasmid ends, 0 otherwise
###############################################################################

sub plasmid_ends {

    my $align1 = shift;
    my $align2 = shift;
    my %params = @_;

    my $min_distance = defined ($params{-min_distance}) ? 
                            $params{-min_distance} : 1500;
    my $max_distance = defined ($params{-max_distance}) ? 
                            $params{-max_distance} : 7000;

    my $hit1 = $align1->hit();
    my $hit2 = $align2->hit();

    if ($hit1 eq $hit2) # on the same contig
    {
        my $comp1 = $align1->comp();
        my $comp2 = $align2->comp();

        if ($comp1 eq $comp2)
        {
             return 0;
        }
        else # are they the right distance?
        {
            my $distance;
            if ($comp2 eq 'C') # second read complemented
            {
                 $distance = $align2->hit_start() + $align2->query_start() - 2
                           - $align1->hit_start() + $align1->query_start();
            } 
            else
            {
                 $distance = $align1->hit_start() + $align1->query_start() - 2
                           - $align2->hit_start() + $align2->query_start();
            }
            my $query = $align1->query();
            if (($distance > $min_distance) && ($distance <= $max_distance))
            {
                 return 1;
            }
            else
            {
                 return 0;
            }
        }
    }
    else # two different contigs (are they at ends directed outward?)
    {
        my $distance = 0;
        foreach my $align ($align1, $align2)
        {
            if ($align->comp() eq 'U')
            {
                $distance += $align->hit_remaining() + $align->query_end() - 1;
            }
            else
            {
                $distance += $align->hit_start() + $align->query_start() - 1;
            }
        }
        my $query1 = $align1->query();
        if ($distance <= $max_distance)
        {
             return 1;
        }
        else
        {
             return 0;
        }
    }

} ## end plasmid_ends

###############################################################################
# Subroutine to return a SAM-formatted CIGAR string for this alignment
#
# INPUT: NISC::Assembly::Alignment object, 
#         -hard_clip - with a true value for hard clipped alignments
#         -max_bases - only report cigar string for first max_bases bases.
# OUTPUT: cigar string (scalar string)
###############################################################################

sub cigar_string {

    my $self = shift;
    my %params = @_;
    my $hard_clip = $params{'-hard_clip'};
    my $max_bases = $params{'-max_bases'};

    my @operations = (); # will store CIGAR operations

    my $query_string = $self->query_string();
    my $hit_string = $self->hit_string();

    my ($in_match, $in_insertion, $in_deletion) = (0, 0, 0);
    my ($match_bases, $ins_bases, $del_bases) = (0, 0, 0);

    while ($query_string)
    {
        my $query_base = chop $query_string;
        my $hit_base = chop $hit_string;

        if ($query_base eq '*') # deletion from reference
        {
            $del_bases++;
            $in_deletion = 1;
            if ($in_match)
            {
                unshift @operations, $match_bases."M";
                $in_match = 0;
                $match_bases = 0;
            }
            elsif ($in_insertion)
            {
                unshift @operations, $ins_bases."I";
                $in_insertion = 0;
                $ins_bases = 0;
            }
        }
        elsif ($hit_base eq '*') # insertion to reference
        {
            $ins_bases++;
            $in_insertion = 1;
            if ($in_match)
            {
                unshift @operations, $match_bases."M";
                $in_match = 0;
                $match_bases = 0;
            }
            elsif ($in_deletion)
            {
                unshift @operations, $del_bases."D";
                $in_deletion = 0;
                $del_bases = 0;
            }

        }
        else # match or mismatch
        {
            $match_bases++;
            $in_match = 1;
            if ($in_insertion)
            {
                unshift @operations, $ins_bases."I";
                $in_insertion = 0;
                $ins_bases = 0;
            }
            elsif ($in_deletion)
            {
                unshift @operations, $del_bases."D";
                $in_deletion = 0;
                $del_bases = 0;
            }
        }
    }

    # a few last operations!

    unshift @operations, ($in_match) ? $match_bases."M" :
                      ($in_insertion) ? $ins_bases."I" :
                      ($in_deletion) ? $del_bases."D" : '';

    my $beginning_bases = $self->query_start() - 1;
    if ($beginning_bases)
    {
        if ($self->comp() eq 'U')
        {
            unshift @operations, ($hard_clip) ? $beginning_bases."H" : $beginning_bases."S";
        }
        else
        {
            push @operations, ($hard_clip) ? $beginning_bases."H" : $beginning_bases."S";
        }
    }
    my $end_bases = $self->query_remaining();
    if ($end_bases)
    {
        if ($self->comp() eq 'U')
        {
            push @operations, ($hard_clip) ? $end_bases."H" : $end_bases."S";
        }
        else
        {
            unshift @operations, ($hard_clip) ? $end_bases."H" : $end_bases."S";
        }
    }

    if ($max_bases) # only include first max_bases operations
    {
        my @new_ops = ();
        my $new_bases = 0;
        my $comp = $self->comp();
        while ($new_bases < $max_bases && @operations)
        {
            my $next_op = ($comp eq 'U') ? shift @operations : pop @operations;
            my ($next_op_bases, $next_op_char) = ($next_op =~ /^(\d+)([A-Z])$/) ? ($1, $2) : ('NA', 'NA');
            if ($next_op_char eq 'D') # don't count deletion bases
            {
                push @new_ops, $next_op if ($comp eq 'U');
                unshift @new_ops, $next_op if ($comp eq 'C');
            }
            elsif ($next_op_bases + $new_bases > $max_bases)
            {
                my $new_op_bases = $max_bases - $new_bases;
                push @new_ops, "$new_op_bases$next_op_char" if ($comp eq 'U');
                unshift @new_ops, "$new_op_bases$next_op_char" if ($comp eq 'C');
                last;
            }
            else
            {
                push @new_ops, $next_op if ($comp eq 'U');
                unshift @new_ops, $next_op if ($comp eq 'C');
                $new_bases += $next_op_bases;
            }
        }
        @operations = @new_ops;
    }

    my $cigar_string = join '', @operations;

    return $cigar_string;

} ## end cigar_string

###############################################################################
# Subroutine to enter this Alignment into a database table.
#
# INPUT: NISC::Assembly::Alignment object, 
#         -clone_id specifies clone_id value (or
#        in the case of WGS alignments, the species),
#        -dbh passes database handle (not optional yet!)
#        -table inserts into the specified table (default is zoo_alignment)
# OUTPUT: new id value
###############################################################################

sub db_insert {

    my $self = shift;
    my %params = @_;
    my $clone_id = $params{'-clone_id'};
    my $species = $params{'-clone_id'};
    my $table = $params{'-table'} || 'zoo_alignment';
    my $dbh = $params{'-dbh'};

    my $query_file = $self->query_file();
    my $query = $self->query();
    my $query_start = $self->query_start();
    my $query_end = $self->query_end();
    my $query_remaining = $self->query_remaining();
    my $hit_file = $self->hit_file();
    my $hit = $self->hit();
    my $hit_start = $self->hit_start();
    my $hit_end = $self->hit_end();
    my $hit_remaining = $self->hit_remaining();
    my $comp = $self->comp();
    my $expt_id = $self->expt_id();

    my $id_select = qq! SELECT MAX(id)
                        FROM nhansen.$table !;

    my $id_sth = $dbh->prepare($id_select)
         or croak "Couldn\'t prepare statement $id_select!\n";
    my $id_rv = $id_sth->execute()
         or croak "Couldn\'t execute statement $id_select!\n";

    my $ra_id = $id_sth->fetchrow_arrayref();

    my $id = ($ra_id->[0]) ? $ra_id->[0] + 1 : 1;

    my $insert;

    if ($table eq 'alignment')
    {
        $insert = qq! INSERT INTO nhansen.$table
                     (id, query_file, query, query_start,
                      query_end, query_remaining, hit_file,
                      hit, hit_start, hit_end, hit_remaining,
                      comp, expt_id) VALUES
                     ($id, '$query_file', '$query', $query_start,
                      $query_end, $query_remaining, '$hit_file',
                      '$hit', $hit_start, $hit_end, $hit_remaining,
                      '$comp', $expt_id) ! ;
    }
    elsif ($table ne 'wgs_alignment')
    {
        $insert = qq! INSERT INTO nhansen.$table
                     (id, clone_id, query, query_start,
                      query_end, query_remaining, hit,
                      hit_start, hit_end, hit_remaining,
                      comp, expt_id) VALUES
                     ($id, $clone_id, '$query', $query_start,
                      $query_end, $query_remaining, '$hit',
                      $hit_start, $hit_end, $hit_remaining,
                      '$comp', $expt_id) ! ;
    }
    else
    {
        $insert = qq! INSERT INTO nhansen.$table
                     (id, species, query, query_start,
                      query_end, query_remaining, hit,
                      hit_start, hit_end, hit_remaining,
                      comp, expt_id) VALUES
                     ($id, '$species', '$query', $query_start,
                      $query_end, $query_remaining, '$hit',
                      $hit_start, $hit_end, $hit_remaining,
                      '$comp', $expt_id) ! ;
    }

    my $sth = $dbh->prepare($insert)
         or croak "Couldn\'t prepare statement $insert!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $insert!\n";

    return $id;

} ## end db_insert

###############################################################################
# Subroutine to retrieve Alignment objects from the niscprod zoo_alignment 
#     table.
#
# INPUT: Class, -clone_id specifies clone_id value,
#        -expt_id specifies which dataset 
#        -table specifies which table to select from (default zoo_alignment)
#        -dbh passes database handle (not optional yet!)
# OUTPUT: reference to an array of NISC::Assembly::Alignment objects
###############################################################################

sub db_select {

    my $class = shift;
    my %params = @_;
    my $species = $params{'-species'};
    my $clone_id = $params{'-clone_id'};
    my $table = $params{'-table'} || 'zoo_alignment';
    my $expt_id = $params{'-expt_id'};
    my $query_file = $params{'-query_file'};
    my $dbh = $params{'-dbh'};

    if ($table eq 'fs_alignment')
    {
        my $ra_alignments = $class->db_fs_select(-clone_id => $clone_id, 
                       -expt_id => $expt_id, -dbh => $dbh);
        return $ra_alignments;
    }
    elsif ($table eq 'sanger_alignment')
    {
        my $ra_alignments = $class->db_sanger_select(-clone_id => $clone_id, 
                       -expt_id => $expt_id, -dbh => $dbh);
        return $ra_alignments;
    }
    elsif ($table eq 'wgs_alignment')
    {
        my $ra_alignments = $class->db_wgs_select(-species => $species, 
                       -expt_id => $expt_id, -dbh => $dbh);
        return $ra_alignments;
    }
    elsif ($table eq 'zoo_alignment')
    {
        my $ra_alignments = $class->db_zoo_select(-clone_id => $clone_id, 
                       -expt_id => $expt_id, -dbh => $dbh);
        return $ra_alignments;
    }


    my $select = qq! SELECT a.id, a.query_file, a.query, 
                            a.query_start, a.query_end,
                            a.query_remaining, a.hit_file,
                            a.hit, a.hit_start, a.hit_end,
                            a.hit_remaining, a.comp, 
                            d.id, d.query_position,
                            d.query_base, d.query_qual, 
                            d.hit_position, d.hit_base, 
                            d.hit_qual, d.type
                     FROM nhansen.alignment a, nhansen.disc d
                     WHERE  a.id = d.alignment_id (+) !;

    if (defined ($expt_id))
    {
        $select .= qq! AND a.expt_id = $expt_id !;
    }

    if (defined ($query_file))
    {
        $select .= qq! AND a.query_file = '$query_file' !;
    }

    my $sth = $dbh->prepare($select)
         or croak "Couldn\'t prepare statement $select!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $select!\n";

    my %discreps;
    my %aligns;
    while (my $ra_align = $sth->fetchrow_arrayref())
    {
        my ($id, $query_file, $query, $query_start, $query_end,
             $query_remaining, $hit_file, $hit, $hit_start, 
             $hit_end, $hit_remaining, $comp, $disc_id, 
             $query_position, $query_base, $query_qual,
             $hit_position, $hit_base, $hit_qual, $type) = @{$ra_align};

        if (!defined ($aligns{$id}))
        {
            $aligns{$id} = $class->new(-query_file => $query_file,
                                       -query => $query, 
                                       -query_start => $query_start,
                                       -query_end => $query_end, 
                                       -query_remaining => $query_remaining,
                                       -hit_file => $hit_file, 
                                       -hit => $hit, 
                                       -hit_start => $hit_start, 
                                       -hit_end => $hit_end,
                                       -hit_remaining => $hit_remaining, 
                                       -comp => $comp,
                                       -expt_id => $expt_id);
        }

        if (defined ($disc_id))
        {
             my $discrep = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Comp::Disc->new(
                             -id => $disc_id,
                             -alignment_id => $id,
                             -query_position => $query_position,
                             -query_base => $query_base,
                             -query_qual => $query_qual,
                             -hit_position => $hit_position,
                             -hit_base => $hit_base,
                             -hit_qual => $hit_qual,
                             -type => $type );

             push @{$discreps{$id}}, $discrep;
        }
    }

    my @alignments = ();
    foreach my $id (keys %aligns)
    {
         my $ra_discreps = $discreps{$id} || [];
         $aligns{$id}->discrepancies($ra_discreps);
         push @alignments, $aligns{$id};
    }

    return \@alignments;

} ## end db_select

###############################################################################
# Subroutine to retrieve Alignment objects from the niscprod zoo_alignment 
#     table.
#
# INPUT: Class, -clone_id specifies clone_id value,
#        -expt_id specifies which dataset 
#        -dbh passes database handle (not optional yet!)
# OUTPUT: reference to an array of NISC::Assembly::Alignment objects
###############################################################################

sub db_zoo_select {

    my $class = shift;
    my %params = @_;
    my $clone_id = $params{'-clone_id'};
    my $expt_id = $params{'-expt_id'};
    my $dbh = $params{'-dbh'};

    if (!defined ($dbh))
    {
        my ($lda) = Opendb::OpenDB();
        $dbh = ${$lda};
    }

    my $select = qq! SELECT a.id, a.clone_id, a.query, a.query_start, a.query_end,
                            a.query_remaining, a.hit, a.hit_start, a.hit_end,
                            a.hit_remaining, a.comp, d.zg_contig, d.zg_position,
                            d.zg_base, d.zg_qual, d.fn_position, d.fn_base,
                            d.fn_qual, a.expt_id
                     FROM nhansen.zoo_alignment a, nhansen.discrepancy d
                     WHERE  a.id = d.zoo_alignment_id (+)
                     AND    a.clone_id = $clone_id 
                     AND    a.expt_id = $expt_id !;

    my $sth = $dbh->prepare($select)
         or croak "Couldn\'t prepare statement $select!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $select!\n";

    my %discreps;
    my %aligns;
    while (my $ra_align = $sth->fetchrow_arrayref())
    {
        my ($id, $clone_id, $query, $query_start, $query_end,
             $query_remaining, $hit, $hit_start, $hit_end, $hit_remaining,
             $comp, $zg_contig, $zg_position, $zg_base, $zg_qual, $fn_position,
             $fn_base, $fn_qual, $expt_id) = @{$ra_align};

        if (!defined ($aligns{$id}))
        {
            $aligns{$id} = $class->new(-query => $query, -query_start => $query_start,
                             -query_end => $query_end, -query_remaining => $query_remaining,
                             -hit => $hit, -hit_start => $hit_start, -hit_end => $hit_end,
                             -hit_remaining => $hit_remaining, -comp => $comp,
                             -expt_id => $expt_id);
        }

        if ($zg_position)
        {
             my $discrep = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Comp::Discrepancy->new(-clone_id => $clone_id, -zg_contig => $zg_contig,
                             -zg_position => $zg_position, -zg_base => $zg_base,
                             -zg_qual => $zg_qual, -fn_position => $fn_position,
                             -fn_base => $fn_base, -fn_qual => $fn_qual);

             push @{$discreps{$id}}, $discrep;
        }
    }

    my @alignments = ();
    foreach my $id (keys %aligns)
    {
         my $ra_discreps = $discreps{$id} || [];
         $aligns{$id}->discrepancies($ra_discreps);
         push @alignments, $aligns{$id};
    }

    return \@alignments;

} ## end db_zoo_select

###############################################################################
# Subroutine to retrieve Alignment objects from the niscprod fs_alignment 
#     table.
#
# INPUT: Class, -clone_id specifies clone_id value,
#        -expt_id specifies which dataset 
#        -dbh passes database handle (not optional yet!)
# OUTPUT: reference to an array of NISC::Assembly::Alignment objects
###############################################################################

sub db_fs_select {

    my $class = shift;
    my %params = @_;
    my $clone_id = $params{'-clone_id'};
    my $expt_id = $params{'-expt_id'};
    my $dbh = $params{'-dbh'};

    if (!defined ($dbh))
    {
        my ($lda) = Opendb::OpenDB();
        $dbh = ${$lda};
    }

    my $select = qq! SELECT a.id, a.clone_id, a.query, a.query_start, a.query_end,
                            a.query_remaining, a.hit, a.hit_start, a.hit_end,
                            a.hit_remaining, a.comp, d.fs_contig, d.fs_position,
                            d.fs_base, d.fs_qual, d.fn_position, d.fn_base,
                            d.fn_qual, a.expt_id, d.confirm
                     FROM nhansen.fs_alignment a, nhansen.fs_discrepancy d
                     WHERE  a.id = d.fs_alignment_id (+)
                     AND    a.expt_id = $expt_id !;

    if (defined ($clone_id))
    {
        $select .= qq! AND a.clone_id = $clone_id !;
    }

    my $sth = $dbh->prepare($select)
         or croak "Couldn\'t prepare statement $select!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $select!\n";

    my %discreps;
    my %aligns;
    while (my $ra_align = $sth->fetchrow_arrayref())
    {
        my ($id, $clone_id, $query, $query_start, $query_end,
             $query_remaining, $hit, $hit_start, $hit_end, $hit_remaining,
             $comp, $fs_contig, $fs_position, $fs_base, $fs_qual, $fn_position,
             $fn_base, $fn_qual, $expt_id, $confirm) = @{$ra_align};

        if (!defined ($aligns{$id}))
        {
            $aligns{$id} = $class->new(-query => $query, -query_start => $query_start,
                             -query_end => $query_end, -query_remaining => $query_remaining,
                             -hit => $hit, -hit_start => $hit_start, -hit_end => $hit_end,
                             -hit_remaining => $hit_remaining, -comp => $comp,
                             -expt_id => $expt_id);
        }

        if ($fs_position)
        {
             my $discrep = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Comp::FSDiscrepancy->new(-clone_id => $clone_id, -fs_contig => $fs_contig,
                             -fs_position => $fs_position, -fs_base => $fs_base,
                             -fs_qual => $fs_qual, -fn_position => $fn_position,
                             -fn_base => $fn_base, -fn_qual => $fn_qual, 
                             -confirm => $confirm);

             push @{$discreps{$id}}, $discrep;
        }
    }

    my @alignments = ();
    foreach my $id (keys %aligns)
    {
         my $ra_discreps = $discreps{$id} || [];
         $aligns{$id}->discrepancies($ra_discreps);
         push @alignments, $aligns{$id};
    }

    return \@alignments;

} ## end db_fs_select

###############################################################################
# Subroutine to retrieve Alignment objects from the niscprod sanger_alignment 
#     table.
#
# INPUT: Class, -clone_id specifies clone_id value,
#        -expt_id specifies which dataset 
#        -dbh passes database handle (not optional yet!)
# OUTPUT: reference to an array of NISC::Assembly::Alignment objects
###############################################################################

sub db_sanger_select {

    my $class = shift;
    my %params = @_;
    my $clone_id = $params{'-clone_id'};
    my $expt_id = $params{'-expt_id'};
    my $dbh = $params{'-dbh'};

    my $select = qq! SELECT a.id, a.clone_id, a.query, a.query_start, a.query_end,
                            a.query_remaining, a.hit, a.hit_start, a.hit_end,
                            a.hit_remaining, a.comp, d.fs_contig, d.fs_position,
                            d.fs_base, d.fs_qual, d.fn_position, d.fn_base,
                            d.fn_qual, a.expt_id
                     FROM nhansen.sanger_alignment a, nhansen.sanger_discrepancy d
                     WHERE  a.id = d.sanger_alignment_id (+)
                     AND    a.clone_id = $clone_id 
                     AND    a.expt_id = $expt_id !;

    my $sth = $dbh->prepare($select)
         or croak "Couldn\'t prepare statement $select!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $select!\n";

    my %discreps;
    my %aligns;
    while (my $ra_align = $sth->fetchrow_arrayref())
    {
        my ($id, $clone_id, $query, $query_start, $query_end,
             $query_remaining, $hit, $hit_start, $hit_end, $hit_remaining,
             $comp, $fs_contig, $fs_position, $fs_base, $fs_qual, $fn_position,
             $fn_base, $fn_qual, $expt_id) = @{$ra_align};

        if (!defined ($aligns{$id}))
        {
            $aligns{$id} = $class->new(-query => $query, -query_start => $query_start,
                             -query_end => $query_end, -query_remaining => $query_remaining,
                             -hit => $hit, -hit_start => $hit_start, -hit_end => $hit_end,
                             -hit_remaining => $hit_remaining, -comp => $comp,
                             -expt_id => $expt_id);
        }

        if ($fs_position)
        {
             my $discrep = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Comp::SangerDiscrepancy->new(-clone_id => $clone_id, -fs_contig => $fs_contig,
                             -fs_position => $fs_position, -fs_base => $fs_base,
                             -fs_qual => $fs_qual, -fn_position => $fn_position,
                             -fn_base => $fn_base, -fn_qual => $fn_qual);

             push @{$discreps{$id}}, $discrep;
        }
    }

    my @alignments = ();
    foreach my $id (keys %aligns)
    {
         my $ra_discreps = $discreps{$id} || [];
         $aligns{$id}->discrepancies($ra_discreps);
         push @alignments, $aligns{$id};
    }

    return \@alignments;

} ## end db_sanger_select

###############################################################################
# Subroutine to retrieve Alignment objects from the niscprod wgs_alignment 
#     table.
#
# INPUT: Class, -clone_id specifies clone_id value,
#        -expt_id specifies which dataset 
#        -dbh passes database handle (not optional yet!)
# OUTPUT: reference to an array of NISC::Assembly::Alignment objects
###############################################################################

sub db_wgs_select {

    my $class = shift;
    my %params = @_;
    my $species = $params{'-species'};
    my $expt_id = $params{'-expt_id'} || 1;
    my $dbh = $params{'-dbh'};

    my $select = qq! SELECT a.id, a.species, a.query, a.query_start, a.query_end,
                            a.query_remaining, a.hit, a.hit_start, a.hit_end,
                            a.hit_remaining, a.comp, d.fs_contig, d.fs_position,
                            d.fs_base, d.fs_qual, d.fn_position, d.fn_base,
                            d.fn_qual, a.expt_id
                     FROM nhansen.wgs_alignment a, nhansen.wgs_discrepancy d
                     WHERE  a.id = d.wgs_alignment_id (+)
                     AND    a.species = '$species'
                     AND    a.expt_id = $expt_id !;

    my $sth = $dbh->prepare($select)
         or croak "Couldn\'t prepare statement $select!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $select!\n";

    my %discreps;
    my %aligns;
    while (my $ra_align = $sth->fetchrow_arrayref())
    {
        my ($id, $clone_id, $query, $query_start, $query_end,
             $query_remaining, $hit, $hit_start, $hit_end, $hit_remaining,
             $comp, $fs_contig, $fs_position, $fs_base, $fs_qual, $fn_position,
             $fn_base, $fn_qual, $expt_id) = @{$ra_align};

        if (!defined ($aligns{$id}))
        {
            $aligns{$id} = $class->new(-query => $query, -query_start => $query_start,
                             -query_end => $query_end, -query_remaining => $query_remaining,
                             -hit => $hit, -hit_start => $hit_start, -hit_end => $hit_end,
                             -hit_remaining => $hit_remaining, -comp => $comp,
                             -expt_id => $expt_id);
        }

        if ($fs_position)
        {
             my $discrep = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Comp::WGSDiscrepancy->new(-species => $clone_id,
                             -fs_position => $fs_position, -fs_base => $fs_base,
                             -fs_qual => $fs_qual, -fn_position => $fn_position,
                             -fn_base => $fn_base, -fn_qual => $fn_qual);

             push @{$discreps{$id}}, $discrep;
        }
    }

    my @alignments = ();
    foreach my $id (keys %aligns)
    {
         my $ra_discreps = $discreps{$id} || [];
         $aligns{$id}->discrepancies($ra_discreps);
         push @alignments, $aligns{$id};
    }

    return \@alignments;

} ## end db_wgs_select

1;

__END__

=head1 NAME

NISC::Assembly::Alignment - Perl extension for containing and manipulating information about an alignment.

=head1 SYNOPSIS

use NISC::Assembly::Alignment;

my $align = NISC::Assembly::Alignment->new(-align_type => 'crossmatch', -comp => 'U', ...);

my $align_type = $align->align_type();

=head1 DESCRIPTION

This module is used to store and retrieve information about sequence alignments.

=head1 AUTHOR

Nancy F. Hansen <nhansen@nhgri.nih.gov>

=head1 SEE ALSO

perl(1).

=cut

