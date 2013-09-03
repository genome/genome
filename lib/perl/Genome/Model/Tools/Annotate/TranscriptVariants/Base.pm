package Genome::Model::Tools::Annotate::TranscriptVariants::Base;

use strict;
use warnings;

use Data::Dumper;
use Genome;
use File::Temp;
use List::Util qw/ max min /;
use List::MoreUtils qw/ uniq /;
use Bio::Seq;
use Bio::Tools::CodonTable;
use DateTime;
use Carp;

class Genome::Model::Tools::Annotate::TranscriptVariants::Base {
    has => [
        codon_translator => {
            is => 'Bio::Tools::CodonTable',
            is_constant => 1,
            calculate => q( Bio::Tools::CodonTable->new( -id => 1) ),
        },

        mitochondrial_codon_translator => {
            is => 'Bio::Tools::CodonTable',
            is_constant => 1,
            calculate => q( Bio::Tools::CodonTable->new( -id => 2) ),
        },

        check_variants => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, the reference sequence on variants is checked against our reference',
        },

        get_frame_shift_sequence => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, the entire modifed sequence of the transcript is placed in the output file for frameshift mutations, even if the modification is silent',
        },

        data_directory => {
            is => 'PATH',
            is_optional => 0,
            doc => 'Pathname to the annotation_data of a build containing transcripts.csv and other files',
        },

        transcript_structure_class_name => {
            is_constant => 1,
            value => 'Genome::TranscriptStructure',
        },

        #priorities => { is => __PACKAGE__ . '::AnnotationPriorities', is_constant => 1, id_by => 1 },
        transcript_source_priorities => {  },
        transcript_status_priorities => {  },
        variant_priorities           => {  },
        transcript_error_priorities  => {  },

    ],
};

## These originally lived in Genome::Info::AnnotationPriorities
sub transcript_source_priorities {
    return (
        genbank => 1,
        ensembl => 2,
    );
}

sub variant_priorities  {
    return (
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
}

sub transcript_error_priorities {
    return (
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
}

sub _resolve_intersector_sub_name {
    my ($self) = @_;

    my $structure_class = $self->transcript_structure_class_name;
    return $structure_class->__meta__->data_source_id . '::intersector_sub';
}


# Reverse complement the given sequence
sub reverse_complement {
    my ($self, $seq) = @_;
    return unless defined $seq;

    $seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $rev_com = CORE::reverse $seq;
    return $rev_com;
}

# Prioritizes a list of annotations and returns the highest priority annotation per gene
# I.E. If 6 annotations go in, from 3 different genes, it will select the "best" annotation
# that each gene has, and return 3 total annotations, one per gene
# See the _prioritize_annotations method for details of sorting
# Corresponds to gene filter in Genome::Model::Tools::Annotate::TranscriptVariants
sub prioritized_transcripts {
    my ($self, %variant) = @_;
    my @annotations = $self->transcripts(%variant) or return;

    my @prioritized_annotations = $self->_prioritize_annotations(@annotations);
    my %gene_annotations;
    for my $annotation (@prioritized_annotations) {
        unless (exists $gene_annotations{$annotation->{gene_name}}) {
            $gene_annotations{$annotation->{gene_name}} = $annotation;
        }
    }
    return values %gene_annotations;
}

# Returns the annotation with highest priority
# Does not sort the full list of annotations, just does a single pass to find the highest priority annotation,
# which should have run time O(n) instead of O(nlogn)
# See the _prioritize_annotations method for details for sorting
# Corresponds to top filter in Genome::Model::Tools::Annotate::TranscriptVariants
sub prioritized_transcript{
    my ($self, %variant) = @_;
    my @annotations = $self->transcripts(%variant) or return;

    my $highest = shift @annotations;
    for my $anno (@annotations) {
        my @sorted_list = $self->_prioritize_annotations($highest, $anno);
        $highest = shift @sorted_list;
    }
    return $highest;
}

# Find which substructures intersect the given variant as quickly as possible.
# $self holds an object iterator for the current chromosome and a cache of
# TranscriptStructures it has read from the iterator.
# It does some "bad" things like poke directly into the object's hash instead of
# going through the accessors
use constant TRANSCRIPT_STRUCTURE_ID => 0;
use constant STRUCTURE_START => 3;
use constant STRUCTURE_STOP => 4;
sub _create_iterator_for_variant_intersection {
    my $self = shift;

    my $structure_iterator;
    my $loaded_substructures = [];
    my $last_variant_start = -1;
    my $last_variant_stop = -1;
    my $last_chrom_name = '';
    my $next_substructure;

    my $variant;   # needs to be visible in both closures below

    my $structure_class = $self->transcript_structure_class_name;
    my $intersect_sub_name = $self->_resolve_intersector_sub_name();

    #Because of this hack, the query cache is not actually valid.  The cache remembers that it did
    #a query, but it didn't actually cache everything.  To fix this, delete the template for the
    #query rule each time you do the query.
    my $template_id = UR::BoolExpr::Template::And->_fast_construct($structure_class, ['chrom_name', 'data_directory', '-order_by'], ['structure_start'])->get_normalized_template_equivalent->id;

    # This sub plugs into a hook in the Genome::DataSource::TranscriptStructures loader
    # to reject data that does not intersect the given variation to avoid passing the
    # data up the call stack and creating objects for TranscriptStructures we aren't
    # interested in
    {   no strict 'refs';
        $$intersect_sub_name = sub {
            return 1 unless defined $variant;

            my $struct = $_[0];
            if ( $variant->{'start'} <= $struct->[STRUCTURE_STOP]) {
                return 1;
            } else {
                return 0;
            }
        };
    }

    return sub {
        $variant = $_[0];

        my $variant_start = $variant->{'start'};
        my $variant_stop  = $variant->{'stop'};

        if ($variant->{'chromosome_name'} ne $last_chrom_name
            or
            $variant_start < $last_variant_start
        ) {
            my $chrom_name = $variant->{'chromosome_name'};
            $self->status_message("Resetting iterator for chromosome $chrom_name");

            $loaded_substructures = [];
            $structure_iterator = undef;

            $structure_iterator = $structure_class->create_iterator(
                                      chrom_name => $chrom_name,
                                      data_directory => $self->data_directory,
                                      -order_by => ['structure_start']);
            $next_substructure = $structure_iterator->next;
        }

        my @intersections;

        # Keep structures that cross or are after the variant
        my @keep_structures;
        foreach my $substr ( @$loaded_substructures ) {
            if ($variant_start <= $substr->{'structure_stop'}) {
                push @keep_structures, $substr;
                if ($variant_stop >= $substr->{'structure_start'}) {
                    push @intersections, $substr;
                }
            }
        }

        # fast-forward through to substructures that cross
        while($next_substructure and $next_substructure->{'structure_start'} <= $variant_stop) {
            push @keep_structures, $next_substructure;
            if ($variant_start <= $next_substructure->{'structure_stop'}) {
                push @intersections, $next_substructure;
            }
            $next_substructure = $structure_iterator->next();
        }
        $loaded_substructures = \@keep_structures;

        $last_variant_start = $variant_start;
        $last_variant_stop  = $variant_stop;
        $last_chrom_name    = $variant->{'chromosome_name'};

        delete $UR::Context::all_params_loaded->{$template_id};
        return \@intersections;
    };
}

# Checks that the sequence on the variant is valid
sub is_valid_variant {
    my ($self, $variant) = @_;
    unless ($variant->{type} eq 'DEL') {
        return 0 if $variant->{variant} =~ /\d/;
    }

    unless ($variant->{type} eq 'INS') {
        return 0 if $variant->{reference} =~ /\d/;
    }
    return 1;
}

# Takes in a group of annotations and returns them in priority order
# Sorts on: variant type, transcript error, source, status, amino acid length, transcript name (descending importance)
# Rankings for each category can be found in Genome::Info::AnnotationPriorities
sub _prioritize_annotations {
    my ($self, @annotations) = @_;

    my %transcript_source_priorities = $self->transcript_source_priorities;
    my %transcript_status_priorities = $self->transcript_status_priorities;
    my %variant_priorities = $self->variant_priorities;

    use sort '_mergesort';  # According to perldoc, performs better than quicksort on large sets with many comparisons
    my @sorted_annotations = sort {
        $variant_priorities{$a->{trv_type}} <=> $variant_priorities{$b->{trv_type}} ||
        $self->_highest_priority_error($a) <=> $self->_highest_priority_error($b) ||
        $transcript_source_priorities{$a->{transcript_source}} <=> $transcript_source_priorities{$b->{transcript_source}} ||
        $transcript_status_priorities{$a->{transcript_status}} <=> $transcript_status_priorities{$b->{transcript_status}} ||
        $b->{amino_acid_length} <=> $a->{amino_acid_length} ||
        $a->{transcript_name} cmp $b->{transcript_name}
    } @annotations;
    use sort 'defaults';

    return @sorted_annotations;
}

# Given an annotation, split up the error string and return the highest priority error listed
sub _highest_priority_error {
    my ($self, $annotation) = @_;

    my %transcript_error_priorities = $self->transcript_error_priorities;

    my $error_string = $annotation->{transcript_error};
    my @errors = map { $transcript_error_priorities{$_} } split(":", $error_string);
    my @sorted_errors = sort { $b <=> $a } @errors;
    return $sorted_errors[0];
}

sub _transcript_annotation_for_rna {
    my ($self, $variant, $structure) = @_;

    return (
        c_position => 'NULL',
        trv_type => 'rna',
        amino_acid_length => 'NULL',
        amino_acid_change => 'NULL',
    );
}

sub _transcript_annotation_for_utr_exon {
    my ($self, $variant, $structure) = @_;
    my $coding_position = $self->_determine_coding_position($variant, $structure);

    # TODO Change to use variant start and stop for more accurate annotation
    # TODO Not sure if the behavior when there is no coding region is ideal
    my $position = $variant->{start};
    my $trv_type;
    if ($structure->transcript_has_coding_region) {
        if ($structure->transcript_before_coding_region($position)) {
            $trv_type = '5_prime_untranslated_region';
        }
        else {
            $trv_type = '3_prime_untranslated_region';
        }
    }
    else {
        if ($structure->transcript_strand eq '+1') {
            $trv_type = '3_prime_untranslated_region';
        }
        else {
            $trv_type = '5_prime_untranslated_region';
        }
    }

    return (
        c_position => "c." . $coding_position,
        trv_type => $trv_type,
        amino_acid_change => 'NULL',
    );
}

sub _transcript_annotation_for_flank {
    my ($self, $variant, $structure) = @_;
    my $coding_position = $self->_determine_coding_position($variant, $structure);

    # TODO Change to use variant start and stop
    my $position = $variant->{start};
    my $trv_type;
    if ($structure->transcript_is_before($variant->{start}, $structure->transcript_transcript_start)) {
        $trv_type = '5_prime_flanking_region';
    }
    else {
        $trv_type = '3_prime_flanking_region';
    }

    # TODO Change to use variant start and stop for more accurate distance
    my $distance_to_transcript = $structure->transcript_distance_to_transcript($position);

    return (
        c_position => "c." . $coding_position,
        trv_type => $trv_type,
        amino_acid_change => 'NULL',
        flank_annotation_distance_to_transcript => $distance_to_transcript,
    );
}

1;
