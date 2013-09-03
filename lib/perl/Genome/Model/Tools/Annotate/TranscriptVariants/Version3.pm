package Genome::Model::Tools::Annotate::TranscriptVariants::Version3;

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

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Genome::Model::Tools::Annotate::TranscriptVariants::Base',

    has => [
        reference_sequence_id => {
            is => "Text",
        },
        eids => {
            is_transient => 1,
            is_optional => 1,
            doc => "Temporary variable used for intermediate calculation.",
        },
    ],

    doc => q(Do proper intersections between variations and transcript structures by considering both entities' start and stop positions rather than just the start position of the variation.),

);


sub transcript_status_priorities {
    return (
        reviewed            => 1,
        validated           => 2,
        provisional         => 3,
        predicted           => 4,
        putative            => 4,
        model               => 5,
        inferred            => 6,
        known               => 7,
        annotated           => 8,
        known_by_projection => 9,
        novel               => 10,
        unknown             => 11,
    );
}

# Given a nucleotide sequence, translate to amino acid sequence and return
# The translator->translate method can take 1-3 bp of sequence. If given
# 1 bp, an empty string is returned. If given 2 bp, it attempts to translate
# this ambiguous sequence and returns an empty string if a no translation
# can be made. Translation works as you'd expect with 3bp.
sub translate {
    my ($self, $seq, $chrom_name) = @_;
    return unless defined $seq and defined $chrom_name;

    my $is_mitochondrial = $chrom_name =~ /^MT?/; #we use the mitochondrial_codon_translator if the chromosome is either the M or MT.  Everything else should use the normal translator

    my $translator;
    if ($is_mitochondrial) {
        $translator = $self->mitochondrial_codon_translator;
    }
    else {
        $translator = $self->codon_translator;
    }

    my $translation;
    my $length = length $seq;
    for (my $i=0; $i<=$length-2; $i+=3) {
        my $codon=substr($seq, $i, 3);
        $codon =~ s/N/X/g;
        my $aa = $translator->translate($codon);
        $aa="*" if ($aa eq 'X');
        $translation.=$aa;
    }
    return $translation;
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


# Annotates all transcripts affected by a variant
# Corresponds to none filter in Genome::Model::Tools::Annotate::TranscriptVariants
sub transcripts {
    my ($self, %variant) = @_;

    if (!defined $self->{_cached_chromosome} or $self->{_cached_chromosome} ne $variant{chromosome_name}) {
        Genome::InterproResult->unload();
        $self->transcript_structure_class_name->unload();
        
        Genome::InterproResult->get(
            data_directory => $self->data_directory,
            chrom_name => $variant{chromosome_name},
        );

        if (!defined $self->{_cached_chromosome}) {
            Genome::ExternalGeneId->get(
                data_directory => $self->data_directory,
                reference_build_id => $self->reference_sequence_id,
                id_type => "ensembl",
            );
            Genome::ExternalGeneId->get(
                data_directory => $self->data_directory,
                reference_build_id => $self->reference_sequence_id,
                id_type => "ensembl_default_external_name",
            );
            Genome::ExternalGeneId->get(
                data_directory => $self->data_directory,
                reference_build_id => $self->reference_sequence_id,
                id_type => "ensembl_default_external_name_db",
            );
        }
        $self->{_cached_chromosome} = $variant{chromosome_name};
    }

    my $variant_start = $variant{'start'};
    my $variant_stop = $variant{'stop'};

    $variant{'type'} = uc($variant{'type'});

    # Make sure variant is set properly
    unless (defined($variant_start) and defined($variant_stop) and defined($variant{variant})
            and defined($variant{reference}) and defined($variant{type})
            and defined($variant{chromosome_name}))
    {
        print Dumper(\%variant);
        confess "Variant is not fully defined: chromosome name, start, stop, variant, reference, and type must be defined.\n";
    }

    # Make sure the sequence on the variant is valid. If not, display a warning and set transcript error
    unless ($self->is_valid_variant(\%variant)) {
        $self->warning_message("The sequence on this variant is not valid! Reference: $variant{reference} Variant: $variant{variant}");

        return { transcript_error => 'invalid_sequence_on_variant' }
    }

    my $windowing_iterator = $self->{'_windowing_iterator'};
    unless ($windowing_iterator) {
        $windowing_iterator = $self->{'_windowing_iterator'} = $self->_create_iterator_for_variant_intersection();
    }
    my $crossing_substructures = $windowing_iterator->(\%variant);

    # Hack to support the old behavior of only annotating against the first structure
    # of a transcript.  We need to keep a list of all the other structures for later
    # listing them in the deletions column of the output
#    my %transcript_substructures;
#    {
#        my @less;
#        foreach my $substructure ( @$crossing_substructures ) {
#            my $transcript_id = $substructure->transcript_transcript_id;
#            if ($substructure->{'structure_start'} <= $variant_start and $substructure->{'structure_stop'} >= $variant_start) {
#                push @less, $substructure;
#            }
#            $transcript_substructures{$transcript_id} ||= [];
#            push @{$transcript_substructures{$transcript_id}}, $substructure;
#        }
#        $crossing_substructures = \@less;
#    }

    return unless @$crossing_substructures;

    my @annotations;
    my $variant_checked = 0;

    Genome::Model::Tools::Annotate::LookupConservationScore->class(); # get it loaded so we can call lookup_sequence

    foreach my $substruct ( @$crossing_substructures ) {
        # If specified, check that the reference sequence stored on the variant correctly matches our reference
        if ($self->check_variants and not $variant_checked) {
            unless ($variant{reference} eq '-') {
                my $chrom = $variant{chromosome_name};
                my $species = $substruct->transcript_species;
                my $ref_seq = Genome::Model::Tools::Sequence::lookup_sequence(
                                  chromosome => $chrom,
                                  start => $variant_start,
                                  stop => $variant_stop,
                                  species => $species,
                              );

                unless ($ref_seq eq $variant{reference}) {
                    $self->warning_message("Sequence on variant on chromosome $chrom between $variant_start and $variant_stop does not match $species reference!");
                    $self->warning_message("Variant sequence : " . $variant{reference});
                    $self->warning_message("$species reference : " . $ref_seq);
                    return;
                }
                $variant_checked = 1;
            }
        }
        
        my %annotation = $self->_transcript_substruct_annotation($substruct, %variant) or next;

#        # Continuation of the hack above about annotating a deletion
#        if ($variant{'type'} eq 'DEL') {
#            my @del_strings = map { $_->structure_type . '[' . $_->structure_start . ',' . $_->structure_stop . ']' }
#                                  @{$transcript_substructures{$substruct->transcript_transcript_id}};
#            $annotation{'deletion_substructures'} = '(deletion:' . join(', ', @del_strings) . ')';
#        }
        push @annotations, \%annotation;
    }
    return @annotations;
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

# Annotates a single transcript-substructure/variant pair
sub _transcript_substruct_annotation {
    my ($self, $substruct, %variant) = @_;
    # Just an FYI... using a copy of variant here instead of a reference prevents reverse complementing
    # the variant twice, which would occur if the variant happened to touch two reverse stranded transcripts

    # TODO This will need to be coupled with splitting up a variant so it doesn't extend beyond a structure
    # TODO There are various hacks in intron and exon annotation to fix the side effects of only annotating the
    # structure at variant start position that will also need removed once this is fixed, and there is still
    # a definite bias for variant start over variant stop throughout this module

#    # If the variant extends beyond the current substructure, it needs to be resized
#    if ($variant{start} < $substruct->{structure_start}) {
#       my $diff = $substruct->{structure_start} - $variant{start};
#       $variant{start} = $variant{start} + $diff;
#       unless ($variant{type} eq 'DEL') {
#           $variant{variant} = substr($variant{variant}, $diff);
#       }
#       unless ($variant{type} eq 'INS') {
#           $variant{reference} = substr($variant{reference}, $diff);
#       }
#    }
#    elsif ($variant{stop} > $substruct->{structure_stop}) {
#        my $diff = $variant{stop} - $substruct->{structure_stop};
#        $variant{stop} = $variant{stop} - $diff;
#        unless ($variant{type} eq 'DEL') {
#            $variant{variant} = substr($variant{variant}, 0, length($variant{variant}) - $diff);
#        }
#        unless ($variant{type} eq 'INS') {
#            $variant{reference} = substr($variant{reference}, 0, length($variant{reference}) - $diff);
#        }
#    }

    # All sequence stored on the variant is forward stranded and needs to be reverse
    # complemented if the transcript is reverse stranded.
    my $strand = $substruct->transcript_strand;
    if ($strand eq '-1') {
        my ($new_variant, $new_reference);
        unless ($variant{type} eq 'DEL') {
            $new_variant = $self->reverse_complement($variant{variant});
            $variant{variant} = $new_variant;
        }
        unless ($variant{type} eq 'INS') {
            $new_reference = $self->reverse_complement($variant{reference});
            $variant{reference} = $new_reference;
        }
    }

    my $structure_type = $substruct->structure_type;
    my $method = '_transcript_annotation_for_' . $structure_type;
    my %structure_annotation = $self->$method(\%variant, $substruct) or return;
    
    my $conservation = $self->_ucsc_conservation_score(\%variant, $substruct);

    my $gene_name = $substruct->transcript_gene_name;
    my $dumper_string = $substruct->id;
    unless ($gene_name) {
        $self->warning_message("Gene name missing for substruct: $dumper_string");
        my $gene = Genome::Gene->get(data_directory => $substruct->data_directory,
                                     id => $substruct->transcript_gene_id,
                                     reference_build_id => $self->reference_sequence_id);
        $gene_name = $gene->name;
    }

    my ($default_gene_name, $ensembl_gene_id, $gene_name_source);
    unless ($self->eids) {
        my %new;
        $self->eids(\%new);
    }
    if ($self->eids and $self->eids->{$substruct->transcript_gene_id}) {
        ($default_gene_name,$ensembl_gene_id,$gene_name_source) = split(/,/, $self->eids->{$substruct->transcript_gene_id});
    }
    else {
        my @e1 = Genome::ExternalGeneId->get(data_directory => $substruct->data_directory,
            gene_id => $substruct->transcript_gene_id,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl_default_external_name");
        if ($e1[0]) {
            $default_gene_name = $e1[0]->id_value;
        }
        unless ($default_gene_name) {
            $self->warning_message("Ensembl gene name missing for substruct: $dumper_string");
            $default_gene_name = "Unknown";
        }
        my @e2 = Genome::ExternalGeneId->get(data_directory => $substruct->data_directory,
            gene_id => $substruct->transcript_gene_id,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl");
        if ($e2[0]) {
            $ensembl_gene_id = $e2[0]->id_value;
        }
        unless ($ensembl_gene_id) {
            $self->warning_message("Ensembl stable gene id missing for substruct: $dumper_string");
            $ensembl_gene_id = "Unknown";
        }
        my @e3 = Genome::ExternalGeneId->get(data_directory => $substruct->data_directory,
            gene_id => $substruct->transcript_gene_id,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl_default_external_name_db");
        if ($e3[0]) {
            $gene_name_source = $e3[0]->id_value;
        }
        unless ($gene_name_source) {
            $self->warning_message("Ensembl gene name source missing for substruct: $dumper_string");
            $gene_name_source = "Unknown";
        }
        $self->eids->{$substruct->transcript_gene_id} = join(',',$default_gene_name, $ensembl_gene_id, $gene_name_source);
    }

    return (
        %structure_annotation,
        transcript_error => $substruct->transcript_transcript_error,
        transcript_name => $substruct->transcript_transcript_name,
        transcript_status => $substruct->transcript_transcript_status,
        transcript_source => $substruct->transcript_source,
        transcript_species=> $substruct->transcript_species,
        transcript_version => $substruct->transcript_version,
        strand => $strand,
        gene_name  => $gene_name,
        amino_acid_length => $substruct->transcript_amino_acid_length,
        ucsc_cons => $conservation,
        default_gene_name => $default_gene_name,
        gene_name_source => $gene_name_source,
        ensembl_gene_id => $ensembl_gene_id,
    );
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

sub _transcript_annotation_for_intron {
    my ($self, $variant, $structure) = @_;
    my $coding_position = $self->_determine_coding_position($variant, $structure);

    # Need shortest distance between variant and the edge of the next structure (one bp past edge of intron)
    my $dist_to_structure_start = $variant->{start} - ($structure->{structure_start} - 1);
    my $dist_to_structure_stop = ($structure->{structure_stop} + 1) - $variant->{stop};
    my $distance_to_edge = $dist_to_structure_start <= $dist_to_structure_stop ? 
        $dist_to_structure_start : $dist_to_structure_stop;

    # Number of basepairs between start of intron and variant (takes into account strand)
    my $intron_position;
    if ($structure->transcript_strand eq '+1') {
        $intron_position = $variant->{start} - $structure->{structure_start} + 1;
    }
    else {
        $intron_position = $structure->{structure_stop} - $variant->{stop} + 1;
    }

    # If variant occurs within 2bp of a neighboring structure, it's splice site
    # If variant occurs within 3-10bp, it's splice region
    # Otherwise, it's intronic
    my $trv_type;
    if ($variant->{type} eq 'INS' or $variant->{'type'} eq 'DEL') {
        $trv_type = "intronic";
        $trv_type = "splice_region_" . lc $variant->{type} if $distance_to_edge <= 10;
        $trv_type = "splice_site_" . lc $variant->{type} if $distance_to_edge <= 2;
    }
    else {
        $trv_type = "intronic";
        $trv_type = "splice_region" if $distance_to_edge <= 10;
        $trv_type = "splice_site" if $distance_to_edge <= 2;
    }

    # If splice site or splice region variant, include ordinal of neighboring exon and how many base pairs away
    # if there are any coding exons in the transcript
    my $amino_acid_change = 'NULL';
    my $has_exons = ($structure->{cds_exons_before} != 0 or $structure->{cds_exons_after} != 0);
    if ($has_exons) {
        # For indels, the distance to the edge of the bordering exon can be negative if the indel starts in the intron
        # and carries over to the exon. This should eventually be fixed by splitting up variants in the _transcript_annotation
        # method above. For now, set the distance to edge to 1 if its negative
        # This is actually a pretty big bug, since currently the only structure annotated is the one at the variant's start. This
        # needs to be changed such that all structures within the range of the variant are annotated.
        # TODO Remove once variants spanning structures are handled correctly
        $distance_to_edge = 1 if $distance_to_edge <= 0;

        # If position in the intron is equal to the distance to the closest edge, 
        # this splice site is just after an exon
        if ($intron_position eq $distance_to_edge) {
            $amino_acid_change = "e" . $structure->{cds_exons_before} . "+" . $distance_to_edge;
        }
        else {
            $amino_acid_change = "e" . ($structure->{cds_exons_before} + 1) . "-" . $distance_to_edge;
        }
    }

    return (
        c_position => "c." . $coding_position,
        trv_type => $trv_type,
        amino_acid_change => $amino_acid_change,
        intron_annotation_substructure_ordinal => $structure->ordinal,
        intron_annotation_substructure_size => $structure->length,
        intron_annotation_substructure_position => $intron_position
    );
}

sub _transcript_annotation_for_cds_exon {
    my ($self, $variant, $structure) = @_;
    

    # If the variant continues beyond the stop position of the exon, then the variant sequence
    # needs to be modified to stop at the exon's stop position. The changes after the exon's stop
    # affect the intron, not the coding sequence, and shouldn't be annotated here. Eventually,
    # it's possible that variants may always be split up so they only touch one structure, but for 
    # now this will have to do.
    # TODO This can be removed once variants spanning structures are handled properly
    unless ($self->{'get_frame_shift_sequence'}) {
        # If we're inspecting the entire sequence, don't chop the variant down...
        if ($variant->{stop} > $structure->structure_stop and $variant->{type} eq 'DEL') {
            my $bases_beyond_stop = $variant->{stop} - $structure->structure_stop;
            my $new_variant_length = (length $variant->{reference}) - $bases_beyond_stop;
            $variant->{reference} = substr($variant->{reference}, 0, $new_variant_length);
            $variant->{stop} = $variant->{stop} - $bases_beyond_stop;
        }
    }

    my $coding_position = $self->_determine_coding_position($variant, $structure);

    # Grab and translate the codons affected by the variation
    my ($original_seq, $mutated_seq, $protein_position) = $self->_get_affected_sequence($structure, $variant);
    my $chrom_name = $structure->transcript_chrom_name;
    my $original_aa = $self->translate($original_seq, $chrom_name);
    my $mutated_aa = $self->translate($mutated_seq, $chrom_name);

    my ($trv_type, $protein_string);
    my ($reduced_original_aa, $reduced_mutated_aa, $offset) = ($original_aa, $mutated_aa, 0);
    if ($variant->{type} eq 'INS') {
        ($reduced_original_aa, $reduced_mutated_aa, $offset) = $self->_reduce($original_aa, $mutated_aa);
        $protein_position += $offset;

        my $indel_size = length $variant->{variant};
        if ($indel_size % 3 == 0) {
            $trv_type = "in_frame_" . lc $variant->{type};
            $protein_string = "p." . $reduced_original_aa . $protein_position . $trv_type . $reduced_mutated_aa;
        }
        else {
            # In this case, the inserted sequence does not change the amino acids on either side of it (if original is
            # MT and mutated is MRPT, then original sequence is reduced to nothing). The first changed amino acid should
            # be set to the first original amino acid that does not occur in the mutated sequence moving 5' to 3'. 
            # In the above example, first changed amino acid would be T.
            $trv_type = "frame_shift_" . lc $variant->{type};
            if ($self->{get_frame_shift_sequence}) {
                my $aa_after_indel = $self->_apply_indel_and_translate($structure, $variant);
                $protein_string = "p." . $aa_after_indel . $protein_position . "fs";
            }
            else {
                if ($reduced_original_aa eq "") {
                    $protein_position -= $offset;
                    for (my $i = 0; $i < (length $original_aa); $i++) {
                        my $original = substr($original_aa, $i, 1);
                        my $mutated = substr($mutated_aa, $i, 1);
                        $protein_position++ and next if $original eq $mutated;
                        $reduced_original_aa = $original;
                        last;
                    }
                    
                    # If the original sequence is STILL reduced to nothing (insertion could occur after original amino acids),
                    # then just use the last amino acid in the unreduced original sequence
                    $reduced_original_aa = substr($original_aa, -1) if $reduced_original_aa eq "";
                }

                $reduced_original_aa = substr($reduced_original_aa, 0, 1);
                $protein_string = "p." . $reduced_original_aa . $protein_position . "fs";
            }
        }
    }
    elsif ($variant->{type} eq 'DEL') {
        ($reduced_original_aa, $reduced_mutated_aa, $offset) = $self->_reduce($original_aa, $mutated_aa);
        $protein_position += $offset;

        my $indel_size = length $variant->{reference};
        if ($indel_size % 3 == 0) {
            $trv_type = "in_frame_" . lc $variant->{type};
            $protein_string = "p." . $reduced_original_aa . $protein_position . $trv_type;
            $protein_string .= $reduced_mutated_aa if defined $reduced_mutated_aa;
        }
        else {
            $trv_type = "frame_shift_" . lc $variant->{type};
            if ($self->{get_frame_shift_sequence}) {
                my $aa_after_indel = $self->_apply_indel_and_translate($structure, $variant);
                $protein_string = "p." . $aa_after_indel . $protein_position . "fs";
            }
            else {
                $reduced_original_aa = substr($reduced_original_aa, 0, 1); 
                $protein_string = "p." . $reduced_original_aa . $protein_position . "fs";
            }
        }
    }
    elsif ($variant->{type} eq 'DNP' or $variant->{type} eq 'SNP') {
        if (!defined $mutated_aa or !defined $original_aa) {
            $trv_type = 'silent';
            $protein_string = "NULL";
        }
        elsif ($mutated_aa eq $original_aa) {
            $trv_type = 'silent';
            $protein_string = "p." . $original_aa . $protein_position;
        }
        else {
            ($reduced_original_aa, $reduced_mutated_aa, $offset) = $self->_reduce($original_aa, $mutated_aa);
            $protein_position += $offset;

            if (index($reduced_mutated_aa, '*') != -1) {
                $trv_type = 'nonsense';
            }
            elsif (index($reduced_original_aa, '*') != -1) {
                $trv_type = 'nonstop';
            }
            else {
                $trv_type = 'missense';
            }
            $protein_string = "p." . $reduced_original_aa . $protein_position . $reduced_mutated_aa;
        }
    }
    else {
        $self->warning_message("Unknown variant type " . $variant->{type} .
                " for variant between " . $variant->{start} . " and " . $variant->{stop} . 
                " on transcript " . $structure->transcript_transcript_name . 
                ", cannot continue annotation of coding exon");
        return;
    }

    # Need to create an amino acid change string for the protein domain method
    my ($protein_domain, $all_protein_domains) = $self->_protein_domain(
        $structure, $variant, $protein_position
    );

    return (
            c_position => "c." . $coding_position,
            trv_type => $trv_type,
            amino_acid_change => $protein_string,
            domain => $protein_domain,
            all_domains => $all_protein_domains,
           );
}


# Taken from Genome::Transcript
# Given a version and species, find the imported reference sequence build
sub get_reference_build_for_transcript {
    my($self, $structure) = @_;

    my ($version) = $structure->transcript_version =~ /^\d+_(\d+)[a-z]/;
    my $species = $structure->transcript_species;

    unless ($self->{'_reference_builds'}->{$version}->{$species}) {
        my $build = Genome::Model::Build->get($self->reference_sequence_id);
        confess "Could not get build version $version" unless $build;

        $self->{'_reference_build'}->{$version}->{$species} = $build;
    }
    return $self->{_reference_build}->{$version}->{$species};
}


#given a transcript, structure, and variant, calculate all of the bases after and including the indel 
#plus a number of bases before the indel to create correct, complete codons.  Then translate to an 
#amino acid sequence and return it. 
#If this is used near the centromere or at the end of a chromosome, it will be unpredictable.  Use this at those positions at your own risk!
sub _apply_indel_and_translate{
    my ($self, $structure, $variant) = @_;

    my $chrom_name = $structure->chrom_name;

    #
    # WARNING - GIANT HACK AHEAD!!!!!
    #
    # Above in _create_iterator_for_variant_intersection(), it adds a function ref to
    # $Genome::DataSource::TranscriptStructures::intersector_sub which is a hook into 
    # the data loader for dynamically changing filters while an iterator is in progress.
    #
    # For the get() immediately below, we actually want to get TranscriptStructures without
    # that extra filtering.  So, to temporarily turn off the hook, set that subref to undef
    #

    no strict 'refs';
    my $intersect_sub_name = $self->_resolve_intersector_sub_name();
    local($$intersect_sub_name);
    my $structures_class = $self->transcript_structure_class_name;
    my @sibling_structures = $structures_class->get(transcript_transcript_id => $structure->transcript_transcript_id,
                                                    chrom_name => $chrom_name,
                                                    data_directory => $structure->data_directory,
                                                    'structure_start <' => $structure->transcript_transcript_stop + 50000,
                                                    # structure_type => 'intron',
                                                  );
    if ($structure->transcript_strand eq '+1') {
        @sibling_structures = sort { $a->structure_start <=> $b->structure_start } @sibling_structures;
    } else {
        @sibling_structures = sort { $b->structure_start <=> $a->structure_start } @sibling_structures;
    }

                                                       
    my @structures = ($structure);

    for my $substructure (@sibling_structures) {
        if($structure->transcript_is_after($substructure->structure_start, $structure->structure_start)){
            if($substructure->structure_type ne 'intron'){
                push @structures, $substructure;
            }
        }
    }
    my $sequence = "";
    $sequence = $structure->phase_bases_before if $structure->phase_bases_before ne 'NULL';
    Genome::Model::Tools::Sequence->class();
    for my $substructure (@structures) {
        if ($substructure->structure_type eq 'flank') {
            my $flank_sequence = Genome::Model::Tools::Sequence::lookup_sequence(
                    chromosome => $chrom_name,
                    start => $substructure->structure_start,
                    stop => $substructure->structure_stop, 
                    build => $self->get_reference_build_for_transcript($substructure),
                    ); 
            die "Unsuccessfully executed sequence fetch" unless $flank_sequence;
            if ($structure->transcript_strand eq '-1'){
                $flank_sequence = $self->reverse_complement($flank_sequence);
            }
            $sequence .= $flank_sequence;
        }
        else {
            $sequence .= $substructure->nucleotide_seq;
        }
    }

    my ($start, $stop) = ($variant->{start}, $variant->{stop});
    ($start, $stop) = ($stop, $start) if $structure->transcript_strand eq '-1';
    my $length = $stop - $start + 1;
    
    my ($relative_variant_start, $borrowed_bases_start) = $structure->sequence_position($start);

    my ($codon, $codon_position) = $structure->codon_at_position($start); 

    my $first_affected_codon_start = $relative_variant_start - $codon_position;

    $sequence = substr($sequence, $first_affected_codon_start);

    my $mutated_seq;
    if ($variant->{type} eq 'DEL') {
        my $first = shift @structures;
        #my $deleted_bases = abs(min($first->stop_with_strand, $stop) - $start) + 1;
        my $deleted_bases = min(abs($first->stop_with_strand - $start), abs($stop - $start)) + 1;
        for my $substructure (@structures) {
            # last if $stop < $structure->structure_start;
            last if $structure->transcript_is_after($substructure->start_with_strand, $stop); 
            #if ($stop < $substructure->structure_stop) {
            if($structure->transcript_is_after($substructure->stop_with_strand, $stop)){
                $deleted_bases += abs($substructure->start_with_strand - $stop);
            }
            else {
                $deleted_bases += $substructure->length;
            }
        }
        $mutated_seq = substr($sequence, 0, $codon_position) . substr($sequence, $codon_position + $deleted_bases);
    }
    elsif ($variant->{type} eq 'INS') {
        $mutated_seq = substr($sequence, 0, $codon_position + 1) . $variant->{variant} . substr($sequence, $codon_position + 1);
        $mutated_seq = substr($mutated_seq, 3) if $codon_position == 2; #When codon_position == 2, there's an extra leading codon that needs to get hacked off here so we don't mislead the analysts into thinking an extra codon changed
    }
    else {
        confess "Malformed variant type: ";
    }
    
    my $aa_sequence = $self->translate($mutated_seq, $chrom_name); #if there's no stop, we don't care.
    $aa_sequence =~ s/\*.*/\*/; #remove everything after the amino_acid stop
    return $aa_sequence;
}

# Given a transcript, structure, and position, returns all of the bases in the
# coding region of the transcript 3' of the position
sub _coding_bases_after_position {
    my ($self, $transcript, $structure, $position) = @_;
    my ($sequence_position, $codon_position) = $structure->sequence_position($position);
    my $coding_bases_before = $structure->coding_bases_before;
    
    my $coding_sequence = Genome::TranscriptCodingSequence->get(
        transcript_id => $transcript->id,
        data_directory => $transcript->data_directory,
    );

    my $sequence = substr($coding_sequence->sequence, $coding_bases_before + $sequence_position - $codon_position - 1);
    return $sequence;
}       
    
# Find the UCSC conservation score for a piece of chromosome
sub _ucsc_conservation_score {
    require Genome::Model::Tools::Annotate::TranscriptVariants::Version1;
    Genome::Model::Tools::Annotate::TranscriptVariants::Version1::_ucsc_conservation_score(@_);
}

# Find the domains affected by the variant and all domains for the transcript
sub _protein_domain {
    my ($self, $structure, $variant, $protein_position) = @_;
    return 'NULL', 'NULL' unless defined $structure and defined $variant;

    my @all_domains = Genome::InterproResult->get(
        transcript_name => $structure->transcript_transcript_name,
        data_directory => $structure->data_directory,
        chrom_name => $variant->{chromosome_name},
        );
    return 'NULL', 'NULL' unless @all_domains;

    my @variant_domains;
    my @all_domain_names;
    for my $domain (@all_domains) {
        if ($protein_position >= $domain->{start} and $protein_position <= $domain->{stop}) {
            push @variant_domains, $domain->{name};
        }
        push @all_domain_names, $domain->{name};
    }

    return 'NULL', join(",", uniq @all_domain_names) unless @variant_domains;
    return join(",", uniq @variant_domains), join(",", uniq @all_domain_names);
}


# For full description of this positioning convention, see http://www.hgvs.org/mutnomen/recs-DNA.html
# and/or http://www.hgmd.cf.ac.uk/docs/mut_nom.html
#
# Here's a summary: Positions 5' ("before") of the coding region of the transcript start with '-', positions 3'
# ("after") of the coding region start with '*'. The A of the ATG start codon (for human) is nucleotide 1, exons 
# count up from there. Introns, if the position is closer to the previous exon, are <previous_exon_bp>+
# <position_in_intron>, if closer to the next exon it's <next_exon_bp>-<position_downstream_in_intron>.
# Those + and - signs are the string literal, not the arithmetic operation.
#
# So, a UTR exon 50 bp before the beginning of the first cds exon would be -50. A UTR exon 50bp after the
# last cds exon would be *50. An intron 5bp after exon 1 would be <length of exon 1>+5, an intron 5bp before 
# exon 2 would be <length of exon 1 + 1>-5. Basepair 50 of exon 3 would be the sum of (length exon 1, length exon 2, 50). 
# Coding regions over a range are separated by a _, so intron would be <previous_exon_bp>+<position_in_intron>_<position_in_intron>
#
# Transcripts without any coding exons (for example, an RNA transcript or a transcript with only utr exons) are
# treated differently. Flanking regions use the transcript start/stop instead of coding region start/stop and
# utr exons use variant start. This behavior is undocumented and I've coded it just so it matches old logic.
#
# TODO This can probably be refactored and simplified, there's a lot of redundant logic here
sub _determine_coding_position {
    my ($self, $variant, $structure) = @_;
    my $cds_region_start = $structure->transcript_coding_region_start;
    my $cds_region_stop = $structure->transcript_coding_region_stop;
    my $structure_type = $structure->structure_type;
    my $start = $variant->{start};
    my $stop = $variant->{stop};
    my $c_position = 'NULL';
    my $coding_bases_before = $structure->coding_bases_before;
    my $strand = $structure->transcript_strand;

    if ($structure_type eq 'rna') {
        $c_position = 'NULL';
    }
    elsif ($structure_type eq 'flank') {
        if ($structure->transcript_has_coding_region) {
            my $distance_to_coding = $structure->transcript_distance_to_coding_region($start);
            if ($structure->transcript_before_coding_region($start)) {
                $c_position = "-" . $distance_to_coding; 
            }
            elsif ($structure->transcript_after_coding_region($start)) {
                $c_position = "*" . $distance_to_coding;
            }
        }
        else {
            my $distance_to_transcript = $structure->transcript_distance_to_transcript($start);
            if ($structure->transcript_is_before($start, $structure->transcript_transcript_start)) {
                $c_position = "-" . $distance_to_transcript;
            }
            else {
                $c_position = "*" . $distance_to_transcript;
            }
        }
    }
    elsif ($structure_type eq 'utr_exon') {
        if ($structure->transcript_has_coding_region) {
            my $distance_to_coding = $structure->transcript_distance_to_coding_region($start);
            if ($structure->transcript_before_coding_region($start)) {
                $c_position = "-" . $distance_to_coding;
            }
            elsif ($structure->transcript_after_coding_region($start)) {
                $c_position = "*" . $distance_to_coding;
            }
        }
        else {
            # I have found no documentation for this particular case and am coding this to match the old output.
            # If transcript strand is +1, then utr exon is always 3' and prefixed with a '*'. Otherwise, the
            # utr exon is considered 5' and prefixed with a '-'. The position is simply the variant start position
            # TODO ... why is the logic like this? Is this what we want?
            if ($structure->transcript_strand eq '+1') {
                $c_position = "*" . $start;
            }
            else {
                $c_position = '-' . $start;
            }
        }
    }
    elsif ($structure_type eq 'cds_exon') {
        my ($start_seq_position) = $structure->sequence_position($start);
        my ($stop_seq_position) = $structure->sequence_position($stop);
        my $num_phase_bases = $structure->num_phase_bases_before;
        $c_position = $coding_bases_before + $start_seq_position - $num_phase_bases + 1;
        $c_position .= "_" . ($coding_bases_before + $stop_seq_position - $num_phase_bases + 1) if $start != $stop;
    }
    elsif ($structure_type eq 'intron') {
        ($start, $stop) = ($stop, $start) if $strand eq '-1';
        my $distance_to_start = $structure->distance_to_start($start) + 1;
        my $distance_to_stop = $structure->distance_to_stop($start) + 1;
        my $start_to_coding_distance = $structure->transcript_distance_to_coding_region($start);
        my $stop_to_coding_distance = $structure->transcript_distance_to_coding_region($stop);

        if ($structure->transcript_has_coding_region) {
            if ($structure->transcript_within_coding_region($start)) {
                if ($distance_to_start <= $distance_to_stop) {
                    $c_position = ($coding_bases_before) . "+" . $distance_to_start;
                    if ($start != $stop) {
                        my $dist = $distance_to_start - abs($stop - $start);
                        $dist = 1 if $dist <= 0;
                        $c_position .= "_" . ($coding_bases_before) . "+" . ($dist);
                    }
                }
                else {
                    $c_position = ($coding_bases_before + 1) . "-" . $distance_to_stop;
                    if ($start != $stop) {
                        my $dist = $distance_to_stop - abs($stop - $start);
                        $dist = 1 if $dist <= 0;
                        $c_position .= "_" . ($coding_bases_before + 1) . "-" . ($dist);
                    }
                }
            }
            else {
                if ($structure->transcript_before_coding_region($start)) {
                    $c_position = "1";
                }
                else {
                    $c_position = $structure->coding_bases_before;
                }

                if ($distance_to_start <= $distance_to_stop) {
                    $c_position .= ("+" . $distance_to_start);
                    $c_position .= "_" . ($distance_to_start + ($stop - $start)) if $start != $stop;
                }
                else {
                    $c_position .= ("-" . $distance_to_stop);
                    $c_position .= "_" . ($distance_to_stop - ($stop - $start)) if $start != $stop;
                }
            }
        }
    }
    else {
        $self->error_message("Unexpected substructure type $structure_type, " .
                "cannot calculate coding sequence position. Returning NULL!");
        $c_position = 'NULL';
    }

    return $c_position;
}

# Given a transcript, exon, and variant, determine the original and mutated sequence
# and the location of the first changed amino acid of the protein
sub _get_affected_sequence {
    my ($self, $structure, $variant) = @_;

    # Occasionally, the first amino acid changed by a variation occurs beyond the codons actually changed. For example,
    # CCC CTG CTG codes for PLL. If the last base of the first codon to the last base of the middle codon are deleted 
    # (CCTG removed, leaving CC CTG), then the first amino acid expressed is still P. The remaining bases TG aren't 
    # enough to see what the next amino acid would be so we still don't know what the first changed amino acid is! 
    # Looking a little further helps prevent this, but the below number is arbitrary. Extra codons are placed before
    # and after the variant so the extra sequence is available for both forward and reverse stranded transcript
    my $indel_extra_codons = 3;
    my $dnp_extra_codons = 1;

    # Determine protein position of amino acid at variant start. This will be modified below for indels to account
    # for extra codons placed at the beginning of the original sequence and further modified during reduction
    my ($variant_position, $phase_bases);
    if ($structure->transcript_strand eq '-1') {
        ($variant_position, $phase_bases) = $structure->sequence_position($variant->{stop});
    }
    else {
        ($variant_position, $phase_bases) = $structure->sequence_position($variant->{start});
    }
    my $coding_bases_before = $structure->coding_bases_before;
    my $protein_position = int(($coding_bases_before + $variant_position - $phase_bases) / 3) + 1;

    my ($orig_seq, $orig_codon_pos, $mutated_seq, $relative_start, $relative_stop);
    if ($variant->{type} eq 'SNP') {
        ($orig_seq, $orig_codon_pos) = $structure->codon_at_position($variant->{start});
        $mutated_seq = substr($orig_seq, 0, $orig_codon_pos) .
            $variant->{variant} . 
            substr($orig_seq, $orig_codon_pos + 1);
    }
    elsif ($variant->{type} eq 'DNP') {
        ($orig_seq, $relative_start, $relative_stop) = $structure->codons_in_range($variant->{start}, $variant->{stop}, $dnp_extra_codons);

        my $codon_position = ($relative_start - 1) % 3;
        my $bases_before = $relative_start - $codon_position;
        my $codons_before = int($bases_before / 3);
        $protein_position -= $codons_before;

        $mutated_seq = substr($orig_seq, 0, $relative_start - 1) .
            $variant->{variant} .
            substr($orig_seq, $relative_stop);
    }
    else {
        # Get codons between variant start and stop as well as extra codons to either side.
        ($orig_seq, $relative_start, $relative_stop) = $structure->codons_in_range($variant->{start}, $variant->{stop}, $indel_extra_codons);

        # Adjusting the protein position to point at the first amino acid in the original sequence requires more than just
        # using the extra_codons value. The biggest example of this is when the variation occurs at the beginning of the
        # exon. If there are fewer than 3 codons before the variation, then only what's available is concatenated. 
        # In addition, sequence from the end of the previous exon is also added (unless there are no previous exons!). 
        # So, the best way to determine what's been added is to see how much sequence exists before the variation.
        my $codon_position = ($relative_start - 1) % 3;
        my $bases_before = $relative_start - $codon_position;
        my $codons_before = int($bases_before / 3);
        $protein_position -= $codons_before;

        #it is possible that the variant goes off the end of the transcript.  In this case,
        #we need to adjust the relative stop.
        if ($relative_stop > length($orig_seq)) {
            $relative_stop = length($orig_seq);
        }

        if ($variant->{type} eq 'DEL') {
            $mutated_seq = substr($orig_seq, 0, $relative_start - 1) .
                           substr($orig_seq, $relative_stop);
        }
        elsif ($variant->{type} eq 'INS') {
            $protein_position-- if $codon_position == 2; # This little fix is due to the variant start of insertions being
                                                         # the base before the insertion instead of the first base of variation
                                                         # as is the case with deletions. 
            $mutated_seq = substr($orig_seq, 0, $relative_start) .
                           $variant->{variant} .
                           substr($orig_seq, $relative_start);
        }
    }

    return (
            $orig_seq,
            $mutated_seq,
            $protein_position
           );
}

# Takes two amino acid sequences, one an original and the other mutated, and removes
# any amino acids that haven't been changed as a result of the mutation
# Returns the reduced sequences and the offset that should be applied to protein position
#
# For example, if the original sequence is MRPRST and the mutated sequence is MRTRPT, then the
# reduced original would be PRS and the reduced mutated would be TRP. The leading MR on both is
# removed as is the trailing T, since these do not change as a result of the mutation.
#
# TODO There's got to be a more elegant way to do this
sub _reduce {
    my ($self, $orig, $mutated) = @_;
    return ($orig, $mutated, 0) unless defined $orig and defined $mutated;

    my ($front_orig, $front_mutated) = ("", ""); # Holds aa sequence after screening forward
        my $offset = 0;
    my $differ = 0;

# Compare sequences forward to back, updating protein position until first difference is found
# Since protein position should point to first changed amino acid, it shouldn't be updated
# after a difference is found. 
# The tricky bit is stop codons. If a stop codon is found in the mutated sequence, the loop 
# stops because all sequence after that point wouldn't be translated. However, original sequence
# from that point forward still needs to be kept.
# Also, if the lengths of mutated sequence and original sequence differ, the extra will need
# to be kept as well.
    my $i;
    my $last_mutated_aa = "";
    for ($i = 0; $i < length $orig and $i < length $mutated; $i++) {
        my $curr_orig = substr($orig, $i, 1);
        my $curr_mutated = substr($mutated, $i, 1);
        $last_mutated_aa = $curr_mutated;
        if ($curr_orig eq $curr_mutated and $differ == 0) {
            $offset++;
        }
        else {
            $differ = 1;
            $front_mutated .= $curr_mutated;
            last if $curr_mutated eq '*'; 
            $front_orig .= $curr_orig;
        }
    }
    if (length $orig > length $mutated or $last_mutated_aa eq '*') {
        $front_orig .= substr($orig, $i);
    }
    if (length $mutated > length $orig and $last_mutated_aa ne '*') {
        $front_mutated .= substr($mutated, $i);
    }

# Now need to compare sequences starting from the end, removing amino acids that are the same until
# a difference is found. In this case, no updating of the protein position is required.
    $differ = 0;
    my @orig_backward = ();
    my @mutated_backward= ();
    while ($front_orig and $front_mutated) {
        my $curr_orig = chop $front_orig;
        my $curr_mutated = chop $front_mutated;
        if ($curr_orig ne $curr_mutated or $differ) {
            push @orig_backward, $curr_orig;
            push @mutated_backward, $curr_mutated;
            $differ = 1;
        }
    }

    $front_orig = "" unless $front_orig;       # Prevents warnings about undefined strings
        $front_mutated = "" unless $front_mutated;
    my $new_orig = $front_orig . join("", reverse @orig_backward);
    my $new_mutated = $front_mutated . join("", reverse @mutated_backward);

    return ($new_orig, $new_mutated, $offset);
}

1;

=pod
=head1 Name

Genome::Transcript::VariantAnnotator

=head1 Synopsis

Given a variant, all transcripts affected by that variant are annotated and returned

=head1 Usage

# Variant file tab delimited, columns are chromosome, start, stop, reference, variant
# Need to infer variant type (SNP, DNP, INS, DEL) as well
my $variant_file = variants.tsv;
my @headers = qw/ chromosome_name start stop reference variant /;
my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $variant_file,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
        );

my $model = Genome::Model->get(name => 'NCBI-human.combined-annotation');
my $build = $model->build_by_version('54_36p');
my $iterator = $build->transcript_iterator;
my $window = Genome::Utility::Window::Transcript->create(
        iterator => $iterator,
        range => 50000,
        );
my $annotator = Genome::Transcript::VariantAnnotator->create(
        transcript_window => $window
        );

while (my $variant = $reader->next) {
    my @annotations = annotator->transcripts($variant);
}

=head1 Methods

=head2 transcripts

=over

=item I<Synopsis>   gets all annotations for a variant

=item I<Arguments>  variant (hash; see 'Variant Properites' below)

=item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

=back

=head2 prioritized_transcripts

=over

=item I<Synopsis>   Gets one prioritized annotation per gene for a variant(snp or indel)

    =item I<Arguments>  variant (hash; see 'Variant properties' below)

    =item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

    =back

    =head2 prioritized_transcript

    =over

    =item I<Snynopsis>  Gets the highest priority transcript affected by variant

    =item I<Arguments>  variant (hash, see 'Variant properties' below)

    =item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

    =back

    =head1 Variant Properties

    =over

    =item I<chromosome_name>  The chromosome of the variant

    =item I<start>            The start position of the variant

    =item I<stop>             The stop position of the variant

    =item I<variant>          The snp base

    =item I<reference>        The reference base at the position

    =item I<type>             snp, dnp, ins, or del

    =back

    =head1 Annotation Properties

    =over

    =item I<transcript_name>    Name of the transcript

    =item I<transcript_source>  Source of the transcript

    =item I<strand>             Strand of the transcript

    =item I<c_position>         Relative position of the variant

    =item I<trv_type>           Called Classification of variant

=item I<priority>           Priority of the trv_type (only from get_prioritized_annotations)

    =item I<gene_name>          Gene name of the transcript

    =item I<intensity>          Gene intenstiy

    =item I<detection>          Gene detection

    =item I<amino_acid_length>  Amino acid length of the protein

    =item I<amino_acid_change>  Resultant change in amino acid in snp is in cds_exon

    =item I<variations>         Hashref w/ keys of known variations at the variant position

    =item I<type>               snp, ins, or del

    =back

    =head1 See Also

    B<Genome::Model::Command::Report>

    =head1 Disclaimer

    Copyright (C) 2008 Washington University Genome Sequencing Center

    This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

    Core Logic:

    B<Xiaoqi Shi> I<xshi@genome.wustl.edu>

    Optimization:

    B<Eddie Belter> I<ebelter@watson.wustl.edu>

    B<Gabe Sanderson> l<gsanders@genome.wustl.edu>

    B<Adam Dukes l<adukes@genome.wustl.edu>

    B<Brian Derickson l<bdericks@genome.wustl.edu>

    =cut

#$HeadURL$
#$Id$
