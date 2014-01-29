package Genome::Transcript;
#:adukes short term: move data directory into id_by, but this has to be done in parallel w/ rewriting all file-based data sources.  It might be better to wait until long term: custom datasource that incorporates data_dir, possibly species/source/version, eliminating the need for these properties in the id, and repeated multiple times in the files

use strict;
use warnings;

use Genome;
use Carp;
use Bio::Tools::CodonTable;
use Bio::Seq;

class Genome::Transcript {
    table_name => 'TRANSCRIPT',
    id_by => [
        chrom_name => {
            is => 'Text',
        },
        transcript_start => {
            is => 'NUMBER',
            is_optional => 1
        },
        transcript_stop => {
            is => 'NUMBER',
            is_optional => 1,
        },
        species => { is => 'varchar',
            is_optional => 1,
        },
        source => { is => 'VARCHAR',
            is_optional => 1,
        },
        version => { is => 'VARCHAR',
            is_optional => 1,
        },
        transcript_id => {
            is => 'NUMBER',
        },
    ],
    has => [
        reference_build_id => {
            is => 'Text',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
        },
        gene_id => {
            is => 'Text',
        },
        gene_name => {
            is => 'Text',
            is_optional => 1,
        },
        transcript_name => {
            is => 'VARCHAR',
            is_optional => 1,
        },
        transcript_status => { is => 'VARCHAR',
            is_optional => 1,
            valid_values => ['reviewed', 'unknown', 'model', 'validated', 'predicted', 'inferred', 'provisional', 'unknown', 'known', 'novel', 'putative'],
        },
        strand => { is => 'VARCHAR',
            is_optional => 1,
            valid_values => ['+1', '-1', 'UNDEF'],
        },
        sub_structures => {
            #is_constant => 1,
            calculate_from => [qw/ id  data_directory chrom_name/],
            calculate => q|
            Genome::TranscriptStructure->get(chrom_name => $chrom_name, transcript_id => $id, data_directory => $data_directory);
            |,
        },
        protein => {
            calculate_from => [qw/ id data_directory reference_build_id/],
            calculate => q|
            Genome::Protein->get(transcript_id => $id, data_directory => $data_directory, reference_build_id => $reference_build_id);
            |,
        },
        gene => {
            calculate_from => [qw/ gene_id data_directory reference_build_id/],
            calculate => q|
            Genome::Gene->get(id => $gene_id, data_directory => $data_directory, reference_build_id => $reference_build_id);
            |,
        },
        data_directory => {
            is => "Path",
        },
        transcript_error => {
            is => 'Text',
            is_optional => 1,
            doc => 'Describes any error with the transcript, which affects its priority when compared with other transcripts',
        },
        coding_region_start => {
            is => 'Number',
            is_optional => 1,
            doc => 'The first nucleotide of the transcript\'s coding region, always less than stop (in other words, ignores strand)',
        },
        coding_region_stop => {
            is => 'Number',
            is_optional => 1,
            doc => 'The last nucleotide of the transcript\'s coding region, always greater than stop (ignores strand)',
        },
        amino_acid_length => {
            is => 'Number',
            is_optional => 1,
            doc => 'The length of the amino acid to which the coding region of the transcript translates',
        },
    ],
    schema_name => 'files',
    data_source => 'Genome::DataSource::Transcripts',
};

# Returns 1 if the position is before the coding region of the transcript
sub before_coding_region {
    my ($self, $position) = @_;
    return 0 if $self->coding_region_start eq 'NULL' or $self->coding_region_stop eq 'NULL';

    if ($self->strand eq '+1') {
        return 1 if $position < $self->coding_region_start;
        return 0;
    }
    elsif ($self->strand eq '-1') {
        return 1 if $position > $self->coding_region_stop;
        return 0;
    }
    else {
        $self->error_message("Transcript " . $self->transcript_name . " has an invalid strand: " . $self->strand);
        confess "Invalid strand on transcript " . $self->transcript_name;
    }
}

# Returns 1 if the position is after the coding region of the transcript
sub after_coding_region {
    my ($self, $position) = @_;
    return 0 if $self->coding_region_start eq 'NULL' or $self->coding_region_stop eq 'NULL';

    if ($self->strand eq '+1') {
        return 1 if $position > $self->coding_region_stop;
        return 0;
    }
    elsif ($self->strand eq '-1') {
        return 1 if $position < $self->coding_region_start;
        return 0;
    }
    else {
        $self->error_message("Transcript " . $self->transcript_name . " has an invalid strand: " . $self->strand);
        confess "Invalid strand on transcript " . $self->transcript_name;
    }
}

# Returns 1 if the position lies within the coding region of the transcript
sub within_coding_region {
    my ($self, $position) = @_;
    return 0 if $self->coding_region_start eq 'NULL' or $self->coding_region_stop eq 'NULL';
    return 1 if $position >= $self->coding_region_start and $position <= $self->coding_region_stop;
    return 0;
}

# Returns 1 if the provided position is within the start/stop range of the transcript, EXCLUDES flanking regions
sub within_transcript {
    my ($self, $position) = @_;
    if ($position >= $self->{transcript_start} and $position <= $self->{transcript_stop}) {
        return 1;
    }
    return 0;
}

# Returns 1 if the provided position is within the transcript with flanking structures included
sub within_transcript_with_flanks {
    my ($self, $position, $flanking_distance) = @_;
    if ($position >= ($self->{transcript_start} - $flanking_distance) and $position <= ($self->{transcript_stop} + $flanking_distance)) {
        return 1;
    }
    return 0;
}

# Returns the distance between the position and the transcript's start or stop
sub distance_to_transcript {
    my ($self, $position) = @_;
    my $start = $self->{transcript_start};
    my $stop = $self->{transcript_stop};
    if ($position < $start) {
        return $start - $position;
    }
    elsif ($position > $stop) {
        return $position - $stop;
    }
    else {
        return 0;
    }
}

sub distance_to_coding_region {
    my ($self, $position) = @_;
    return 0 if (! $self->has_coding_region or (index($self->transcript_error,'rna_with_coding_region') >= 0));
    my $coding_start = $self->coding_region_start;
    my $coding_stop = $self->coding_region_stop;
    return 0 if $coding_start eq 'NULL' or $coding_stop eq 'NULL';
    return 0 if $position >= $coding_start and $position <= $coding_stop;

    my $distance;
    if ($self->before_coding_region($position)) {
        $distance = abs($coding_start - $position) if $self->strand eq '+1';
        $distance = abs($coding_stop - $position) if $self->strand eq '-1';
    }
    else {
        $distance = abs($coding_stop - $position) if $self->strand eq '+1';
        $distance = abs($coding_start - $position) if $self->strand eq '-1';
    }

    return $distance;
}

# Returns 1 if first position is 5' of the second position
sub is_before {
    my ($self, $p1, $p2) = @_;
    my $flanking_distance = 50000;
    return 0 unless $self->within_transcript_with_flanks($p1, $flanking_distance) and $self->within_transcript_with_flanks($p2, $flanking_distance);
    if ($self->strand eq '+1') {
        return 1 if $p1 < $p2;
        return 0;
    }
    elsif ($self->strand eq '-1') {
        return 1 if $p1 > $p2;
        return 0;
    }
    else {
        confess "Invalid strand " . $self->strand . " on transcript " . $self->transcript_name;
    }
}

# Returns 1 if the first position is 3' of the second position
sub is_after {
    my ($self, $p1, $p2) = @_;
    my $flanking_distance = 50000;
    return 0 unless $self->within_transcript_with_flanks($p1, $flanking_distance) and $self->within_transcript_with_flanks($p2, $flanking_distance);
    if ($self->strand eq '+1') {
        return 1 if $p1 > $p2;
        return 0;
    }
    elsif ($self->strand eq '-1') {
        return 1 if $p1 < $p2;
        return 0;
    }
    else {
        confess "Invalid strand " . $self->strand . " on transcript " . $self->transcript_name;
    }
}

sub structure_at_position {
    my ($self, $position) = @_;
    my $flanking_distance = 50000;
    return unless $self->within_transcript_with_flanks($position, $flanking_distance);

    my @structures = $self->ordered_sub_structures;
    unless (@structures) {
        $self->debug_message("No sub-structures for transcript " . $self->transcript_name . ", cannot find structure at position $position");
        return;
    }

    for my $struct ( @structures ) {
        return $struct if $position >= $struct->structure_start
            and $position <= $struct->structure_stop;
    }

    $self->debug_message("No substructure found for transcript " . $self->transcript_name . " at position $position");
    return;
}

# Returns the transcript substructures (if any) that are between the given start and stop position
sub structures_in_range {
    my ($self, $start, $stop) = @_;
    ($start,$stop) = sort {$a <=> $b} ($start, $stop);
    my @structures = sort {$a->structure_start <=> $b->structure_start} $self->ordered_sub_structures;
    unless (@structures){
        $self->warning_message("No substructures for transcript " . $self->transcript_name);
        return;
    }

    my @structures_in_range;
    for my $structure (@structures){
        my $ss_start = $structure->structure_start;
        my $ss_stop = $structure->structure_stop;

        if (($ss_start >= $start and $ss_start <= $stop ) or
            ( $ss_stop >= $start and $ss_stop <= $stop ) or
            ( $ss_start <= $start and $ss_stop >=$stop ))
        {
            push @structures_in_range, $structure;
        }
    }
    return @structures_in_range;
}

# Gets the structure at the specified position and returns the structure(s) on either side of it
sub structures_flanking_structure_at_position {
    my ($self, $position) = @_;
    my $structure = $self->structure_at_position($position);
    return unless $structure;

    my @structures = $self->ordered_sub_structures;
    my $previous_structure;
    while (@structures) {
        my $s = shift @structures;
        if ($s->id eq $structure->id) {
            return ($previous_structure, shift @structures);
        }
        $previous_structure = $s;
    }

    return;
}

# Orders the sub structures using strand
sub ordered_sub_structures {
    my $self = shift;
    unless (exists $self->{'_ordered_sub_structures'}) {
        my @subs;
        if ($self->strand eq '+1') {
            @subs = sort { $a->structure_start <=> $b->structure_start } $self->sub_structures;
        }
        elsif ($self->strand eq '-1') {
            @subs = sort { $b->structure_start <=> $a->structure_start } $self->sub_structures;
        }
        $self->{'_ordered_sub_structures'} = \@subs;
    }
    return @{$self->{'_ordered_sub_structures'}};
}

# Translates basepair sequence into amino acid sequence
sub translate_to_aa {
    my $self = shift;
    my $seq = shift;
    my $length = length $seq;
    my $translator = $self->get_codon_translator;
    my $translation;
    for (my $i=0; $i<=$length-2; $i+=3) {
        my $codon=substr($seq, $i, 3);
        $codon =~ s/N/X/g;
        my $aa = $translator->translate($codon);
        $aa="*" if ($aa eq 'X');
        $translation.=$aa;
    }
    return $translation;
}

# Returns codon table object
sub get_codon_translator {
    my $self = shift;
    my $translator;
    if ($self->chrom_name =~ /^MT/){
        $translator = $self->{'_mitochondrial_codon_translator'};
        unless ($translator) {
            $translator = Bio::Tools::CodonTable->new(-id => 2);
            $self->{'_mitochondrial_codon_translator'} = $translator;
        }
    }
    else {
        $translator = $self->{'_codon_translator'};
        unless ($translator) {
            $translator = Bio::Tools::CodonTable->new(-id => 1);
            $self->{'_codon_translator'} = $translator;
        }
    }

    unless ($translator) {
        $self->error_message("Could not get codon table object");
        die;
    }

    return $translator;
}

# Performs a few checks on the transcript to determine if its valid
# All errors are joined by a : and set to the transcript_error field
# TODO Look into just returning the string here instead of setting the
# transcript error field, which can result in loss of information (for
# example, pseudogene is set during import using information that is only
# available during import, so the internal_stop_codon call wouldn't catch
# it. This change needs to be coupled with changing how the importer
# calls this method, which should be trivial.
sub is_valid {
    my $self = shift;

    my @errors;
    unless ($self->check_start_codon) {
        push @errors, 'no_start_codon';
    }
    unless ($self->substructures_are_contiguous) {
        push @errors, 'gap_between_substructures';
    }
    unless ($self->internal_stop_codon) {
        push @errors, 'pseudogene';
    }
    unless ($self->cds_region_has_stop_codon) {
        push @errors, 'no_stop_codon';
    }
    unless ($self->has_coding_region or $self->is_rna) {
        push @errors, 'no_coding_region';
    }
    unless ($self->check_intron_size) {
        push @errors, 'overly_large_intron';
    }
    unless ($self->correct_bp_length_for_exons) {
        push @errors, 'bad_bp_length_for_coding_region';
    }
    unless ($self->exon_seq_matches_genome_seq) {
        push @errors, 'mismatch_between_exon_seq_and_reference';
    }
    unless ($self->rna_with_no_coding_region) {
        push @errors, 'rna_with_coding_region';
    }

    if (@errors) {
        my $error_string = join(":", @errors);
        $self->transcript_error($error_string);
        return 0;
    }

    $self->transcript_error('no_errors');
    return 1;
}

# If the transcript is RNA and has a coding region, it's invalid
sub rna_with_no_coding_region {
    my $self = shift;
    return 1 unless $self->is_rna;

    if (defined $self->cds_full_nucleotide_sequence) {
        $self->warning_message("Transcript " . $self->transcript_name .
            " is RNA with coding regions!");
        return 0;
    }
    return 1;
}

# Checks that the transcript has a coding region
sub has_coding_region {
    my $self = shift;

    return 0 if $self->coding_region_start eq 'NULL' or $self->coding_region_stop eq 'NULL';
    return 1;

    ################ not needed.... ###################
    my $seq = $self->cds_full_nucleotide_sequence;
    return 1 if defined $seq and length $seq > 0;
    return 0;

}

# Checks if the transcript represents rna
sub is_rna {
    my $self = shift;
    my @subs = $self->sub_structures;
    for my $sub (@subs) {
        return 1 if $sub->structure_type eq 'rna';
    }
    return 0;
}

# Checks that the transcript has an associated gene
sub has_associated_gene {
    my $self = shift;
    my $gene = $self->gene;
    return 0 unless defined $gene;
    return 1;
}

# Checks that the transcript has either +1 or -1 for strand
sub has_valid_strand {
    my $self = shift;
    my $strand = $self->strand;
    return 1 if $strand eq '+1' or $strand eq '-1';
    return 0;
}

# Compares the protein aa sequence with the transcript's translated aa sequence
sub transcript_translation_matches_aa_seq {
    my $self = shift;
    return 1 unless $self->has_coding_region;
    return 1 if $self->is_rna;
    my $protein_aa_seq = $self->protein->amino_acid_seq;
    my $transcript_translation = $self->translate_to_aa($self->cds_full_nucleotide_sequence);
    unless ($protein_aa_seq eq $transcript_translation) {
        return 0;
    }
    return 1;
}

# Make sure that coding region starts with a start codon
sub check_start_codon {
    my $self = shift;

    return 1 unless $self->has_coding_region;
    return 1 if $self->is_rna;                 # Rna doesn't have a coding region, and so doesn't have a start codon
    return 1 unless $self->species eq 'human'; # Start codon check only works for human

    my $seq = $self->cds_full_nucleotide_sequence;
    return 1 unless defined $seq;

    my $translator = $self->get_codon_translator;
    return $translator->is_start_codon(substr($seq, 0, 3))
}

# Ensures that no introns are larger than 900kb
sub check_intron_size {
    my $self = shift;
    my @introns = $self->introns;
    foreach my $intron (@introns) {
        return 0 if $intron->length > 900000
    }
    return 1;
}

# Checks that each exon's nucleotide sequence matches the reference sequence
sub exon_seq_matches_genome_seq {
    my $self = shift;
    my @exons = $self->cds_exons;
    my $ref_build = $self->get_reference_build;

    foreach my $exon (@exons) {
        my $ref_seq = $ref_build->sequence(
            $self->chrom_name, $exon->structure_start, $exon->structure_stop
        );

        unless ($ref_seq) {
            # Some of the imported reference sequences don't have mitochondrial data, which isn't the
            # transcript's fault, so don't mark it as having an error
            if ($self->chrom_name =~ /^[MN]T*/i) {
                return 1;
            }
            else {
                $self->warning_message("Could not get sequence from chromosome " . $self->chrom_name .
                    " between positions " . $exon->structure_start . " and " . $exon->structure_stop .
                    " for transcript " . $self->transcript_name . " and structure " . $exon->transcript_structure_id);
                return 0;
            }
        }

        if ($self->strand eq '-1') {
            $ref_seq = $self->reverse_complement($ref_seq);
        }

        my $exon_seq = $exon->nucleotide_seq;
        return 0 unless $ref_seq eq $exon_seq;
    }
    return 1;
}

# Ensures that the coding region has a stop codon at the end
sub cds_region_has_stop_codon {
    my $self = shift;

    return 1 if $self->is_rna;
    return 1 unless $self->has_coding_region;

    my $seq = $self->cds_full_nucleotide_sequence;
    return 1 unless defined $seq and length $seq > 0;

    my $aa = $self->translate_to_aa($seq);
    if (substr($aa, -1) eq "*") {
        return 1;
    }
    else {
        return 0;
    }
}

# Checks for stop codons in the middle of the coding region
sub internal_stop_codon {
    my $self = shift;

    return 1 if $self->is_rna;
    return 1 unless $self->has_coding_region;

    my $seq = $self->cds_full_nucleotide_sequence;
    return 1 unless defined $seq and length $seq > 0;

    my $aa = $self->translate_to_aa($seq);
    my $stop = index($aa, "*");
    unless ($stop == -1 or $stop == length($aa) - 1) {
        return 0;
    }
    return 1;
}

# Checks that the coding region base pairs are correctly grouped into 3bp codons
sub correct_bp_length_for_exons {
    my $self = shift;
    return 1 unless $self->has_coding_region;
    return 1 if $self->is_rna;

    my $seq = $self->cds_full_nucleotide_sequence;
    return 1 unless defined $seq and length $seq > 1;
    if ($self->chrom_name =~ /^MT/) {
        return 1 if (length $seq) % 3 == 2;
    }
    return 1 if (length $seq) % 3 == 0;
    return 0;
}

# Checks that all exons (coding and untranslated) and introns are contiguous
sub substructures_are_contiguous {
    my $self = shift;
    my @ss = $self->ordered_sub_structures;
    my $prev;
    while (my $s = shift @ss){
        if ($prev) {
            if ($self->strand eq '+1') {
                return 0 unless $s->{structure_start} == $prev->{structure_stop} + 1
            }
            else {
                return 0 unless $s->{structure_stop} == $prev->{structure_start} - 1;
            }
        }
        $prev = $s;
    }

    return 1;
}

# Takes a sequence and returns its reverse complement
sub reverse_complement {
    my $self = shift;
    my $seq = shift;
    return unless defined $seq;

    my $s = Bio::Seq->new(-display_id => "junk", -seq => $seq);
    my $rev_com = $s->revcom->seq;
    unless ($rev_com) {
        $self->error_message("Could not create reverse complement for sequence");
        confess;
    }

    return $rev_com;
}

# Given a version and species, find the imported reference sequence build
sub get_reference_build {
    my $self = shift;
    return $self->reference_build;
}

# Returns all coding exons associated with this transcript
sub cds_exons {
    my $self = shift;

    my @ex = grep { $_->structure_type eq 'cds_exon' } $self->ordered_sub_structures;
    return @ex;
}

# Returns all introns associated with this transcript
sub introns {
    my $self = shift;

    my @int = grep { $_->structure_type eq 'intron' } $self->ordered_sub_structures;
    return @int;
}

# Returns all untranslated exons associated with this transcript
sub utr_exons{
    my $self = shift;

    my @utr_ex = grep { $_->structure_type eq 'utr_exon' } $self->ordered_sub_structures;
    return @utr_ex;
}

# Returns all flanking regions associated with this transcript
sub flanks{
    my $self = shift;

    my @flanks = grep { $_->structure_type eq 'flank'} $self->ordered_sub_structures;
    return @flanks
}

#returns the utr exon sequence 3' of the coding sequence
sub utr_exon_sequence_after_coding_sequence{
    my $self = shift;

    my @utr_exons_after_cds = grep { $_->coding_bases_after == 0 } $self->utr_exons;
    my $utr_sequence_after_cds = "";
    for my $utr (@utr_exons_after_cds){
        $utr_sequence_after_cds .=  $utr->nucleotide_seq;
    }
    return $utr_sequence_after_cds; #this coming out in the proper order depends on utr_exons calling ordered_sub_structures
}

#returns the flank region sequence 3' of the coding sequence
sub flank_after_coding_sequence{
    my $self =shift;

    my @flanks = $self->flanks;
    my @flanks_after_cds;
    for my $flank (@flanks){
        if($self->strand eq '+1'){
            push(@flanks_after_cds, $flank) if $flank->structure_start >= $self->transcript_stop;
        }else{
            push(@flanks_after_cds, $flank) if $flank->structure_stop <= $self->transcript_start;
        }
    }
    return @flanks_after_cds; #this coming out in the proper order depends on utr_exons calling ordered_sub_structures
}

# Returns the start position of the first exon and the stop position the last exon on the transcript
sub cds_exon_range {
    my $self = shift;

    my @cds_exons = $self->cds_exons or return;
    @cds_exons = sort { $a->{structure_start} <=> $b->{structure_stop} } @cds_exons;

    return ($cds_exons[0]->structure_start, $cds_exons[$#cds_exons]->structure_stop);
}

# Returns the total length of all exons before the position (does not include structure at position)
sub length_of_cds_exons_before_structure_at_position {
    my ($self, $position) = @_;
    return $self->_length_of_cds_exons_before_or_after_position($position, 1);
}

# Returns the total length of all exons after the position (does not include structure at position)
sub length_of_cds_exons_after_structure_at_position {
    my ($self, $position) = @_;
    return $self->_length_of_cds_exons_before_or_after_position($position, 0);
}

sub _length_of_cds_exons_before_or_after_position {
    my ($self, $position, $look_before) = @_;
    return unless defined $position and defined $look_before;
    my @cds_exons = $self->cds_exons or return 0;
    my $structure = $self->structure_at_position($position);
    unless ($structure) {
        $self->warning_message("No structure at position $position for length_of_cds_exons_after_structure_at_position!");
        return;
    }

    my $length = 0;
    for my $e (@cds_exons) {
        if ($look_before and $self->is_before($e->{structure_start}, $structure->{structure_start})) {
            $length += $e->length;
        }
        elsif (not $look_before and $self->is_after($e->{structure_start}, $structure->{structure_start})) {
            $length += $e->length;
        }
    }

    return $length;
}

sub structures_after_position {
    my ($self, $position) = @_;
    my @ss = $self->ordered_sub_structures;

    while (@ss) {
        if ($self->strand eq '+1' and $position < $ss[0]->structure_start) {
            return @ss;
        }
        elsif ($self->strand eq '-1' and $position > $ss[0]->structure_start) {
            return @ss;
        }
        shift @ss;
    }

    return;
}

# Grab only those coding exons that have ordinal defined
sub cds_exon_with_ordinal {
    my ($self, $ordinal) = @_;

    foreach my $cds_exon ( $self->cds_exons ) {
        return $cds_exon if $cds_exon->ordinal == $ordinal;
    }

    return;
}

# Full base pair sequence of coding regions
sub cds_full_nucleotide_sequence{
    my $self = shift;
    my $seq;
    foreach my $cds_exon ( sort { $a->ordinal <=> $b->ordinal} $self->cds_exons ) {
        $seq .= $cds_exon->nucleotide_seq;
    }
    return $seq;
}


# Returns name of associated gene
#sub gene_name
#{
#    my $self = shift;
#
#    my $gene = $self->gene;
#    my $gene_name = $gene->name($self->source);;
#
#    return $gene_name;
#}

# Returns strand as either + or - (used for bed string below)
sub strand_string {
    my $self = shift;
    my $strand = '.';
    if ($self->strand eq '+1') {
        $strand = '+';
    } elsif ($self->strand eq '-1') {
        $strand = '-';
    }
    return $strand;
}

# Returns string containing transcript info in bed format
sub bed_string {
    my $self = shift;
    # BED string should only be written if BED12 format is adopted, otherwise it's just one line per sub structure
    return;
    # BED format uses zero-based start positions
    my $bed_string = $self->chrom_name ."\t". ($self->transcript_start - 1) ."\t". $self->transcript_stop ."\t". $self->transcript_name ."\t0\t". $self->strand_string;
    return $bed_string ."\n";
}

# Base string for gff format
sub _base_gff_string {
    my $self = shift;
    return $self->chrom_name ."\t". $self->source .'_'. $self->version ."\t". 'transcript' ."\t". $self->transcript_start ."\t". $self->transcript_stop ."\t.\t". $self->strand_string ."\t.";
}

# Returns string containing transcript info in gff file format
sub gff_string {
    my $self = shift;
    return $self->_base_gff_string ."\t". $self->gene->name ."\n";
}

# Returns string containing transcript info in gff3 file format
sub gff3_string {
    my $self = shift;
    return $self->_base_gff_string ."\tID=".$self->transcript_id ."; NAME=". $self->transcript_name ."; PARENT=". $self->gene->gene_id .';' ."\n";
}

sub gtf_string {
    my $self = shift;
    my @sub_structure = grep {$_->structure_type ne 'flank' && $_->structure_type ne 'intron'} $self->ordered_sub_structures;
    my %exon_sub_structures;
    my $i = 1;
    for my $ss (@sub_structure){
        if ($ss->ordinal) {
            push @{$exon_sub_structures{$ss->ordinal}}, $ss;
        } else {
            push @{$exon_sub_structures{$i++}}, $ss;
        }
    }
    my $string;
    for my $ordinal ( sort {$a <=> $b} keys %exon_sub_structures ) {
        my @cds_strings;
        my $exon_start;
        my $exon_stop;
        my $exon_strand;
        for my $ss (@{$exon_sub_structures{$ordinal}}) {
            my $type = $ss->structure_type;
            if ($type =~ /intron/ || $type =~ /flank/) {
                next;
            } elsif ($type eq 'cds_exon') {
                $type = 'CDS';
                # These are just ignored by cufflinks and are duplicates, only useful if exon was comprised of utr_exon and cds_exon
                my $gene_name = $self->gene->name("ensembl_default_external_name");
                my $gene_id = $self->gene->name("ensembl");
                unless ($gene_id) { $gene_id = $gene_name }
                push @cds_strings, $ss->chrom_name ."\t". $ss->source .'_'. $ss->version ."\t". $type ."\t". $ss->structure_start ."\t". $ss->structure_stop ."\t.\t". $ss->strand ."\t". $ss->frame ."\t" .' gene_name "'. $gene_name .'"; gene_id "'. $gene_id .'"; transcript_id "'. $self->transcript_name .'"; exon_number "'. $ordinal .'";';
            } elsif ($type eq 'rna') {
                # Should RNA even be included as if it's coding, just as an exon or included at all...
                $type = 'RNA';
                my $frame = $ss->frame;
                if (!defined($frame) || $frame =~ /^\s*$/) {
                    # No frame attribute for rna
                    $frame = '.';
                }
                #It seems cufflinks now fails if multiple annotations (ie. RNA and exon) exist in the GTF file for a single transcript_id
                #push @cds_strings, $ss->chrom_name ."\t". $ss->source .'_'. $ss->version ."\t". $type ."\t". $ss->structure_start ."\t". $ss->structure_stop ."\t.\t". $ss->strand ."\t". $frame ."\t".' gene_id "'. $ss->gene_name .'"; transcript_id "'. $ss->transcript_name .'"; exon_number "'. $ordinal .'";';
            }
            unless ($exon_start && $exon_stop) {
                $exon_start = $ss->structure_start;
                $exon_stop = $ss->structure_stop;
                $exon_strand = $ss->strand;
            } else {
                if ($ss->structure_start < $exon_start) {
                    # Should never happen since ss are ordered
                    $exon_start = $ss->structure_start;
                }
                if ($ss->structure_stop > $exon_stop) {
                    $exon_stop = $ss->structure_stop;
                }
                if ($ss->strand ne $exon_strand) {
                    #This should never happen
                    die('Inconsistent strand on transcript '. $self->transcript_name);
                }
            }
        }
        unless ($exon_start && $exon_stop) {
            die(Data::Dumper::Dumper($self));
        }
        my $gene_name = $self->gene->name('ensembl_default_external_name');
        my $gene_id = $self->gene->name('ensembl');
        unless ($gene_id) { $gene_id = $gene_name }
        $string .= $self->chrom_name ."\t". $self->source .'_'. $self->version ."\texon\t". $exon_start ."\t". $exon_stop ."\t.\t". $exon_strand ."\t.\t".' gene_name "'. $gene_name .'"; gene_id "'. $gene_id .'"; transcript_id "'. $self->transcript_name .'"; exon_number "'. $ordinal .'";' ."\n";
        if (scalar(@cds_strings)) {
            $string .= join("\n", @cds_strings) ."\n";
        }
    }
    return $string;
}

1;

#TODO
=pod


=cut

