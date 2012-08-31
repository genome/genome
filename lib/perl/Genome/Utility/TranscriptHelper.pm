package Genome::Utility::TranscriptHelper;

use strict;
use warnings;
use Genome;

class Genome::Utility::TranscriptHelper{
    is => 'UR::Object',
    has => [
        model_name => {
            is => 'Text',
            is_optional => 1,
            default => 'NCBI-human.combined-annotation',
        },
        version => {
            is => 'Text',
            is_optional => 1,
            default => '54_36p_v2',
        }
    ],
};

sub _get_transcript_window {
    # Given a chromosome name, model, and version, returns a transcript window
    my $self = shift;
    my $chromosome = shift;

    my $build = Genome::Model->get(name=>$self->model_name)->build_by_version($self->version);
    my $transcript_iterator = $build->transcript_iterator(chrom_name => $chromosome);

    unless (defined $transcript_iterator) {
        $self->error_message("Could not create transcript iterator for chromosome $chromosome!");
        return;
    }

    my $transcript_window = Genome::Utility::Window::Transcript->create(iterator => $transcript_iterator);

    unless (defined $transcript_window) {
        $self->error_message("Could not create transcript window for chromosome $chromosome!");
        return;
    }

    return $transcript_window;
}

sub annotate_sv {
    # Given a SV type, two chromosomes, and two breakpoints, returns a hash containing information about
    # each breakpoint and, if SV type is DEL, a list of all genes that have been deleted between
    # breakpoint a and breakpoint b
    # Expecting self, type, chrom_a, break_a, chrom_b, break_b
    my $self = shift;
    my ($type, $chroma, $posa, $chromb, $posb) = @_;

    my %anno;
    $anno{type} = $type;
    $anno{breakpoint_a} = {chromosome => $chroma, position => $posa};
    $anno{breakpoint_b} = {chromosome => $chromb, position => $posb};

    if ($posa > $posb and $chroma eq $chromb) {
        $self->error_message("Breakpoint A should be less than breakpoint B on same chromosome");
        return;
    }
    
    unless (grep {$type eq $_} ("DEL" , "INS" , "INV" , "CTX") ) {
        $self->error_message("$type is an invalid event type, valid values are DEL, INS, INV, CTX");

        return;
    }

    if ( ($chroma ne $chromb) and ($type ne "CTX") ) {
        $self->error_message("Only CTX can occur on two different chromosomes");
        return;
    }

    my $transcript_window = $self->_get_transcript_window($chroma);
    return unless defined $transcript_window;
    
    my $chrom_for_window;
    my $last_pos;
    for my $break_anno ( $anno{breakpoint_a}, $anno{breakpoint_b} ) {
        my $chrom = $break_anno->{chromosome};
        my $pos = $break_anno->{position};
        
        $transcript_window = $self->_get_transcript_window($chrom) unless $chrom_for_window and $chrom_for_window eq $chrom;

        $self->error_message("couldn't get transcript window for chromosome $chrom") and return unless defined $transcript_window;
        if ($type eq "DEL" and $chrom_for_window){
            my @deleted_transcripts = $transcript_window->scroll($last_pos+1, $pos-1);
            my %genes;
            for my $t (@deleted_transcripts) {
                my $gene = $t->gene;
                my $name = $gene->name() if defined $gene;
                $genes{$name}++ if defined $name;
            }
            $anno{'deleted_genes'} = [keys %genes];
        }

        my @transcripts = $transcript_window->scroll($pos);
        $break_anno->{'transcripts'} = @transcripts? \@transcripts : undef;
        $chrom_for_window = $chrom;
        $last_pos = $pos;
    }
    
    return %anno;
}

sub gene_names_from_transcripts {
    # Given a list of transcripts, returns a hash of the unique corresponding genes and the number of times those genes appear
    my $self = shift;
    my @transcripts = @_;
    return unless @transcripts;
    my %genes;
    for my $t (@transcripts) {
        my $gene = $t->gene;
        my $name = $gene->name() if defined $gene;
        $genes{$name}++ if defined $name;
    }
    return %genes;
}

sub transcripts_at_position {
    # Returns all transcripts at the given position (or range) on the given chromosome
    my $self = shift;
    my ($chrom_name, $start, $stop) = @_;
    unless (defined $chrom_name and defined $start) {
        return;
    }

    if (defined $stop and $start > $stop) {
        $self->error_message("Start position must be less than stop position");
        return;
    }

    my $transcript_window = $self->_get_transcript_window($chrom_name);
    return unless $transcript_window;
    return $transcript_window->scroll($start, $stop);
}

1;

