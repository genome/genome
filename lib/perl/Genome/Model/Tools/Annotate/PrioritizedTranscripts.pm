package Genome::Model::Tools::Annotate::PrioritizedTranscripts;

use strict;
use warnings;

use Genome;
use Genome::Info::AnnotationPriorities;

use Carp qw(confess);
use IO::File;

class Genome::Model::Tools::Annotate::PrioritizedTranscripts {
    is => 'Genome::Model::Tools::Annotate',
    has => [
        regions_file => {
            is => 'Path',
            doc => 'Path to a regions files, with chromosome, start, and stop columns (tab delimited)',
        },
        output_file => {
            is => 'Path',
            doc => 'File to write output to',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Reference sequence build with the coordinates used in regions file',
        }
    ],
    has_optional => [
        reference_transcripts => {
            is => 'String',
            default => 'NCBI-human.combined-annotation/54_36p_v2',
            doc => 'Annotation version to use, format is <annotation model name>/<version>',
        },
        filter => {
            is => 'String',
            default => 'none',
            valid_values => ['none', 'gene', 'top'],
            doc => 'Similar to filters for annotator. If none, then all transcripts are put in output. ' .
                   'If set to gene, only the highest priority transcript per gene is placed in output. ' .
                   'If set to top, only the top transcript per given region is placed in output.',
        },
    ],
};

sub help_detail {
    return <<EOS 
This tool is given a bunch of regions and spits out a prioritized list of transcripts and genes for
each region. Prioritization is done identically to how the annotator does it, minus the variant.
EOS
}

sub help_brief {
    return "Finds a prioritized list of transcripts for each given region";
}

sub help_synopsis {
    return "Finds a prioritized list of transcripts for each given region";
}

sub execute {
    my $self = shift;
    my $reference_build_id = $self->reference_build->id;

    unless (-e $self->regions_file) {
        confess "No file found at " . $self->regions_file;
    }
    if (-s $self->output_file) {
        confess "File exists and has size at " . $self->output_file . ", will not overwrite!";
    }
    
    my ($name, $version) = split("/", $self->reference_transcripts);
    unless (defined $name and defined $version) {
        confess "Invalid reference transcripts string: " . $self->reference_transcripts;
    }

    my $model = Genome::Model::ImportedAnnotation->get(name => $name);
    unless ($model) {
        confess "Could not get imported annotation model with name $name";
    }

    my $build = $model->build_by_version($version);
    unless ($build) {
        confess "Could not get build version $version from model $name";
    }

    my @data_dirs = $build->determine_data_directory;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        headers => ['chrom_name', 'start', 'stop'],
        separator => "\t",
        is_regex => 1,
        ignore_extra_columns => 1,
        input => $self->regions_file,
    );
    unless ($reader) {
        confess "Could not make separated value reader for " . $self->regions_file;
    }

    my $output_fh = IO::File->new($self->output_file, "w");

    while (my $line = $reader->next) {
        my @transcripts = Genome::Transcript->get(
            chrom_name => $line->{chrom_name},
            transcript_start => {operator => '<=', value => $line->{stop}},
            transcript_stop => {operator => '>=', value => $line->{start}},
            reference_build_id => $reference_build_id,
            data_directory => \@data_dirs,
        );
        next unless @transcripts;

        my @prioritized_transcripts = $self->prioritize_by_filter(@transcripts);
        for my $t (@prioritized_transcripts) {
            $output_fh->print(join("\t", $line->{chrom_name}, $line->{start}, $line->{stop}, $t->transcript_start, 
                    $t->transcript_stop, $t->gene_name, $t->transcript_name) . "\n");
        }
    }

    $output_fh->close;
    return 1;
}

sub highest_priority_error {
    my ($self, $transcript) = @_;
    my %transcript_error_priorities = Genome::Info::AnnotationPriorities->transcript_error_priorities;
    my $error_string = $transcript->{transcript_error};
    my @errors = map { $transcript_error_priorities{$_} } split(":", $error_string);
    my @sorted_errors = sort { $b <=> $a } @errors;
    return $sorted_errors[0];
}

sub prioritize_by_filter {
    my ($self, @transcripts) = @_;
    my $filter = $self->filter;

    my @prioritized_transcripts = $self->prioritize_transcripts(@transcripts);
    if ($filter eq 'none') {
        return @prioritized_transcripts;
    }
    elsif ($filter eq 'top') {
        return shift @prioritized_transcripts;
    }
    elsif ($filter eq 'gene') {
        my %top_per_gene;
        for my $t (@prioritized_transcripts) {
            $top_per_gene{$t->gene_name} = $t unless exists $top_per_gene{$t->gene_name};
        }
        return values %top_per_gene;
    }
    else {
        confess "Invalid filter value $filter!";
    }
}

sub prioritize_transcripts {
    my ($self, @transcripts) = @_;

    my %transcript_source_priorities = Genome::Info::AnnotationPriorities->transcript_source_priorities;
    my %transcript_status_priorities = Genome::Info::AnnotationPriorities->transcript_status_priorities;

    my @prioritized_transcripts = sort {
        $self->highest_priority_error($a) <=> $self->highest_priority_error($b) ||
        $transcript_source_priorities{$a->{source}} <=> $transcript_source_priorities{$b->{source}} ||
        $transcript_status_priorities{$a->{transcript_status}} <=> $transcript_status_priorities{$b->{transcript_status}} ||
        $b->{amino_acid_length} <=> $a->{amino_acid_length} ||
        $a->{transcript_name} cmp $b->{transcript_name} } @transcripts;
    
    return @prioritized_transcripts;
}

1;

