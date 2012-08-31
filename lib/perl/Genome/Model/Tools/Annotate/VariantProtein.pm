package Genome::Model::Tools::Annotate::VariantProtein;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Annotate::VariantProtein {
    is => ['Command'],
    has_input => [
        input_tsv_file => {
            is => 'Text',
            doc => 'A tab separated input file from the annotator',
        },
        output_tsv_file => {
            is => 'Text',
            doc => 'A tab separated output file with the amino acid sequence, wiletype and mutant.',
        },
    ],
    has_optional_input => [
        anno_db => {
            is => 'Text',
            doc => 'The name of the annotation database.  Example: NCBI-human.combined-annotation',
        },
        version => {
            is => 'Text',
            doc => 'The version of the annotation database. Example: 54_36p_v2',
        },
    ],
};

sub create {
    my $class = shift;
    my $self  = $class->SUPER::create(@_);
    unless ($self) { return; }
    if (defined($self->anno_db) || defined($self->version)) {
        unless (defined($self->anno_db) && defined($self->version)) {
            die('Please define both anno_db and version!');
        }
    }
    return $self;
}

sub execute {
    my $self = shift;
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input =>  $self->input_tsv_file,
        separator => "\t",
    );
    my $headers = $reader->headers;
    my @new_headers = @{$headers};
    push @new_headers, 'wildtype_amino_acid_sequence';
    #push @new_headers, 'mutant_amino_acid_sequence';
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_tsv_file,
        separator => "\t",
        headers => \@new_headers,
    );
    my %binned_by_chr;
    my %sources;
    my %versions;
    my %species;
    while (my $data = $reader->next) {
        push @{$binned_by_chr{$data->{chromosome_name}}}, $data;
        unless ($data->{transcript_species} eq '-') {
            $species{$data->{transcript_species}} = 1;
        }
        unless ($data->{transcript_source} eq '-') {
            $sources{$data->{transcript_source}} = 1;
        }
        unless ($data->{transcript_version} eq '-') {
            $versions{$data->{transcript_version}} = 1;
        }
    };
    my $anno_db = $self->anno_db;
    unless ($anno_db) {
        my @sources = keys %sources;
        my $source;
        if (@sources ne 1) {
            $self->warning_message('Attempting to use the combined-annotation db.  Multiple sources found: '. Data::Dumper::Dumper(@sources));
            $source = 'combined-annotation';
        } else {
            $source = $sources[0];
        }
        my @species = keys %species;
        if (@species ne 1) {
            $self->error_message('More than one species found: '. Data::Dumper::Dumper(@species));
            die($self->error_message);
        }
        my $species = $species[0];
        $anno_db = 'NCBI-'. $species .'.'. $source;
    }
    my $model = Genome::Model->get(name => $anno_db);
    unless ($model) {
        $self->error_message('Failed to find annotation model by name '. $anno_db);
        die($self->error_message);
    }
    my $version = $self->version;
    unless ($version) {
        my @versions = keys %versions;
        if (@versions ne 1) {
            $self->error_message('Multiple versions found: '. Data::Dumper::Dumper(@versions));
            die($self->error_message);
        }
        $version = $versions[0];
    }
    my $build = $model->build_by_version($version);
    unless ($build) {
        $self->error_message('Failed to find annotation build by version '. $version);
        die($self->error_message);
    }
    for my $chr (keys %binned_by_chr) {
        my $ti = $build->transcript_iterator(chrom_name => $chr);
        my $transcript_window =  Genome::Utility::Window::Transcript->create(iterator => $ti);
        for my $data (sort {$a->{start} <=> $b->{start}} @{$binned_by_chr{$chr}}) {
            for my $t ($transcript_window->scroll($data->{start},$data->{stop})){
                unless ($t->transcript_name eq $data->{transcript_name}) { next; }
                my $protein = $t->protein;
                if ($protein) {
                    $data->{wildtype_amino_acid_sequence} = $protein->amino_acid_seq;
                }
            }
            unless (defined($data->{wildtype_amino_acid_sequence})) {
                if ($data->{amino_acid_change} eq 'NULL' || $data->{amino_acid_change} eq '-') {
                    $data->{wildtype_amino_acid_sequence} = $data->{amino_acid_change};
                } else {
                    $self->warning_message('Failed to find protein for variant: '. Data::Dumper::Dumper($data));
                    $data->{wildtype_amino_acid_sequence} = '';
                }
            }
            $writer->write_one($data);
        }
    }
    return 1;
}

1;
