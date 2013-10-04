package Genome::Gene;
#:adukes short term: move data directory into id_by, but this has to be done in parallel w/ rewriting all file-based data sources.  It might be better to wait until long term: custom datasource that incorporates data_dir, possibly species/source/version, eliminating the need for these properties in the id, and repeated multiple times in the files

use strict;
use warnings;

use Genome;

class Genome::Gene {
    type_name => 'genome gene',
    table_name => 'GENE',
    id_by => [
        gene_id => { 
            is => 'Text' 
        },
        species => { is => 'Text',
            is_optional => 1,
        },
        source => { is => 'Text',
            is_optional => 1,
        },
        version => { is => 'Text',
            is_optional => 1,
        },
    ],
    has => [
        hugo_gene_name => { 
            is => 'Text',
            is_optional => 1,
        },
        strand => {
            is => 'Text',
            valid_values => ['+1', '-1', 'UNDEF'],
        },
        data_directory => {
            is => "Path",
        },
        reference_build_id => {
            is => 'Text',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
        },
    ],
    has_many => [
        transcripts => { 
            calculate_from => [qw/ id data_directory reference_build_id /],
            calculate => q|
                Genome::Transcript->get(gene_id => $id,  data_directory => $data_directory, reference_build_id => $reference_build_id);
            |,
        },
        external_ids => { 
            calculate_from => [qw/ id data_directory reference_build_id /],
            calculate => q|
                Genome::ExternalGeneId->get(gene_id => $id, data_directory => $data_directory, reference_build_id => $reference_build_id);
            |,
        },
    ],
    schema_name => 'files',
    data_source => 'Genome::DataSource::Genes',
};

sub name
{
    my ($self, $source) = @_;

    my @egis;
    
    if ( $source )
    {
        if ( $source eq "genbank")
        {
            $source = 'entrez';
        }
        @egis = grep { $_->id_type() eq $source } $self->external_ids();
    }
    else
    {
        my $name = $self->hugo_gene_name;

        return $name if $name;
        @egis = $self->external_ids;
    }

    unless ($egis[0]) {
        return '';
    }

    return $egis[0]->id_value;
}


#- EXPRESSIONS -#
sub expressions_by_intensity
{
    my $self = shift;

    # Sort by decrementing intensity
    my @expressions = sort { $b->expression_intensity <=> $a->expression_intensity }
                           $self->expressions;
    return @expressions;
}

sub chrom_name {
    my $self = shift;
    my @t = $self->transcripts;
    my $chrom_name;
    for my $t (@t) {
        unless ($chrom_name) {
            $chrom_name = $t->chrom_name;
        } else {
            if ($chrom_name ne $t->chrom_name) {
                die('Expected chrom '. $chrom_name .' but found chrom '. $t->chrom_name .' for transcript '. $t->transcript_name . ' of gene '. $self->name);
            }
        }
    }
    return $chrom_name;
}

sub gene_start {
    my $self = shift;
    my @t = $self->transcripts;
    my $gene_start;
    for my $t (@t) {
        unless ($gene_start) {
            $gene_start = $t->transcript_start;
        } else {
            if ($gene_start > $t->transcript_start) {
                $gene_start = $t->transcript_start;
            }
        }
    }
    return $gene_start;
}

sub gene_stop {
    my $self = shift;
    my @t = $self->transcripts;
    my $gene_stop;
    for my $t (@t) {
        unless ($gene_stop) {
            $gene_stop = $t->transcript_stop;
        } else {
            if ($gene_stop < $t->transcript_stop) {
                $gene_stop = $t->transcript_stop;
            }
        }
    }
    return $gene_stop;
}

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

sub bed_string {
    my $self = shift;
    # BED entries should only be written per sub-structure unless BED12 format is adopted, even then it should be per transcript
    return;
    # BED format uses zero-based start positions
    my $bed_string = $self->chrom_name ."\t". ($self->gene_start - 1) ."\t". $self->gene_stop ."\t". $self->name ."\t0\t". $self->strand_string;
    return $bed_string ."\n";
}

sub _base_gff_string {
    my $self = shift;
    return $self->chrom_name ."\t". $self->source .'_'. $self->version ."\tgene\t". $self->gene_start ."\t". $self->gene_stop ."\t.\t". $self->strand_string ."\t.";
}

sub gff_string {
    my $self = shift;
    return $self->_base_gff_string ."\t". $self->name ."\n";
}

sub gff3_string {
    my $self = shift;
    return $self->_base_gff_string ."\tID=". $self->gene_id .'; NAME='. $self->name .';' ."\n";
}

sub gtf_string {
    my $self = shift;
    return undef;
}

1;

