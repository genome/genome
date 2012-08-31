#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Storable;
use Test::More skip_all => 'out of date, need to update without using models';

my $build    = Genome::Model::ImportedAnnotation->get(name => 'NCBI-human.combined-annotation')->build_by_version(0);
my $iterator = $build->transcript_iterator;

my $storable_data_file = ($ENV{GENOME_TEST_INPUTS} . '/Genome-Transcript/annot-var-5.stor');

my $sd = retrieve("$storable_data_file");
my %data;

for ( 1 .. 5 )
{
    my $transcript       = $iterator->next;
    my $transcript_id    = $transcript->id;
    my @substructures    = $transcript->sub_structures;
    my $gene             = $transcript->gene;
    my $protein          = $transcript->protein;
    $data{$_}{gene} = $gene;
    $data{$_}{substructures} = \@substructures;
    $data{$_}{transcript_id} = $transcript->id;
    $data{$_}{protein} = $protein;

}

foreach my $key ( 1 .. 5 )
{

    is( $data{$key}{transcript_id},
        $sd->{$key}->{transcript_id},
        'transcript id'
    );
    is( $data{$key}{protein}{protein_id},
        $sd->{$key}->{protein}->{protein_id},
        'protein id'
    );
    is( $data{$key}{protein}{protein_name},
        $sd->{$key}->{protein}->{protein_name},
        'protein name'
    );
    is( $data{$key}{protein}{amino_acid_seq},
        $sd->{$key}->{protein}->{amino_acid_seq},
        'amino acid seq'
    );
    is( $data{$key}{gene}{gene_id},
        $sd->{$key}->{gene}->{gene_id},
        'gene id'
    );
    is( $data{$key}{gene}{hugo_gene_name},
        $sd->{$key}->{gene}->{hugo_gene_name},
        'hugo gene name'
    );

    foreach my $item ( 0 .. $#{ $sd->{$key}->{substructures} } )
    {
        is( $data{$key}{substructures}[$item]{structure_type},
            $sd->{$key}->{substructures}->[$item]->{structure_type},
            'transcript substructure types match'
        );

        is( $data{$key}{substructures}[$item]{structure_start},
            $sd->{$key}->{substructures}->[$item]->{structure_start},
            'transcript substructure starts match'
        );

        is( $data{$key}{substructures}[$item]{structure_end},
            $sd->{$key}->{substructures}->[$item]->{structure_end},
            'transcript substructure ends match'
        );

        is( $data{$key}{substructures}[$item]{ordinal},
            $sd->{$key}->{substructures}->[$item]->{ordinal},
            'transcript substructure ordinal match'
        );
    }

}

# _fin_
