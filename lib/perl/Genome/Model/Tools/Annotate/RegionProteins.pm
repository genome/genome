package Genome::Model::Tools::Annotate::RegionProteins;

use strict;
use warnings;

use Genome;

use Genome::Model::Tools::RefCov::ROI::Bed;

class Genome::Model::Tools::Annotate::RegionProteins {
    is => ['Command'],
    has_input => [
        anno_db => {
            default_value => 'NCBI-human.combined-annotation',
        },
        version => {
            default_value => '54_36p_v2',
        },
        bed_file => {

        },
        protein_fasta_file => {

        },
    ],
};

sub execute {
    my $self = shift;
    my $bed = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->bed_file,
        make_objects => 1,
        load_all => 1,
    );
    my @chromosomes = $bed->chromosomes;
    my %region_protein_fastas;
    # TODO: use bioperl to write in true fasta format(is there one... I guess 60bp per line rather than one line per)
    my $fh = IO::File->new($self->protein_fasta_file,'w');
    for my $chromosome (@chromosomes) {
        my @regions = sort{ $a->start <=> $b->start } $bed->chromosome_regions($chromosome);
        my $ti = Genome::Model->get(name => $self->anno_db)->build_by_version($self->version)->transcript_iterator(chrom_name => $chromosome);
        my $transcript_window =  Genome::Utility::Window::Transcript->create(iterator => $ti);
        for my $region (@regions) {
            for my $t ($transcript_window->scroll($region->start,$region->end)){
                my $protein = $t->protein;
                if ($protein) {
                    $region_protein_fastas{$protein->protein_name}{fasta} = $protein->amino_acid_seq;
                    $region->{_id} = $region->{chrom} .':'. $region->{start} .'-'. $region->{end};
                    push @{$region_protein_fastas{$protein->protein_name}{regions}}, $region;
                }
            }
        }
        for my $protein (keys %region_protein_fastas) {
            my @regions = @{$region_protein_fastas{$protein}{regions}};
            my $region_ids = join(',', map { $_->{_id} } @regions);
            print $fh '>'. $protein .' '. $region_ids ."\n";
            print $fh $region_protein_fastas{$protein}{fasta} ."\n";
        }
    }
    $fh->close;
    return 1;
}

1;
