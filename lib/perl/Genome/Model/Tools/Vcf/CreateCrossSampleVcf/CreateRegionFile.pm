package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateRegionFile;

use strict;
use warnings;

use Genome;
use Genome::File::Vcf::Reader;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateRegionFile {
    is  => 'Genome::Command::Base',
    has_input => [
        segregating_sites_vcf_file => {
            is => 'File',
        },
    ],
    has_calculated_output => [
        region_file => {
            is => 'File',
            calculate => q("$segregating_sites_vcf_file.converted"),
            calculate_from => 'segregating_sites_vcf_file',
        },
    ],
};

sub execute {
    my $self = shift;

    my $ofh = Genome::Sys->open_file_for_writing($self->region_file);

    my $reader = Genome::File::Vcf::Reader->new($self->segregating_sites_vcf_file);
    while (my $entry = $reader->next) {
        my $alts = join(',', $entry->{alternate_alleles});
        $ofh->print( join("\t", $entry->{chrom}, $entry->{position},
                                $entry->{position}, $alts) );
    }
    $ofh->close();

    return 1;
}

1;
