package Genome::Model::Build::RunsDV2;

use strict;
use warnings FATAL => 'all';

class Genome::Model::Build::RunsDV2 {
    is_abstract => 1,
};

sub get_detailed_indels_vcf {
    my $self = shift;
    return File::Spec->join($self->variants_directory, "indels.detailed.vcf.gz");
}

sub get_detailed_snvs_vcf {
    my $self = shift;
    return File::Spec->join($self->variants_directory, "snvs.detailed.vcf.gz");
}

sub get_indels_vcf {
    my $self = shift;
    return File::Spec->join($self->variants_directory, "indels.vcf.gz");
}

sub get_snvs_vcf {
    my $self = shift;
    return File::Spec->join($self->variants_directory, "snvs.vcf.gz");
}

sub variants_directory {
    my $self = shift;

    my $expected_directory = File::Spec->join($self->data_directory, 'variants');

    if (-d $expected_directory) {
        return $expected_directory;
    } else {
        #for compatibility with previously existing builds
        for my $dir_name qw(snp_related_metrics  sam_snp_related_metrics  maq_snp_related_metrics  var-scan_snp_related_metrics) {
            my $dir = File::Spec->join($self->data_directory, $dir_name);
            return $dir if -d $dir;
        }
        die $self->error_message("Variants directory does not exist at $expected_directory");
    }
}
