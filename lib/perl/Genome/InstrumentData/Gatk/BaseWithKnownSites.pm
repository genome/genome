package Genome::InstrumentData::Gatk::BaseWithKnownSites;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Gatk::BaseWithKnownSites {
    is => 'Genome::InstrumentData::Gatk::Base',
    has_optional_param => [
        known_sites => {# PROVIDES known_sites_vcf
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
    ],
};

sub known_sites_vcfs {
    my $self = shift;

    return $self->{_known_sites_vcfs} if $self->{_known_sites_vcfs};

    my @known_sites_vcfs;
    for my $known_site ( $self->known_sites ) {
        my $indel_result = $known_site->indel_result;
        if ( not $indel_result ) {
            $self->error_message('No indel result for known indel! '.$known_site->__display_name__);
            return;
        }
        my $source_known_sites_vcf = $indel_result->path.'/indels.hq.vcf';
        if ( not -s $source_known_sites_vcf ) {
            $self->error_message('No indels vcf (indels.hq.vcf) in indel result path! '.$indel_result->path);
            return;
        }
        my $known_sites_vcf = $self->_tmpdir.'/'.$indel_result->id.'.vcf';
        symlink($source_known_sites_vcf, $known_sites_vcf);
        push @known_sites_vcfs, $known_sites_vcf;
    }

    $self->{_known_sites_vcfs} = \@known_sites_vcfs;
    return $self->{_known_sites_vcfs};
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $known_sites_vcfs = $self->known_sites_vcfs;
    return if not $known_sites_vcfs; # undef on error

    return $self;
}

1;

