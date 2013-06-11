package Genome::InstrumentData::Gatk::BaseWithKnownIndels;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Gatk::BaseWithKnownIndels {
    is => 'Genome::InstrumentData::Gatk::Base',
    has_optional_param => [
        known_indels => {# PROVIDES indels_vcf
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
    ],
};

sub known_indels_vcfs {
    my $self = shift;

    return $self->{_indels_vcfs} if $self->{_indels_vcfs};

    my @indels_vcfs;
    for my $known_indel ( $self->known_indels ) {
        my $indel_result = $known_indel->indel_result;
        if ( not $indel_result ) {
            $self->error_message('No indel result for known indel! '.$known_indel->__display_name__);
            return;
        }
        my $source_indels_vcf = $indel_result->path.'/indels.hq.vcf';
        if ( not -s $source_indels_vcf ) {
            $self->error_message('No indels vcf (indels.hq.vcf) in indel result path! '.$indel_result->path);
            return;
        }
        my $indels_vcf = $self->_tmpdir.'/'.$indel_result->id.'.vcf';
        symlink($source_indels_vcf, $indels_vcf);
        push @indels_vcfs, $indels_vcf;
    }

    $self->{_indels_vcfs} = \@indels_vcfs;
    return $self->{_indels_vcfs};
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $known_indels_vcfs = $self->known_indels_vcfs;
    return if not $known_indels_vcfs; # undef on error

    return $self;
}

1;

