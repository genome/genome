package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump;

use strict;
use warnings;

use JSON;
use UR;


my @_PROPERTY_NAMES = qw(
    backfilled_vcf
    bam_file
    pileup_output_file
    sample
    vcf_file
);

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump {
    id_by => \@_PROPERTY_NAMES,
    has => \@_PROPERTY_NAMES,

    data_source => 'UR::DataSource::Default',
};

sub _serialize {
    my $js = new JSON->allow_nonref;
    return $js->canonical->encode(@_);
}

sub _deserialize {
    my $js = new JSON->allow_nonref;
    return $js->decode(@_);
}

sub Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump::Type::get_composite_id_resolver {
    my $class_meta = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->__meta__;
    my @ids = $class_meta->id_property_names;

    return sub {
        my %hash;
        for (my $i = 0; $i < scalar(@ids); $i++) {
            $hash{$ids[$i]} = $_[$i];
        }

        return _serialize(\%hash);
    }
}

sub _get_values_from_id {
    my $id = shift;

    my @property_order;
    if (@_) {
        @property_order = @_;
    } else {
        my $class_meta = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->__meta__;
        @property_order = $class_meta->id_property_names;
    }

    my $hash = _deserialize($id);

    return map {$hash->{$_}} @property_order;
}

sub Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump::Type::get_composite_id_decomposer {
    return sub {
        return _get_values_from_id(@_);
    };
}

sub __load__ {
    my ($class, $bx, $headers) = @_;

    my $id = $bx->value_for('id', $headers);
    my $rows = [
         [_get_values_from_id($id), $id],
    ];
    return $headers, $rows;
}

1;
