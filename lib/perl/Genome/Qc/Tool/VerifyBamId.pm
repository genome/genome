package Genome::Qc::Tool::VerifyBamId;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw(any);

class Genome::Qc::Tool::VerifyBamId {
    is => 'Genome::Qc::Tool::WithVariationListVcf',
    has => [
        default_genotype_vcf_file_id => {
            is => 'Text',
            is_optional => 1,
        },
        default_genotype_vcf_file => {
            is => 'Genome::SoftwareResult::StageableSimple::SingleFile',
            id_by => 'default_genotype_vcf_file_id',
        },
    ],
};


sub cmd_line {
    my $self = shift;

    my $cmd = $self->gmt_class->create($self->gmt_params);
    return $cmd->_get_cmd_list;
}

sub get_metrics {
    my $self = shift;

    my %metrics;
    for my $file_extension ($self->_file_extensions_to_parse) {
        my $file = $self->qc_metrics_file . ".$file_extension";
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            input => $file,
            separator => '\t',
            is_regex => 1,
        );

        while ( my $line = $reader->next ) {
            my $place_to_assign = \%metrics;
            for my $accumulation ('RG') {
                $place_to_assign = $place_to_assign->{$line->{$accumulation}} ||= {};
            }
            %$place_to_assign = %$line;
        }
    }
    return $self->_flatten_metrics_hash(\%metrics);
}

sub gmt_class {
    return 'Genome::Model::Tools::VerifyBamId';
}

sub qc_metrics_file_accessor {
    return 'out_prefix';
}

sub _non_metric_columns {
    return split(' ', '#SEQ_ID RG');
}

sub _file_extensions_to_parse {
    return qw(selfSM selfRG);
}

sub genotype_vcf_file {
    my $self = shift;
    if ($self->qc_genotype_vcf_file) {
        return $self->qc_genotype_vcf_file->file_path;
    }
    elsif ($self->default_genotype_vcf_file) {
        return $self->default_genotype_vcf_file->file_path;
    }
    else {
        $self->fatal_message("No qc genotype vcf file provided.");
    }
}

1;
