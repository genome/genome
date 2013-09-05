package Genome::Site::TGI::IlluminaGenotyping;

use strict;
use warnings;
use Genome;

class Genome::Site::TGI::IlluminaGenotyping {
    table_name => 'ILLUMINA_GENOTYPING',
    data_source => 'Genome::DataSource::Dwrac',
    id_by => [
        seq_id => { is => 'Number' },
    ],
    has => [
        name => { is => 'Text' , column_name => 'DNA_NAME'},
        source_barcode => { is => 'Text' },
        well => { is => 'Text' },
        bead_chip_barcode => { is => 'Text' },
        organism_sample_id => { is => 'Number' },
    ],
    has_optional => [
        replicate_dna_name => { is => 'Text' },
        call_rate => { is => 'Number' },
        replicate_error_rate => { is => 'Number' },
        status => { is => 'Text' },
        analysis_name => { is => 'Text' },
        genome_studio_version => { is => 'Text' },
        sample_concordance => { is => 'Number' },
        replicate_seq_id => { is => 'Number' },
        custom_num_fail_cutoff => { is => 'Number' },
        custom_cutoff => { is => 'Number' },
        num_fail_cutoff => { is => 'Number' },
        cutoff => { is => 'Number' },
        num_of_probe => { is => 'Number' },
    ],
};


sub __display_name__ {
    my $self = shift;
    return $self->source_barcode . ' for sample ' . $self->name;
}


sub meets_default_criteria {
    my $self = shift;

    return unless ($self->status eq 'pass');

    my @snp_concordance = Genome::Site::TGI::SnpConcordance->get(seq_id => $self->seq_id); # will return at least one external comparison and sometimes the internal comparison
    push @snp_concordance, Genome::Site::TGI::SnpConcordance->get(replicate_seq_id => $self->seq_id); # returns the internal comparison when the above doesn't

    my $list_has_external_comparison = grep { $_->is_external_comparison } @snp_concordance;
    my $list_has_internal_comparison = grep { $_->is_internal_comparison } @snp_concordance;
    return unless ($list_has_internal_comparison && $list_has_external_comparison);

    my @good_snp_concordance = grep { $_->match_percent && $_->match_percent > 90 } @snp_concordance;
    return (@snp_concordance == @good_snp_concordance);
}


1;

