package Genome::Model::Tools::Relationship::SequencingQcResult;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Relationship::SequencingQcResult {
    is => 'Genome::SoftwareResult::DiskAllocationStaged',
    has_param => [
        bams => { is => 'Text', is_many => 1, },
        snp_files => { is =>'Text', is_many=>1, },
        reference_fasta=> { is => 'Text', },
        use_sites_found_in_any_sample => { is=>"Boolean", },
        use_1000genomes_build37_sites => { is=>"Boolean", },
        use_this_bed_file=> { is=>"Text", },
        min_coverage => { is => 'Number',},
        ped_file=> { is_optional=>1, },
        parent_relationship_cutoff => { is => 'Number', },
    ],
    has_metric => [
        qc_status => { is => 'Text', valid_values => ['Pass', 'Fail'] },
    ],
};

sub _generate_result {
    my ($self, $staging_directory) = @_;

    my @bams = $self->bams;
    my @snp_files = $self->snp_files;

    my $command = Genome::Model::Tools::Relationship::SequencingQc->create(
        output_dir => $staging_directory,
        use_sites_found_in_any_sample => $self->use_sites_found_in_any_sample,
        use_1000genomes_build37_sites => $self->use_1000genomes_build37_sites,
        use_this_bed_file => $self->use_this_bed_file,
        bams => \@bams,
        snp_files => \@snp_files,
        reference_fasta => $self->reference_fasta,
        min_coverage => $self->min_coverage,
        ped_file => $self->ped_file,
        parent_relationship_cutoff => $self->parent_relationship_cutoff,
    );

    my $return = $command->_generate_data();
    if ($return eq "Pass" || $return eq "Fail") {
        $self->qc_status($return);
    }

    return $return;
}

1;

