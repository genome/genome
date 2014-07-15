package Genome::VariantReporting::BamReadcount::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::File::Vcf::Reader;
use Sys::Hostname;
use IPC::Run qw(run);

class Genome::VariantReporting::BamReadcount::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        aligned_bam_result_id => {
            is => 'Text',
            doc => 'The bam result used to calculate read counts',
        },
    ],
    has_param => [
        version => { is  => 'Version', },
        minimum_mapping_quality => { is => 'Integer', },
        minimum_base_quality => { is => 'Integer', },
        max_count => { is  => 'Integer', },
        per_library => { is  => 'Bool', },
        insertion_centric => { is  => 'Bool', },
    ],
};

sub output_filename {
    return 'bam-readcount-output.tsv';
}

sub aligned_bam_result {
    my $self = shift;
    return Genome::InstrumentData::AlignedBamResult->get($self->aligned_bam_result_id);
}

sub _run {
    my $self = shift;

    my $region_list = $self->make_region_file($self->input_vcf);

    # Sam::Readcount doesn't accept variant_type or test_name
    my %params = $self->param_hash;
    delete $params{variant_type};
    delete $params{test_name};

    Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $self->aligned_bam_result->bam_file,
        reference_fasta => $self->aligned_bam_result->reference_fasta,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
        region_list => $region_list,
        use_version => delete $params{version},
        %params,
    );

    return;
}

sub sample_name {
    my $self = shift;
    return $self->aligned_bam_result->sample_name;
}

sub make_region_file {
    my ($self, $vcf_file) = @_;

    my ($out_fh, $region_list) = Genome::Sys->create_temp_file();
    my $reader = Genome::File::Vcf::Reader->new($vcf_file);

    while (my $entry = $reader->next) {
        my %hash = sort_hash($entry);
        
        for my $p (sort{$a <=> $b} keys %hash) {
            $out_fh->print(join "\t", $entry->{chrom}, $p, $p."\n");
        }
    }

    $out_fh->close;
    return $region_list;
}

sub sort_hash {
    my $entry = shift;
    my %hash;
    my $pos = $entry->{position};

    if ($entry->has_insertion or $entry->has_substitution) {
        $hash{$pos}++;
    }
    if ($entry->has_deletion) {
        $hash{$pos+1}++;
    }
    return %hash;
}


