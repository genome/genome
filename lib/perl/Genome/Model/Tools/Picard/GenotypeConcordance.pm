package Genome::Model::Tools::Picard::GenotypeConcordance;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::GenotypeConcordance {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        truth_vcf => {
            is  => 'String',
            picard_param_name => 'TRUTH_VCF',
            doc => 'The VCF containing the truth sample',
        },
        call_vcf => {
            is  => 'String',
            picard_param_name => 'CALL_VCF',
            doc => 'The VCF containing the call sample',
        },
        output => {
            is  => 'String',
            picard_param_name => 'OUTPUT',
            doc => 'Basename for the two metrics files that are to be written. Resulting files will be <OUTPUT>.summary_metrics.txt and <OUTPUT>.detailed_metrics.txt.',
        },
        truth_sample => {
            is  => 'String',
            picard_param_name => 'TRUTH_SAMPLE',
            doc => 'The name of the truth sample within the truth VCF',
        },
        call_sample => {
            is => 'String',
            picard_param_name => 'CALL_SAMPLE',
            doc => 'The name of the call sample within the call VCF',
        },
        intervals => {
            is => 'String',
            is_many => 1,
            picard_param_name => 'INTERVALS',
            doc => 'One or more interval list files that will be used to limit the genotype concordance. This option may be specified 0 or more times.',
            is_optional => 1,
        },
        intersect_intervals => {
            is => 'Boolean',
            default_value => 1,
            picard_param_name => 'INTERSECT_INTERVALS',
            doc => 'If true, multiple interval lists will be intersected. If false multiple lists will be unioned.',
            is_optional => 1,
        },
        min_gq => {
            is => 'Integer',
            default_value => 0,
            picard_param_name => 'MIN_GQ',
            doc => 'Genotypes below this genotype quality will have genotypes classified as LowGq.',
            is_optional => 1,
        },
        min_dp => {
            is => 'Integer',
            default_value => 0,
            picard_param_name => 'MIN_DP',
            doc => 'Genotypes below this depth will have genotypes classified as LowDp.',
            is_optional => 1,
        },
        output_all_rows => {
            is => 'Boolean',
            default_value => 0,
            picard_param_name => 'OUTPUT_ALL_ROWS',
            doc => 'If true, output all rows in detailed statistics even when count == 0. When false only output rows with non-zero counts.',
            is_optional => 1,
        },
        use_vcf_index => {
            is => 'Boolean',
            default_value => 0,
            picard_param_name => 'USE_VCF_INDEX',
            doc => 'If true, use the VCF index, else iterate over the entire VCF.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Calculates the concordance between genotype data for two samples in two different VCFs.';
}

sub help_detail {
    return <<EOS

Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#GenotypeConcordance

Calculates the concordance between genotype data for two samples in two different VCFs - one being considered the truth (or reference) the other being considered the call. The concordance is broken into separate results sections for SNPs and indels. Summary and detailed statistics are reported

Note that for any pair of variants to compare, only the alleles for the samples under interrogation are considered and MNP, Symbolic, and Mixed classes of variants are not included.

EOS
}

sub minimum_version_required { return 1.122; }
sub _jar_name {
    my $self = shift;
    if ($self) {
        if ($self->version_newer_than('1.123')) {
            return 'picard-'. $self->use_version .'.jar';
        }
    }
    return 'GenotypeConcordance.jar';
}
sub _java_class { 'picard.vcf.GenotypeConcordance'; }

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files => [ $self->truth_vcf, $self->call_vcf ],
        skip_if_output_is_present => 0,
    );
}

sub parse_file_into_metrics_hashref {
    my ($class,$metrics_file) = @_;

    my $is_fh = Genome::Sys->open_file_for_reading($metrics_file);

    my @headers;
    my %data;
    while (my $line = $is_fh->getline) {
        chomp($line);
        if ($line =~ /^## METRICS CLASS/) {
            my $next_line = $is_fh->getline;
            chomp($next_line);
            @headers = split("\t",$next_line);
            next;
        }
        if (@headers) {
            if ($line =~ /^\s*$/) {
                last;
            } else {
                my @values = split("\t",$line);
                $data{$values[0]}{$values[1]}{$values[2]}{$values[3]}{$values[4]} = $values[5];
            }
        }
    }
    return \%data;
}

1;
