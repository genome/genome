package Genome::Model::Tools::Picard::GenotypeConcordance;

use strict;
use warnings;

use Genome;

my %STATES_TO_IGNORE = (
    'MISSING' => 1,
    'NO_CALL' => 1,
    'VC_FILTERED' => 1,
    'LOW_DP' => 1,
    'LOW_GQ' => 1,
    'GT_FILTERED' => 1,
);

class Genome::Model::Tools::Picard::GenotypeConcordance {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        truth_vcf => {
            is  => 'String',
            doc => 'The VCF containing the truth sample',
        },
        call_vcf => {
            is  => 'String',
            doc => 'The VCF containing the call sample',
        },
        output => {
            is  => 'String',
            doc => 'Basename for the two metrics files that are to be written. Resulting files will be <OUTPUT>.summary_metrics.txt and <OUTPUT>.detailed_metrics.txt.',
        },
        truth_sample => {
            is  => 'String',
            doc => 'The name of the truth sample within the truth VCF',
        },
        call_sample => {
            is => 'String',
            doc => 'The name of the call sample within the call VCF',
        },
        intervals => {
            is => 'String',
            is_many => 1,
            doc => 'One or more interval list files that will be used to limit the genotype concordance. This option may be specified 0 or more times.',
            is_optional => 1,
        },
        intersect_intervals => {
            is => 'Boolean',
            default_value => 1,
            doc => 'If true, multiple interval lists will be intersected. If false multiple lists will be unioned.',
            is_optional => 1,
        },
        min_gq => {
            is => 'Integer',
            default_value => 0,
            doc => 'Genotypes below this genotype quality will have genotypes classified as LowGq.',
            is_optional => 1,
        },
        min_dp => {
            is => 'Integer',
            default_value => 0,
            doc => 'Genotypes below this depth will have genotypes classified as LowDp.',
            is_optional => 1,
        },
        output_all_rows => {
            is => 'Boolean',
            default_value => 0,
            doc => 'If true, output all rows in detailed statistics even when count == 0. When false only output rows with non-zero counts.',
            is_optional => 1,
        },
        use_vcf_index => {
            is => 'Boolean',
            default_value => 0,
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

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/GenotypeConcordance.jar';
    unless (-e $jar_path) {
        $self->error_message('Failed to find '. $jar_path .'!  GenotypeConcordance is only available in Picard version 1.122 or greater.  The version requested is '. $self->use_version);
        die($self->error_message);
    }

    my $gc_cmd = $jar_path .' picard.vcf.GenotypeConcordance TRUTH_VCF='. $self->truth_vcf .' CALL_VCF='.$self->call_vcf .' OUTPUT='. $self->output .' TRUTH_SAMPLE='. $self->truth_sample .' CALL_SAMPLE='. $self->call_sample;
    if (defined($self->min_dp)) {
        $gc_cmd .= ' MIN_DP='. $self->min_dp;
    }
    if (defined($self->min_gq)) {
        $gc_cmd .= ' MIN_GQ='. $self->min_gq;
    }
    if (defined($self->intervals)) {
        $gc_cmd .= ' INTERVALS='. $self->intervals;
    }
    if (defined($self->intersect_intervals)) {
        if ($self->intersect_intervals) {
            $gc_cmd .= ' INTERSECT_INTERVALS=TRUE';
        }
    }
    if (defined($self->output_all_rows)) {
        if ($self->output_all_rows) {
            $gc_cmd .= ' OUTPUT_ALL_ROWS=TRUE';
        }
    }
    if (defined($self->use_vcf_index)) {
        if ($self->use_vcf_index) {
            $gc_cmd .= ' USE_VCF_INDEX=TRUE';
        }
    }
    $self->run_java_vm(
        cmd => $gc_cmd,
        input_files => [$self->truth_vcf,$self->call_vcf],
        skip_if_output_is_present => 0,
    );
    
    return 1;
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

sub ignore_state {
    my $class = shift;
    my $state = shift;
    if ( defined($STATES_TO_IGNORE{$state}) ) {
        return $STATES_TO_IGNORE{$state};
    }
    return 0;
}

1;
