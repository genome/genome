package Genome::Model::Tools::Mutect;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;

my $DEFAULT_VERSION = '1.1.4';
my $MUTECT_BASE = 'muTect';
my $MUTECT_COMMAND = "$MUTECT_BASE.jar";

class Genome::Model::Tools::Mutect {
    is => ['Command'],
    has_input => [
        version => {
            is    => 'string',
            doc   => 'version of Mutect to use',
            default => $DEFAULT_VERSION,
        },
        normal_bam => {
            is => 'Text',
            doc => "BAM file for Normal Sample",
            is_optional => 0,
        },
        tumor_bam => {
            is => 'Text',
            doc => "BAM file for Tumor Sample",
            is_optional => 0,
        },
        reference => {
            is => 'String',
            doc => 'Reference sequence',
            is_optional => 0,
        },
        cosmic_vcf => {
            is => 'String',
            doc => 'VCF file containing sites from COSMIC',
            is_optional => 1,
        },
        dbsnp_vcf => {
            is => 'String',
            doc => 'VCF file containing sites from dbSNP',
            is_optional => 1,
        },
        intervals => {
            is => 'String',
            doc => 'Intervals to run mutect on. Should be in format chr:start-end',
            is_many => 1,
            is_optional => 1,
        },
        tumor_sample_name => {
            is => 'String',
            doc => 'name to use for tumor in output files',
            is_optional => 1,
        },
        normal_sample_name => {
            is => 'String',
            doc => 'name to use for normal in output files',
            is_optional => 1,
        },
        only_passing_calls => {
            is => 'Boolean',
            doc => 'whether or not to report non-passing calls',
            is_optional => 1,
        },
    ],
    has_output => [
        vcf => {
            is => 'String',
            doc => 'Optional output in VCF format. Not up to spec...',
            is_optional => 1,
        },
        coverage_file => {
            is => 'String',
            doc => 'Information (in wig format) about which bases were sufficiently covered',
            is_optional => 1,
        },
        output_file => {
            is => 'String',
            doc => 'Native output format file',
            is_optional => 0,
        },
    ],
    has_param => [
        lsf_resource => {
            default => '-M 16777216 rusage[mem=16384] select[type==LINUX64 & mem > 16384] span[hosts=1]',
        },
    ],
};

sub help_brief {
    "the Broad's somatic SNV caller"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

sub mutect_versions {
    my %MUTECT_VERSIONS = (
        'test' => "/gscuser/dlarson/mutect/muTect-1.1.4.jar", # This is temporary until this can be packaged...
        '1.1.4' => Genome::Sys->jar_path($MUTECT_BASE, '1.1.4'),
    );
    return %MUTECT_VERSIONS;
}

sub mutect_path {
    my $self = shift;
    return $self->path_for_mutect_version($self->version);
}

sub available_mutect_versions {
    my $self = shift;
    my %versions = $self->mutect_versions;
    return keys %versions;
}

sub path_for_mutect_version {
    my $class = shift;
    my $version = shift;
    my %versions = $class->mutect_versions;
    if (defined $versions{$version}) {
        return $versions{$version};
    }
    die('No path for mutect version '. $version);
}

sub default_mutect_version {
    my $class = shift;
    my %versions = $class->mutect_versions;
    die "default mutect version: $DEFAULT_VERSION is not valid" unless $versions{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my %versions = $self->mutect_versions;
    if(exists($versions{$version})){
        return 1;
    }
    return 0;
}
sub execute {
    my $self = shift;

    my $cmd = "java -Xmx5g -jar " . $self->path_for_mutect_version($self->version) . " -et NO_ET --analysis_type MuTect";
    $cmd .= " --reference_sequence " . $self->reference;
    $cmd .= " --input_file:normal " . $self->normal_bam;
    $cmd .= " --input_file:tumor " . $self->tumor_bam;
    $cmd .= " --out " . $self->output_file;
    
    #optional arguments
    $cmd .= " --cosmic " . $self->cosmic_vcf if $self->cosmic_vcf;
    $cmd .= " --dbsnp " . $self->dbsnp_vcf if $self->dbsnp_vcf;
    $cmd .= " --tumor_sample_name " . $self->tumor_sample_name if $self->tumor_sample_name;
    $cmd .= " --normal_sample_name " . $self->normal_sample_name if $self->normal_sample_name;
    $cmd .= " --only_passing_calls" if $self->only_passing_calls;
    $cmd .= " --vcf " . $self->vcf if $self->vcf;
    $cmd .= " --coverage_file " . $self->coverage_file if $self->coverage_file;
    $cmd .= " --intervals " . join(" --intervals ", $self->intervals) if($self->intervals); #do we get back a list or an array ref?


    my @output_files = grep {defined $_} ($self->output_file, $self->vcf, $self->coverage_file);

    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
        output_files => \@output_files,
        skip_if_output_is_present => 0,
        allow_zero_size_output_files => 1,
    );
=cut
java -Xmx5g -jar /gscuser/dlarson/mutect/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence /gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa --input_file:normal /gscmnt/gc13001/info/model_data/2889433023/build136498879/alignments/136268587.bam --input_file:tumor /gscmnt/gc9006/info/model_data/2891225638/build136025755/alignments/135842703.bam --vcf aml31relapse_mutect_test100.vcf --coverage_file coverage100.txt --out aml31relapse_call_stats100.txt --cosmic /gscmnt/sata135/info/medseq/dlarson/aml31_mutect_manual/b37_cosmic_v54_120711.vcf --dbsnp /gscmnt/sata135/info/medseq/dlarson/aml31_mutect_manual/gatk_bundle_v2.3_b37/dbsnp_137.b37.vcf --intervals Y:34482907-59373566 --intervals MT:1-16569 
=cut
}

1;

