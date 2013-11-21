package Genome::Model::Tools::Samtools::Mpileup;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $VERSION = '0.01';
use Cwd;
use File::Basename;
use File::Path;

class Genome::Model::Tools::Samtools::Mpileup {
    is => 'Command',
    has_optional_input => [
    ref_fasta => {
        is=>'Text',
        is_optional=>0,
    },
    bed_file => {
        is=>'Text',
        is_optional=>0,
    },
    bams => {
        is_many=>1,
        is=>'Text',
        is_optional=>0,
    },
    output_vcf => {
        is=>'Text',
        is_optional=>0,
        is_output=>1,
    },
    ],
    has_param => [
    lsf_resource => {
        is => 'Text',
        default => "-M 5000000 -R 'select[type==LINUX64 && mem>4500 && tmp>20000] rusage[mem=4500]'",
    },
    lsf_queue => {
        is => 'Text',
        default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
    },
    ],

};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}
#samtools mpileup @bams -uf $ref -Dl $bed_file | bcftools view -g - | bgzip -c > $output_vcf_gz
sub execute {
    $DB::single=1;
    my $self=shift;
    my @bams = $self->bams;
    my $ref = $self->ref_fasta;
    my $bed_file = $self->bed_file;
    my $output_vcf_gz = $self->output_vcf;
    my($file, $path, $suffix) = fileparse($output_vcf_gz, ".vcf.gz");
    unless(-d $path) {
        mkpath($path);
    }
    my $cmd = "samtools mpileup @bams -uf $ref -Dl $bed_file | bcftools view -g - | bgzip -c > $output_vcf_gz";
    my $rv = Genome::Sys->shellcmd( cmd => $cmd, input_files => [$bed_file, $ref, @bams]);
    if($rv != 1) {
        return;
    }
    return 1;
}

1;
