package Genome::Model::Tools::VerifyBamId;

use strict;
use warnings;
use Genome;
use Carp "confess";

class Genome::Model::Tools::VerifyBamId {
    is => 'Command::V2',
    has => [
        vcf => {
            is => "File",
            doc => "input VCF file containing individual genotypes or AF or AC/AN fields in the INFO field. gzipped VCF is also allowed",
        },
        bam => {
            is => "File",
            doc => "a BAM (Binary Alignment Map) file of a sequence reads",
        },
        out_prefix => {
            is => "Text",
            doc => "output prefix of output files - [outPrefix].{selfRG,selfSM,bestRG,bestSM,depthRG,depthSM} will be created.",
        },
        max_depth => {
            is => "Integer",
        },
        precise => {
            is => "Boolean",
            doc => "Recommended for max-depth > 20",
        },
        version => {
            is => "Text",
            doc => "Version of VerifyBamId to run",
        },
    ],
};

sub execute {
    my $self = shift;

    my $cmd = $self->_get_cmd;
    return Genome::Sys->shellcmd(cmd => $cmd, input_files => [$self->vcf, $self->bam]); 
}

sub _get_cmd {
    my $self = shift;
    my $executable = _get_exe_path($self->version);
    my $precise = "";
    if ($self->precise) {
        $precise = "--precise";
    }
    return join(" ", $executable, "--vcf", $self->vcf, "--bam", $self->bam, 
                    "--out", $self->out_prefix, "--maxDepth", $self->max_depth,
                    $precise, "--ignoreRG");
}

sub _get_exe_path {
    my $version = shift;
    my $path = "/usr/bin/verifyBamID$version";
    unless (-x $path) {
        confess("version $version of verifyBamID is not available");
        return;
    }
    return $path;
}

1;

