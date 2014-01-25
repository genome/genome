package Genome::Model::PhenotypeCorrelation::Command::VepWrapper;

use Sort::Naturally qw/nsort/;
use File::Basename qw/basename/;
use File::Temp;
use Genome;
use Workflow::Simple;
use POSIX qw/WIFEXITED/;
use Carp qw/confess/;
use File::Path qw/mkpath/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::VepWrapper {
    is => ["Genome::Command::Base"],
    doc => "Run VEP in parallel by chromosome.",
    has_input => [
        input_vcf => {
            is => "File",
            doc => "The vcf file of variants to annotate",
        },
        ensembl_annotation_build_id => {
            is => 'String',
            doc => 'ID of ImportedAnnotation build with the desired ensembl version.',
        },
        region => {
            is => "String",
            doc => "The region to annotate",
        },
        output_file => {
            is => "File",
            doc => "The final merged output file",
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
This command annotates vcf variants in a particular region with Vep.
EOS
}

sub execute {
    my $self = shift;

    my $tabix = Genome::Model::Tools::Tabix->tabix_path;
    my $vcf = $self->input_vcf;
    my $region = $self->region;
    my $output_file = $self->output_file;

    $self->status_message("VEP processing region $region");
    my $tabix_cmd = "$tabix $vcf $region";

    open(STDIN, "-|", $tabix_cmd);
    my $vep_command = Genome::Db::Ensembl::Command::Vep->create(
        input_file => "-",
        output_file => $output_file,
        ensembl_annotation_build_id => $self->ensembl_annotation_build_id,
        format => "vcf",
        condel => "b",
        polyphen => "b",
        sift => "b",
        quiet => 1,
    );

    return $vep_command->execute;
}

1;
