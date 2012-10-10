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
    my $vep_command = Genome::Db::Ensembl::Vep->create(
        input_file => "-",
        output_file => $output_file,
        format => "vcf",
        condel => "b",
        polyphen => "b",
        sift => "b",
        hgnc => 1,
    );

    return $vep_command->execute;
}

1;
