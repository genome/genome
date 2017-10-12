package Genome::Db::Ensembl::Gnomad;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::Db::Ensembl::Gnomad {
    is => ['Genome::SoftwareResult::StageableSimple'],
    has_param => [
        version => {
            is => 'Text',
            doc => "Version of gnomad data",
        },
        species => {
            is => 'Text',
            doc => 'Species name for the data',
        },
        exac_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'if provided, also pull in this exac file',
        },
        reference_name => {
            is => 'Text',
            doc => 'Ensembl reference build identifier',
        },
    ],
    has_metric => [
        ftp_path => {
            is => 'Text',
            default_value => 'ftp://ftp.ensembl.org/pub/data_files/',
        },
    ],
};

sub _run {
    my $self = shift;


    my $gnomad_pattern = 'gnomad.%s.%s.sites.%s.noVEP.vcf.gz';

    for my $type ('exomes', 'genomes') {
        my $file = sprintf($gnomad_pattern, $type, $self->version, $self->reference_name);
        my $index = "$file.tbi";

        $self->_stage_ftp_file($file);
        $self->_stage_ftp_file($index);
    }

    if ($self->exac_version) {
        my $exac_pattern = 'ExAC.%s.%s.vcf.gz';

        my $file = sprintf($exac_pattern, $self->exac_version, $self->reference_name);
        my $index = "$file.tbi";

        $self->_stage_ftp_file($file);
        $self->_stage_ftp_file($index);
    }

    return 1;
}

sub _stage_ftp_file {
    my $self = shift;
    my $file = shift;

    my $ftp_dir = join('/',
        $self->ftp_path,
        $self->species,
        $self->reference_name,
        'variation_genotype'
    );

    my $staging_dir = $self->temp_staging_directory;

    my @cmd = ('wget', "$ftp_dir/$file", '-O', "$staging_dir/$file");
    Genome::Sys->shellcmd(
        cmd => \@cmd,
        output_files => ["$staging_dir/$file"]
    );

    return 1;
}

sub resolve_allocation_subdirectory {
    my $self = shift;

    return File::Spec->join('model_data', 'genome-db-ensembl-gnomad', $self->id);
}

sub resolve_allocation_disk_group_name {
    Genome::Config::get('disk_group_references');
}


1;
