package Genome::Model::Build::MetagenomicCompositionShotgun;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::MetagenomicCompositionShotgun{
    is => 'Genome::Model::Build',
    has =>[
        _final_metagenomic_bam => {
            is_calculated => 1,
            calculate_from => ['data_directory'],
            calculate => sub {
                my ($data_directory) = @_;
                $data_directory."/metagenomic_alignment.combined.sorted.bam";
            },
        },
        _contamination_screen_alignment_build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            via => 'from_build_links',
            to => 'from_build',
            where => [role => 'contamination_screen_alignment_build'],
        },
        contamination_screen_reference=>{
            is => 'Genome::Model::Build::ReferenceAlignment',
            via => 'model',
        },
        _metagenomic_alignment_builds => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            is_many => 1,
            via => 'from_build_links',
            to => 'from_build',
            where => [role => 'metagenomic_alignment_build'],
        },
        metagenomic_references =>{
            is =>'Genome::Model::Build::ReferenceAlignment',
            via => 'model',
        },
        _unaligned_metagenomic_alignment_build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            is_many=> 1,
            via => 'from_build_links',
            to => 'from_build',
            where => [role => 'unaligned_metagenomic_alignment_build'],
        },
        unaligned_metagenomic_alignment_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            via => 'model',
        },
        _first_viral_verification_alignment_build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            is_many=> 1,
            via => 'from_build_links',
            to => 'from_build',
            where => [role => 'first_viral_verification_alignment_build'],
        },
        first_viral_verification_alignment_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            via => 'model',
        },
        _second_viral_verification_alignment_build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            is_many=> 1,
            via => 'from_build_links',
            to => 'from_build',
            where => [role => 'second_viral_verification_alignment_build'],
        },
        second_viral_verification_alignment_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            via => 'model',
        },
        refcov_output => {
            calculate_from => [qw/ reports_directory /],
            calculate => q| return $reports_directory.'/report_combined_refcov_regions_file.regions.txt'; |,
        },
        read_count_file => {
            calculate_from => [qw/ reports_directory /],
            calculate => q| return $reports_directory.'/read_count_output'; |,
        }
    ],
};

sub calculate_estimated_kb_usage {
    return 50_000_000;
}

sub sra_sample_id {
    my $self = shift;
    my @id = $self->instrument_data;
    die "no instrument data unless id" unless @id;
    my $sra_sample_id = $id[0]->sra_sample_id;
    die "no sra_sample_id specified for instrument data" unless $sra_sample_id;
    return $sra_sample_id;
}

# META REFS
sub metagenomic_reference_hmp_dir {
    my $self = shift;

    my @metagenomic_references = $self->model->metagenomic_references;
    my ($hmp_dir) = grep { -d $_ } map { $_->data_directory.'/hmp' } @metagenomic_references;
    if ( not $hmp_dir ) {
        $self->error_message('No hmp directory found in reference builds: '.join(' ', map { $_->__display_name__ } @metagenomic_references));
        return;
    }

    return $hmp_dir;
}

sub metagenomic_reference_regions_file {
    my $self = shift;

    my $hmp_dir = $self->metagenomic_reference_hmp_dir;
    return if not $hmp_dir;

    return $hmp_dir.'/combined_refcov_regions_file.regions.bed';
}

sub metagenomic_reference_taxonomy_file {
    my $self = shift;

    my $hmp_dir = $self->metagenomic_reference_hmp_dir;
    return if not $hmp_dir;

    return $hmp_dir.'/Bact_Arch_Euky.taxonomy.txt';
}

sub metagenomic_reference_viral_headers_file {
    my $self = shift;

    my $hmp_dir = $self->metagenomic_reference_hmp_dir;
    return if not $hmp_dir;

    return $hmp_dir.'/viruses_nuc.fasta.headers';
}

sub metagenomic_reference_viral_taxonomy_file {
    my $self = shift;

    my $hmp_dir = $self->metagenomic_reference_hmp_dir;
    return if not $hmp_dir;

    return $hmp_dir.'/viruses_taxonomy_feb_25_2010.txt';
}

# BUILD DIFF
sub files_ignored_by_diff {
    return qw(
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded/report.xml
        build.xml
    );
}

sub dirs_ignored_by_diff {
    return qw(
        logs/
    );
}

1;

