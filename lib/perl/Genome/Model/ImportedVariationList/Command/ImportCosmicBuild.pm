package Genome::Model::ImportedVariationList::Command::ImportCosmicBuild;

use warnings;
use strict;

use Genome;

class Genome::Model::ImportedVariationList::Command::ImportCosmicBuild {
   is => 'Command::V2',
   has_input => [
       vcf_file_urls => {
           is => 'Text',
           is_many => 1,
           doc => 'Paths to the VCF file(s) on the Cosmic ftp server'
       },
       version => {
           is => 'Text',
           doc => 'The version of the build to create',
       },
       reference_sequence_build => {
           is => 'Genome::Model::Build::ReferenceSequence',
           doc => 'The reference upon which the Cosmic build will be based',
       },
    ],
    has_output => [
        build => {
            is_optional => 1,
            is => 'Genome::Model::Build::ImportedVariation',
        },
    ],
};

sub kilobytes_requested {
    my $self = shift;
    return 20_971_520;
}

sub execute {
    my $self = shift;
    my $allocation = Genome::Disk::Allocation->create(
        kilobytes_requested => $self->kilobytes_requested, 
        disk_group_name => 'info_genome_models', 
        allocation_path => 'build_merged_alignments/import_cosmic_' . $self->version . '_' . Genome::Sys->md5sum_data(join(",",$self->vcf_file_urls)),
        owner_class_name => 'Genome::Model::ImportedVariationList::Command::ImportCosmicBuild',
        owner_id => $self->id
    );

    local $ENV{'TMPDIR'} = $allocation->absolute_path;

    my $original_file_path = $allocation->absolute_path . '/cosmic.vcf';
    my @urls = $self->vcf_file_urls;
    my $import_vcf = Genome::Db::Cosmic::Command::ImportVcf->create(
        urls => \@urls,
        output_file => $original_file_path,
    );

    unless ($import_vcf->execute()){
        die($self->error_message("VCF file download and merge failed"));
    }

    my $import_cmd = Genome::Model::ImportedVariationList::Command::ImportVariants->create(
        input_path => $original_file_path,
        reference_sequence_build => $self->reference_sequence_build,
        source_name => "cosmic",
        description => 'Imported VCF file from Cosmic ' . join(",",$self->vcf_file_urls),
        variant_type => "snv",
        format => "vcf",
        version => $self->version,
    );

    my $rv = $import_cmd->execute;
    
    $self->build($import_cmd->build);
    $allocation->delete;

    return $rv;
}

1;

