package Genome::Model::ImportedVariationList::Command::ImportDbsnpBuild;

use warnings;
use strict;

use Genome;

class Genome::Model::ImportedVariationList::Command::ImportDbsnpBuild {
   is => 'Command::V2',
   has_input => [
       vcf_file_url => {
           is => 'Text',
           doc => 'Path to the full VCF file on the DBSnp ftp server'
       },
       flat_file_pattern => {
           is => 'Text',
           is_optional => 1,
           doc => 'String representing the pattern that the flat file filenames follow with [X] substituted in for the chromosome number',
       },
       vcf_file_pattern => {
           is => 'Text',
           is_optional => 1,
           doc => 'String representing the pattern that the vcf filenames follow with [X] substituted for the chromosome number',
       },
       version => {
           is => 'Text',
           doc => 'The version of the build to create',
       },
       reference_sequence_build => {
           is => 'Genome::Model::Build::ReferenceSequence',
           doc => 'The reference upon which the DBSnp build will be based'
       },
       contig_names_translation_file => {
           is => 'Path',
           is_optional => 1,
           doc => 'File path that contains translations of contig names',
       },
       from_names_column => {
           is => 'Number',
           is_optional => 1,
           doc => '0-based column number containing names you want to translate from',
       },
       to_names_column => {
           is => 'Number',
           is_optional => 1,
           doc => '0-based column number containing names you want to translate to',
       },
       chromosome_names => {
           is => 'String',
           is_many => 1,
           is_optional => 1,
           default_value => ["1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "MT", "X", "Y", "Multi"],
       },
       reference_coordinates => {
           is => 'String',
           is_optional => 1,
           doc => 'reference_coordinates whose coordinates will be used, regex syntax accepted for matching multiple   patch levels',
       },
       import_vcf => {
           is => 'Boolean',
           default => 1,
           doc => "Whether to import the dbsnp build as a vcf",
       },
       import_bed => {
           is => 'Boolean',
           default => 1,
           doc => "Whether to import the dbsnp build as a bed",
       },
   ],
   has_transient_optional_output => [
       build => {
           is => 'Genome::Model::Build::ImportedVariationList',
           doc => 'Build created by this command'
       },
   ],
   has => [
       flat_url => {
           is_optional => 1,
           is_calculated => 1,
           calculate => sub {
               my @a = split('/', $_[0]);
               my $levels;
               if ($_[1]) {
                $levels = 2;
               }
               else {
                $levels = 3;
               }
               join('/', @a[0..(@a-$levels)], 'ASN1_flat');
           },
          calculate_from => ['vcf_file_url', 'vcf_file_pattern'],
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
        disk_group_name => $ENV{GENOME_DISK_GROUP_MODELS}, 
        allocation_path => 'build_merged_alignments/import_dbsnp_' . $self->version . '_' . Genome::Sys->md5sum_data($self->vcf_file_url),
        owner_class_name => 'Genome::Model::ImportedVariationList::Command::ImportDbsnpBuild',
        owner_id => $self->id
    );

    local $ENV{'TMPDIR'} = $allocation->absolute_path;

    my $original_file_path = $allocation->absolute_path . '/merged_dbsnp.vcf';
    if ($self->import_vcf) {
        my $import_vcf = Genome::Model::Tools::Dbsnp::ImportVcf->create(
            vcf_file_url => $self->vcf_file_url,
            ($self->flat_file_pattern ? (flat_file_pattern => $self->flat_file_pattern) : ()),
            ($self->vcf_file_pattern ? (vcf_file_pattern => $self->vcf_file_pattern) : ()),
            ($self->chromosome_names ? (chromosome_names => [$self->chromosome_names]) : ()),
            output_file_path => $original_file_path,
            flat_url => $self->flat_url,
        );

        unless ($import_vcf->execute()){
            die($self->error_message("VCF file download and merge failed"));
        }
    }

    my $bed_file_path;
    if ($self->import_bed) {
        $bed_file_path = $allocation->absolute_path."/dbsnp.bed";
        my $import_bed = Genome::Model::Tools::Dbsnp::Import->create(
            flat_file_url => $self->flat_url,
            ($self->flat_file_pattern ? (filename_pattern => $self->flat_file_pattern) : ()),
            output_file => $bed_file_path,
            ($self->contig_names_translation_file ? (contig_name_translation_file => $self->contig_names_translation_file):()),
            ($self->from_names_column ? (from_names_column => $self->from_names_column):()),
            ($self->to_names_column ? (to_names_column => $self->to_names_column):()),
            ($self->chromosome_names ? (chromosome_names => [$self->chromosome_names]):()),
            ($self->reference_coordinates ? (reference_coordinates => $self->reference_coordinates) : ()),
        );

        unless ($import_bed->execute()){
            die($self->error_message("Bed file import failed"));
        }
    }

    my $import_cmd;
    if ($self->import_vcf and $self->import_bed) {
        my %params = (
            input_path => $original_file_path,
            reference_sequence_build => $self->reference_sequence_build,
            source_name => "dbsnp",
            description => "this had better work!",
            description => 'Imported VCF file from DBSnp ' . $self->vcf_file_url,
            variant_type => "snv",
            format => "vcf",
            version => $self->version,
            bed_file => $bed_file_path,
        );
        
        $import_cmd = Genome::Model::ImportedVariationList::Command::ImportVariants->create(%params);
    }
    elsif ($self->import_vcf) {
        $import_cmd = Genome::Model::ImportedVariationList::Command::ImportVariants->create(
            input_path => $original_file_path,
            reference_sequence_build => $self->reference_sequence_build,
            source_name => "dbsnp",
            description => "this had better work!",
            description => 'Imported VCF file from DBSnp ' . $self->vcf_file_url,
            variant_type => "snv",
            format => "vcf",
            version => $self->version,
        )
    }
    elsif ($self->import_bed) {
        $import_cmd = Genome::Model::ImportedVariationList::Command::ImportVariants->create(
            input_path => $bed_file_path,
            reference_sequence_build => $self->reference_sequence_build,
            source_name => "dbsnp",
            description => "this had better work!",
            description => 'Imported bed file from DBSnp ' . $self->flat_url,
            variant_type => "snv",
            format => "bed",
            version => $self->version,
        )
    }


    my $rv = $import_cmd->execute;
    $self->build($import_cmd->build);
    
    $allocation->delete;
    return $rv;
}

1;

