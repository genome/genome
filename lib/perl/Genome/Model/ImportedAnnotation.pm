package Genome::Model::ImportedAnnotation;

use strict;
use warnings;

use Data::Dumper;
use Genome;

class Genome::Model::ImportedAnnotation{
    is => 'Genome::ModelDeprecated',
    has_param => [
        annotation_source => {
            is_optional => 0,
            doc => 'Where the annotation comes from (ensembl, genbank, etc.) This value is "combined-annotation" for a combined-annotation model',
        },
        interpro_version => {
            is_optional => 0,
            default_value => 4.5,
            doc => 'Version of interpro used to import interpro results', 
        },
        rna_seq_only => {
            default_value => 0,
            doc => 'True if only a minimal set of annotation files for RNA-Seq pipeline are needed',
        },
    ],
    has =>[
        species_name => {
            is => 'UR::Value',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'species_name', value_class_name => 'UR::Value'],
            is_mutable => 1,
        },
        annotation_import_version => {
            is => 'UR::Value',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'annotation_import_version', value_class_name => 'UR::Value'],
            is_mutable => 1,
        },
        reference_sequence_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'reference_sequence', value_class_name => 'Genome::Model::ImportedReferenceSequence' ],
            is_many => 0,
            is_optional => 1, # TODO: make this non-optional when all data is updated
            is_mutable => 1,
            doc => 'id of the reference sequence model associated with this annotation model',
        },
        reference_sequence => {
            is => 'Genome::Model::ImportedReferenceSequence',
            id_by => 'reference_sequence_id',
        },

    ],
};

sub build_by_version {
    my $self = shift;
    my $version = shift;

    my @builds = grep{ $_->status eq 'Succeeded' || $_->status eq 'Running' } $self->builds;
    @builds =  grep { $_->version eq $version } @builds;
    if (@builds > 1) {
        my $versions_string = join("\n", map { "model_id ".$_->model_id." build_id ".$_->build_id." version ".$_->version } @builds);
        $self->error_message("Multiple builds for version $version of model " . $self->genome_model_id.", ".$self->name."\n".$versions_string."\n");
        die;
    }
    return $builds[0];
}

sub annotation_data_directory {
    my $self = shift;
    my $build = $self->last_complete_build;
    return $build->determine_data_directory;
}

sub annotation_build_for_reference {
    my ($class, $reference) = @_;
    my $build;

    #TODO: Remove this hardcoded crap and come up with an intelligent heuristic

    if($reference->name eq 'NCBI-human-build36'){
        $build = Genome::Model::Build::ImportedAnnotation->get(113049382);
    }
    elsif($reference->name eq 'GRCh37-lite-build37' || $reference->name eq 'g1k-human-build37'){
        $build = Genome::Model::Build::ImportedAnnotation->get(124434505);
    }
    elsif($reference->name eq 'UCSC-mouse-buildmm9'){
        $build = Genome::Model::Build::ImportedAnnotation->get(106410073);
    }
    return $build;
}

sub _resolve_resource_requirements_for_build {
    return "-R 'rusage[mem=32000]' -M 32000000";
}

sub _execute_build {
    my $self = shift;
    my $build = shift;

    my $source = $self->annotation_source;
    unless (defined $source){
        $self->error_message("Could not get imported annotation source!");
        return;
    }

    my $version = $build->version;
    unless (defined $version){
        $self->error_message("Could not get build version!");
        return;
    }

    my $data_directory = $build->data_directory;
    unless (defined $data_directory){
        $self->error_message("Could not get data directory for build!");
        return;
    }
    unless (-d $data_directory){
        Genome::Sys->create_directory($build->data_directory);
        unless (-d $data_directory) {
            $self->error_message("Failed to create new build dir: " . $build->data_directory);
            return;
        }
    }

    unless (-d $build->_annotation_data_directory) {
        Genome::Sys->create_directory($build->_annotation_data_directory);
        unless (-d $build->_annotation_data_directory) {
            $self->error_message("Failed to create annotation directory: ".$build->_annotation_data_directory);
        }
    }

    my $species_name = $build->species_name;
    unless (defined $species_name){
        $self->error_message('Could not get species name!');
        return;
    }

    unless ($self->rna_seq_only) {
        my $name = ucfirst(lc($source));
        my $importer_class_name = join('::', 'Genome', 'Db', $name, 'Command', 'Import', 'Run');
        my $cmd = $importer_class_name->execute(
            data_set => 'Core', 
            imported_annotation_build => $build,
            software_version => $build->annotation_import_version,
        );

        my $tiering_cmd;
        my $annotation_directory = $build->_annotation_data_directory;
        my $bitmasks_directory = $annotation_directory."/tiering_bitmasks";
        unless ( -d $bitmasks_directory) {
            Genome::Sys->create_directory($bitmasks_directory);
            unless (-d $bitmasks_directory) {
                $self->error_message("Failed to create new build dir: " . $bitmasks_directory);
                return;
            }
        }
        my $bed_directory = $annotation_directory."/tiering_bed_files_v3";
        unless ( -d $bed_directory) {
            Genome::Sys->create_directory($bed_directory);
            unless (-d $bed_directory) {
                $self->error_message("Failed to create new build dir: " . $bed_directory);
                return;
            }
        }

        if ($species_name eq 'human' or $species_name eq 'mouse') {
            $tiering_cmd = Genome::Model::Tools::FastTier::MakeTierBitmasks->create(
                output_directory => $annotation_directory."/tiering_bitmasks",
                reference_sequence_build => $build->reference_sequence,
                transcript_version => $build->ensembl_version,
                ucsc_directory => $build->reference_sequence->get_or_create_ucsc_tiering_directory,
                species => $species_name,
                annotation_import_version => $build->annotation_import_version,
            );

            $tiering_cmd->execute;
            foreach my $file (glob $annotation_directory."/tiering_bitmasks/*") {
                my $bed_name = $file;
                $bed_name =~ s/tiering_bitmasks/tiering_bed_files_v3/;
                $bed_name =~ s/bitmask/bed/;
                my $convert_cmd = Genome::Model::Tools::FastTier::BitmaskToBed->create(
                    output_file => $bed_name,
                    bitmask => $file,
                );
                $convert_cmd->execute;
            }
        }

        my $ucsc_directory = $annotation_directory."/ucsc_conservation";
        my $original_ucsc_dir = $build->reference_sequence->get_or_create_ucsc_conservation_directory;
        if ($original_ucsc_dir) {
        Genome::Sys->create_symlink($original_ucsc_dir, $ucsc_directory); 
        }

        #generate the rna seq files
        $self->generate_rna_seq_files($build);

        #Make ROI FeatureList
        $build->get_or_create_roi_bed;

        my $gap_list_command =
            Genome::Model::ImportedAnnotation::Command::FetchUcscGapList->execute(
                annotation_build => $build,
            );
        $build->gap_feature_list($gap_list_command->gap_feature_list);
    }

    return 1;
}

sub get_ensembl_info {
    my $self = shift;
    my $version = shift;
    my ($eversion,$ncbiversion) = split(/_/,$version);

    my $host = defined $ENV{GENOME_DB_ENSEMBL_HOST} ? $ENV{GENOME_DB_ENSEMBL_HOST} : 'mysql1';
    my $user = defined $ENV{GENOME_DB_ENSEMBL_USER} ? $ENV{GENOME_DB_ENSEMBL_USER} : 'mse'; 
    my $password = defined $ENV{GENOME_DB_ENSEMBL_PASSWORD} ? $ENV{GENOME_DB_ENSEMBL_PASSWORD} : undef;

    return ($host, $user, $password);
}

sub generate_rna_seq_files {
    my $self = shift;
    my $build = shift;

    my $cmd = Genome::Model::ImportedAnnotation::Command::CopyRibosomalGeneNames->create(output_file => join('/', $build->_annotation_data_directory, 'RibosomalGeneNames.txt'), species_name => $build->species_name);
    unless($cmd->execute){
        $self->error_message("Failed to generate the ribosomal gene name file!");
        return;
    }

    unless($build->generate_RNA_annotation_files('gtf', $build->reference_sequence_id)){
        $self->error_message("Failed to generate RNA Seq gtf files!");
        return;
    }

    unless($build->generate_RNA_annotation_files('bed', $build->reference_sequence_id,0)){
        $self->error_message("Failed to generate RNA Seq bed files!");
        return;
    }

    unless($build->generate_RNA_annotation_files('bed', $build->reference_sequence_id,1)){
        $self->error_message("Failed to generate RNA Seq squashed bed files!");
        return;
    }

    return 1;
}

sub calculate_snapshot_date {
    my ($self, $genbank_file) = @_;
    my $output = `ls -l $genbank_file`;
    my @parts = split(" ", $output);
    my $date = $parts[5];
    return $date;
}

1;

