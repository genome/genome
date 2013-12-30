package Genome::Model::Tools::Tcga::CreateSubmissionArchive;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;

my $IDF_FILE_EXTENSION = "idf";
my $SDRF_FILE_EXTENSION = "sdrf";
my $MANIFEST_FILE_NAME = "MANIFEST.txt";

class Genome::Model::Tools::Tcga::CreateSubmissionArchive {
    is => 'Command::V2',
    has => [
        models => {
            doc => "Models to include in the archive",
            is => 'Genome::Model::SomaticVariation',
            is_many => 1,
        },       
        output_dir => {
            is => "Text",
            doc => "Directory where the archives should be written",
        },
        archive_name => {
            is => "Text",
            doc => "Name of the archive, e.g genome.wustl.edu_SARC.IlluminaHiSeq_DNASeq",
        },
        archive_version => {
            is => "Text",
            doc => "Version of the archive, e.g. 1.0.0",
        },
        create_archive => {
            is => "Boolean",
            default_value => 0,
            doc => "Create the final tar.gz and md5 files",
        },
        cghub_id_file => {
            is => "Text",
            doc => "A tab-delimited file that maps bam files (column header: BAM_path) to CGHub ids (column header: CGHub_ID).",
        },
    ],
    has_optional => [
        somatic_maf_file => {
            is => "Text",
            doc => "Somatic maf file to be included in the archive",
        },
        protected_maf_file => {
            is => "Text",
            doc => "Germline maf file to be included in the archive",
        },
    ],
};

sub execute {
    my $self = shift;
    my @sdrf_rows;
    my $idf = Genome::Model::Tools::Tcga::Idf->create;
    
    for my $somatic_model ($self->models) {
        $idf->add_pp_protocols($somatic_model->last_succeeded_build->processing_profile);
    }

    my $sdrf = Genome::Model::Tools::Tcga::Sdrf->create(idf => $idf, 
                                                        cghub_id_file => $self->cghub_id_file,
                                                        archive_name => $self->complete_archive_name("Level_2"));

    my $vcf_archive_dir = $self->output_dir."/".$self->complete_archive_name("Level_2");
    Genome::Sys->create_directory($vcf_archive_dir);
    my $magetab_archive_dir = $self->output_dir."/".$self->complete_archive_name("mage-tab");
    Genome::Sys->create_directory($magetab_archive_dir);

    my %patient_ids;
    for my $somatic_model ($self->models) {
        my $somatic_build = $somatic_model->last_succeeded_build;
        unless($somatic_build) {
            $self->error_message("Could not resolve build from model ".$somatic_model->__display_name__);
            return;
        }
        my $normal_build = $somatic_build->normal_build;
        my $tumor_build = $somatic_build->tumor_build;

        my $patient_id = $self->resolve_patient_id($somatic_build);
        my $patient_id_counter = ++$patient_ids{$patient_id};
        
        my $snvs_vcf = $self->construct_vcf_name("snv", $patient_id, $patient_id_counter);
        my $indels_vcf = $self->construct_vcf_name("indel", $patient_id, $patient_id_counter);

        for my $variant_type (qw(snv indel)) {
            my $local_file = $somatic_build->data_directory."/variants/".$variant_type."s_tcga/".$variant_type."s_tcga.vcf";
            die "Tcga compliant $variant_type vcf not found for build ".$somatic_build->id unless(-s $local_file);
            Genome::Sys->copy_file($local_file, "$vcf_archive_dir/".$self->construct_vcf_name($variant_type, $patient_id, $patient_id_counter));
        }

        my $vcf_sample_info = $self->get_sample_info_from_vcf("$vcf_archive_dir/$snvs_vcf");
        unless (defined $vcf_sample_info) {
            die $self->error_message("No SAMPLE vcf header in $snvs_vcf for build ".$somatic_build->id);
        }

        for my $build(($normal_build, $tumor_build)) {
            my $sample_info = $self->get_info_for_sample($build->subject->extraction_label, $vcf_sample_info);
            
            for my $vcf($snvs_vcf, $indels_vcf) {
                push @sdrf_rows, $sdrf->create_vcf_row($build, $somatic_build, $vcf, $sample_info);
            }

            for my $maf_type (qw(somatic protected)) {
                my $maf_accessor = $maf_type."_maf_file";
                if ($self->$maf_accessor) {
                    my $maf_name = $self->complete_archive_name.".$maf_type.maf";
                    unless(-s $vcf_archive_dir."/".$maf_name) {
                        Genome::Sys->copy_file($self->$maf_accessor, $vcf_archive_dir."/".$maf_name);
                    }
                    push @sdrf_rows, $sdrf->create_maf_row($build, $somatic_build, $maf_name, $sample_info);
                }
            }
        }
    }

    my $sdrf_name = $self->complete_archive_name.".".$SDRF_FILE_EXTENSION.".txt";
    my $idf_name = $magetab_archive_dir."/".$self->complete_archive_name.".".$IDF_FILE_EXTENSION.".txt";
    $idf->sdrf_file($sdrf_name);
    $idf->print_idf($idf_name);
    $sdrf->print_sdrf($magetab_archive_dir."/".$sdrf_name, @sdrf_rows);

    $self->print_manifest($vcf_archive_dir);
    $self->print_manifest($magetab_archive_dir);

    if ($self->create_archive) {
        $self->tar_and_md5_dir($vcf_archive_dir);
        $self->tar_and_md5_dir($magetab_archive_dir);
    }

    return 1;
}

sub complete_archive_name {
    my $self = shift;
    my $infix = shift;

    my @parts;
    push @parts, $self->archive_name;
    if (defined $infix) {push @parts, $infix;}
    push @parts, $self->archive_version;
    return join(".", @parts);
}

sub construct_vcf_name {
    my $self = shift;
    my $variant_type = shift;
    my $patient_id = shift;
    my $patient_id_counter = shift;
    my $snvs_vcf = "genome.wustl.edu.$patient_id.$variant_type.".$patient_id_counter.".vcf";
}

sub print_manifest {
    my $self = shift;
    my $directory = shift;
    my $cd_cmd = "cd $directory ; md5sum * > $MANIFEST_FILE_NAME";
    Genome::Sys->shellcmd(cmd => $cd_cmd, output_files => [$directory."/$MANIFEST_FILE_NAME"]);
}

sub get_info_for_sample {
    my $self = shift;
    my $desired_sample = shift;
    my $sample_info_collection = shift;
    for my $sample (@$sample_info_collection) {
        if ($sample->{"ID"}->{content} eq $desired_sample) {
            return $sample;
        }
    }
    die "Info for sample $desired_sample was not available";
}

sub get_sample_info_from_vcf {
    my $self = shift;
    my $vcf_file = shift;
    my $vcf_reader = new Genome::File::Vcf::Reader($vcf_file);
    return $vcf_reader->header->metainfo->{"SAMPLE"};
}

sub tar_and_md5_dir {
    my $self = shift;
    my $dir = shift;

    my $tar_file = "$dir.tar.gz";
    my $tar_cmd = "tar -czf $tar_file $dir";
    Genome::Sys->shellcmd(cmd => $tar_cmd, output_files => [$tar_file]);
    my $md5 = Genome::Sys->md5sum($tar_file);
    my $md5_cmd = "echo $md5 > $tar_file.md5";
    Genome::Sys->shellcmd(cmd => $md5_cmd, output_files => ["$tar_file.md5"]);
    return 1;
}

sub resolve_patient_id {
    my $self = shift;
    my $build = shift;
    my $patient_id = $build->subject->source->upn;
    die "Could not resolve patient_id for build ".$build->id unless (defined $patient_id);
    return $patient_id;
}




1;

