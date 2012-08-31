package Genome::Model::RnaSeq::DetectFusionsResult::BreakFusionResult;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::BreakFusionResult{
    is => "Genome::Model::RnaSeq::DetectFusionsResult",
    has => [
        detector_params => {
            #TODO
            doc => 'params to use for sub steps in this detector breakdancer params : tigra-sv params : break annot params'
        }
    ]
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_output_directory();
    $self->_prepare_staging_directory();

    my ($breakdancer_params, $tigra_params, $break_annot_params) = split(/:/, $self->detector_params);
    my $breakdancer_cfg_file = $self->temp_staging_directory . "/breakdancer.cfg";
    my $cfg_cmd = $self->_path_for_command($self->version, "bam2cfg.RNAseq.pl");
    $cfg_cmd .= " " . $self->alignment_result->bam_file . " > " . $breakdancer_cfg_file;

    Genome::Sys->shellcmd(
        cmd => $cfg_cmd,
        input_files => [$self->alignment_result->bam_file],
        output_files => [$breakdancer_cfg_file],
    );

    my $breakdancer_output = $self->temp_staging_directory . "/breakdancer";

    my $breakdancer_cmd = Genome::Model::Tools::DetectVariants2::Breakdancer->create(
        config_file => $breakdancer_cfg_file,
        version => "1.2",
        params => ":" . $breakdancer_params,
        reference_build_id => $self->alignment_result->reference_build->id,
        output_directory => $breakdancer_output,
        aligned_reads_input => $self->alignment_result->bam_file,
        aligned_reads_sample => ($self->alignment_result->instrument_data)[0]->sample->name,
    );

    unless($breakdancer_cmd->execute()){
        die($self->error_message("Running breakdancer failed!"));
    }

    my $index_result = Genome::Model::RnaSeq::DetectFusionsResult::BreakFusionResult::Index->get_or_create(
        reference_build => $self->alignment_result->reference_build
    );

    my $break_annot_output_file =  $self->temp_staging_directory . "/svs.hq.anno";
    $self->_run_break_annot_cmd( $breakdancer_output . "/svs.hq" , $break_annot_output_file, $break_annot_params, $index_result);

    my $alignment_result_fasta =  $self->alignment_result->reference_build->full_consensus_path("fa");
    my $tigra_output_file = $self->temp_staging_directory . "/svs.hq.tigra-sv.fa";

    #my $tigra_cmd = $self->_path_for_command($self->version,"tigra-sv.static.linux-x86_64");
    my $tigra_cmd ="tigra-sv.static.linux-x86_64";
    $tigra_cmd .= " -R " . $alignment_result_fasta;
    $tigra_cmd .= " -b -o " . $tigra_output_file . " ";
    $tigra_cmd .= $break_annot_output_file . " " . $self->alignment_result->bam_file;

    Genome::Sys->shellcmd(
        cmd => $tigra_cmd,
        output_files => [$tigra_output_file],
        input_files => [$self->alignment_result->bam_file, $alignment_result_fasta, $break_annot_output_file ]
    );

    my $blat_sv_output_file = $tigra_output_file . ".BLAT.csv";

    my $blat_sv_contig_cmd = $self->_path_for_command($self->version, "BlatSVContig.pl");
    $blat_sv_contig_cmd .= " -r " . $index_result->twobit_file .  " $tigra_output_file > $blat_sv_output_file";

    Genome::Sys->shellcmd(
        cmd => $blat_sv_contig_cmd,
        output_files => [$blat_sv_output_file],
        input_files => [$tigra_output_file, $index_result->twobit_file]
    );

    my $break_annot_blat_output_file = $blat_sv_output_file . "anno";
    $self->_run_break_annot_cmd($blat_sv_output_file, $break_annot_blat_output_file, $break_annot_params, $index_result);

    $self->_promote_data();
    $self->_remove_staging_directory();
    $self->_reallocate_disk_allocation();

    return $self;
}

sub _staging_disk_usage {
    return 60 * 1024 * 1024;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/break-fusion/' . $self->id;
}

sub _run_break_annot_cmd {
    my ($self, $input_file , $output_file, $break_annot_params, $index_result) = @_;

    my $break_annot_cmd = $self->_path_for_command($self->version, "BreakAnnot.pl");
    $break_annot_cmd .= $break_annot_params . " -g "  . $index_result->output_dir ."/refGene.txt -C " . $index_result->output_dir . "/chainSelf.tab ";
    $break_annot_cmd .= $input_file . " > " . $output_file;

    return  Genome::Sys->shellcmd(
        cmd => $break_annot_cmd,
        input_files => [$index_result->output_dir . "/refGene.txt", $input_file],
        output_files => [$output_file]
    );
}

1;
