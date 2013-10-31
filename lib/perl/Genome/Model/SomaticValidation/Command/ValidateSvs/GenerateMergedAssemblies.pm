package Genome::Model::SomaticValidation::Command::ValidateSvs::GenerateMergedAssemblies;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ValidateSvs::GenerateMergedAssemblies {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Text',
            doc => 'ID of the somatic validation build upon which to run',
            is_output => 1,
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        output_dir => {
            is_output => 1,
            is => 'Text',
            doc => 'Place where the output goes',
            is_calculated => 1,
            calculate_from => 'build',
            calculate => q{ return join("/", $build->data_directory, "validation/sv"); }
        },
        skip => {
            is_output => 1,
            default_value => 0,
            doc => 'If there is no data to validate, signal the remaining steps not to run.',
        },
    ],


    has_optional_transient => [
        sv_call_file => {
            is => 'Text',
            doc => 'The existing SV calls to validate',
        },
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            doc => 'Original data now undergoing validation',
        },
    ],

    doc => 'generates a merged callset and FASTA for validating by alignment',
};

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;
    my $build = $self->build;

    unless($build) {
        die $self->error_message('No build found.');
    }

    Genome::Sys->create_directory($self->output_dir);

    my $sv_file = $self->_resolve_svs_input();
    unless( $sv_file ) {
        #This is okay for "discovery" and "extension" validation runs.
        $self->status_message('Skipping SV validation due to lack of inputs.');
        $self->skip(1);
        return 1;
    }

    #gather up the pieces we need
    my $patient_id = $build->tumor_sample->name;
    my $ref_seq_build = $build->reference_sequence_build;
    my $reference_fasta = $ref_seq_build->full_consensus_path('fa');
    my $tumor_val_bam = $build->tumor_bam;
    my $normal_val_bam = $build->normal_bam;

    unless($normal_val_bam) {
        $self->status_message('Skipping SV validation due to lack of normal BAM.');
        $self->skip(1);
        return 1;
    }

    my ($merged_output_file, $merged_fasta_file) = $self->_generate_merged_callset();
    unless (-s $merged_fasta_file) {
        $self->status_message('Skipping SV validation due to empty merged fasta file.');
        $self->skip(1);
        return 1;
    }

    my $readcount_output = "$merged_output_file.readcounts";
    my $validation_remap_cmd = Genome::Model::Tools::Sv::AssemblyPipeline::RemapReads->create(
        assembly_file => $merged_fasta_file,
        sv_file => $merged_output_file,
        tumor_bam => $tumor_val_bam,
        normal_bam => $normal_val_bam,
        patient_id => "VAL.$patient_id",
        output_file => $readcount_output,
        reference_fasta => $ref_seq_build->full_consensus_path('fa'),
    );
    unless($validation_remap_cmd->execute) {
        die $self->error_message('Failed to run remap-reads on validation data');
    }

    my $classify_cmd = Genome::Model::Tools::Sv::AssemblyPipeline::ClassifyEvents->create(
        readcount_file => $readcount_output
    );
    unless($classify_cmd->execute) {
        die $self->error_message('Failed to classify events');
    }

    my $somatic_sv_file = "$readcount_output.somatic"; # from classify_cmd
    unless (_found_somatic_SVs($somatic_sv_file)) {
        $self->debug_message("Found no somatic SVs in $somatic_sv_file, skipping SV detection.");
        $self->skip(1);
        return 1;
    }

    if($self->somatic_variation_build) {
        my $variation_build = $self->somatic_variation_build;

        if($self->somatic_variation_build->can("tumor_bam")) {
            my $tumor_wgs_bam = $variation_build->tumor_bam;
            my $normal_wgs_bam = $variation_build->normal_bam;

            my $wgs_remap_cmd = Genome::Model::Tools::Sv::AssemblyPipeline::RemapReads->create(
                assembly_file => $merged_fasta_file,
                sv_file => $somatic_sv_file,
                tumor_bam => $tumor_wgs_bam,
                normal_bam => $normal_wgs_bam,
                patient_id => "WGS.$patient_id",
                output_file => "$readcount_output.somatic.wgs_readcounts",
                reference_fasta => $ref_seq_build->full_consensus_path('fa'),
            );
            unless($wgs_remap_cmd->execute) {
                die $self->error_message('Failed to run remap-reads on wgs data');
            }

            $self->_process_wgs_readcounts($wgs_remap_cmd->output_file, $wgs_remap_cmd->patient_id);
        }
    }

    return 1;
}

sub _found_somatic_SVs {
    my ($sv_file) = @_;
    # The header is the first line...
    # an optimization could be to line_count 'head -n 2 <sv_file>'
    # since all we want to know is if it has more than one line.
    if (Genome::Sys->line_count($sv_file) > 1) {
        return 1;
    }
    return 0;
}

sub _resolve_svs_input {
    my $self = shift;
    my $build = $self->build;

    my @sv_input_files;
    if(my $sv_list = $build->sv_variant_list) {
        $self->somatic_variation_build(Genome::Model::Build->get($sv_list->source_build_id));
        my $sv_file = join("/", $sv_list->output_dir, "svs.hq");
        if(-s $sv_file) {
            push @sv_input_files, $sv_file;
        }
    }

    if($build->sv_detection_strategy) {
        push @sv_input_files, $build->data_set_path('variants/svs', undef, 'hq');
    }

    return unless @sv_input_files;

    if(@sv_input_files > 1) {
        my $merged_file = join("/", $self->output_dir, "assembly_input");
        my $merge_cmd = Genome::Model::Tools::Breakdancer::MergeFiles->create(
            input_files => join(',', @sv_input_files),
            output_file => $merged_file,
        );
        unless($merge_cmd->execute) {
            die $self->error_message('Failed to generate merged SV call file');
        }

        $self->sv_call_file($merged_file);
    } else {
        $self->sv_call_file($sv_input_files[0]);
    }

    return 1;
}

sub _generate_merged_callset {
    my $self = shift;
    my $build = $self->build;

    my $assembly_input_file = join("/", $self->output_dir, "assembly_input");
    Genome::Sys->create_symlink($self->sv_call_file, $assembly_input_file) if $self->sv_call_file ne $assembly_input_file;
    my $assembly_output_file = join("/", $self->output_dir, "assembly_output.csv");
    my $assembly_output_fasta = join("/", $self->output_dir, "assembly_output.fasta");
    my $assembly_output_cm = join("/", $self->output_dir, "assembly_output.cm");

    my $tumor_val_bam = $build->tumor_bam;
    my $normal_val_bam = $build->normal_bam;
    my $ref_seq_build = $build->reference_sequence_build;
    my $reference_fasta = $ref_seq_build->full_consensus_path('fa');

    my $assembly_cmd = Genome::Model::Tools::Sv::AssemblyValidation->create(
        bam_files => join(",", $tumor_val_bam, $normal_val_bam),
        output_file => $assembly_output_file,
        sv_file => $assembly_input_file,
        asm_high_coverage => 1,
        min_size_of_confirm_asm_sv => 10,
        breakpoint_seq_file => $assembly_output_fasta,
        cm_aln_file => $assembly_output_cm,
        reference_file => $reference_fasta,
    );
    unless($assembly_cmd->execute) {
        die $self->error_message('Failed to execute assembly-validation command.');
    }

    #merge assembled callsets requires an "index" of the files to use
    my $index_file = $assembly_output_file . '.index';
    Genome::Sys->write_file($index_file, join("\t", "calls", $assembly_output_file, $assembly_output_fasta));

    my $merged_output_file = "$assembly_output_file.merged";
    my $merged_fasta_file = "$assembly_output_fasta.merged";
    my $merge_cmd = Genome::Model::Tools::Sv::MergeAssembledCallsets->create(
        index_file => $index_file,
        output_file => $merged_output_file,
        output_fasta => $merged_fasta_file,
    );
    unless($merge_cmd->execute) {
        die $self->error_message('Failed to merge callsets');
    }

    return ($merged_output_file, $merged_fasta_file);
}

sub _process_wgs_readcounts {
    my $self = shift;
    my $wgs_readcounts_file = shift;
    my $wgs_patient_id = shift;

    my $normal_wgs_reads_cutoff = 0; #Possibly a future PP param?

    my $somatic_file = "$wgs_readcounts_file.somatic";
    my $somatic_fh = Genome::Sys->open_file_for_writing($somatic_file);

    my $readcount_fh = Genome::Sys->open_file_for_reading($wgs_readcounts_file);

    while (my $line = $readcount_fh->getline) {
        if ( $line =~ /^#/ ) { print $somatic_fh $line; next; }
        if ( $line =~ /no\s+fasta\s+sequence/ ) { next; }
        if ( $line =~ /$wgs_patient_id.normal.svReadCount\:(\d+)/i ) {
            my ($normal_sv_readcount) = $line =~ /$wgs_patient_id.normal.svReadCount\:(\d+)/i;
            if ($normal_sv_readcount > $normal_wgs_reads_cutoff) { next; }
            else { print $somatic_fh $line; next; }
        }
    }
    $readcount_fh->close;
    $somatic_fh->close;

    return $somatic_file;
}


1;
