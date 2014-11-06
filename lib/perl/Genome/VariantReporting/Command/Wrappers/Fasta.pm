package Genome::VariantReporting::Command::Wrappers::Fasta;

use strict;
use warnings;
use Genome;
use File::Basename qw(dirname);

class Genome::VariantReporting::Command::Wrappers::Fasta {
    is => 'Command::V2',
    has_input => [
        model => {
            is => "Genome::Model::SomaticVariation",
        },
        output_directory => {
            is => 'Path',
        },
    ],
    has_calculated_optional => [
        build => {
            calculate_from => [qw(model)],
            calculate => q/$model->last_succeeded_build/,
        },
        translations_file => {
            calculate_from => [qw(output_directory)],
            calculate => q/File::Spec->join($output_directory, "resources.yaml")/,
        },
        combined_fasta_file => {
            calculate_from => [qw(output_directory)],
            calculate => q/File::Spec->join($output_directory, "all.fa")/,
        },
        final_fimo_output_file => {
            calculate_from => [qw(output_directory)],
            calculate => q/File::Spec->join($output_directory, "combined_fimo_output.txt")/,
        },
        meme_directory => {
            calculate_from => [qw(output_directory)],
            calculate => q/File::Spec->join($output_directory, "fimo_output")/,
        },
        fimo_output_file => {
            calculate_from => [qw(meme_directory)],
            calculate => q/File::Spec->join($meme_directory, "fimo.txt")/,
        },
    ],
};

sub execute {
    my $self = shift;
    $self->generate_translations_file;
    $self->run_reports;
    $self->combine_snvs_and_indels;
    $self->run_fimo;
    $self->combine_fimo_output;
}

sub generate_translations_file {
    my $self = shift;
    my %translations;
    my @aligned_bams;
    push @aligned_bams, $self->build->tumor_build->merged_alignment_result->id;
    $translations{aligned_bam_result_id} = \@aligned_bams;

    $translations{reference_fasta} = $self->build->tumor_build->reference_sequence_build->full_consensus_path("fa");
    my %feature_list_ids;
    $translations{feature_list_ids} = \%feature_list_ids;
    $translations{tumor_sample_name} = $translations{tumor};
    $translations{tumor} = $self->build->tumor_build->subject->name;
    YAML::DumpFile(File::Spec->join($self->translations_file), \%translations);
}

sub run_reports {
    my $self = shift; 
    for my $variant_type (qw(snvs indels)) {
        Genome::Sys->create_directory($self->log_directory($variant_type));
        Genome::Sys->create_directory($self->report_directory($variant_type));
        my %params = (
            input_vcf => $self->input_vcf($variant_type),
            variant_type => $variant_type,
            output_directory => $self->report_directory($variant_type),
            plan_file => $self->plan_file($variant_type),
            translations_file => $self->translations_file,
            log_directory => $self->log_directory($variant_type),
        );
        Genome::VariantReporting::Command::CreateReport->execute(%params);
    }
}

sub run_fimo {
    my $self = shift;
    Genome::Sys->create_directory($self->meme_directory);
    my $cmd = "/gscuser/aregier/scratch/meme/meme_4.10.0/meme/bin/fimo --text -oc ".
        $self->meme_directory." /gscuser/aregier/scratch/tfbs/motifs.meme ".$self->combined_fasta_file." > ".$self->fimo_output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$self->fimo_output_file],
    );
}

sub combine_fimo_output {
    my $self = shift;
    my $fimo_file = $self->fimo_output_file;
    my $fimo = Genome::Sys->open_file_for_reading($fimo_file);
    my ($wt, $wt_file) = Genome::Sys->create_temp_file;
    my ($mut, $mut_file) = Genome::Sys->create_temp_file;
    my $header = <$fimo>;
    while (my $line = <$fimo>) {
        chomp $line;
        my @fields = split("\t", $line);
        if ($fields[1] =~ /-wt$/) {
            $fields[1] =~ s/-wt$//;
            print $wt join("\t", @fields)."\n";
        }
        elsif ($fields[1] =~ /-mut$/) {
            $fields[1] =~ s/-mut$//;
            print $mut join("\t", @fields)."\n";
        }
    }
    $wt->close;
    $mut->close;
    my $wt_out = _run_awk_cmd($wt_file);
    my $mut_out = _run_awk_cmd($mut_file);
    my $echo_cmd = 'echo "pattern-name sequence-name start stop strand wt-score mut-score wt-pvalue mut-pvalue wt-sequence mut-sequence" > '.
        $self->final_fimo_output_file;
    Genome::Sys->shellcmd(
        cmd => $echo_cmd,
        output_files => [$self->final_fimo_output_file],
    );
    my $cmd = "join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,2.7,1.8,2.8,1.9,2.9 $wt_out $mut_out >> ".$self->final_fimo_output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$wt_out, $mut_out, $self->final_fimo_output_file],
    );
}

sub _run_awk_cmd {
    my $in = shift;
    my $out = Genome::Sys->create_temp_file_path;
    my $cmd = qq(
        awk '{printf("%s:%s:%s:%s:%s %s %s %s %s %s %s %s %s %s\\n", \$1, \$2, \$3, \$4, \$5, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);}' $in | sort > $out);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$in],
        output_files => [$out],
    );
    return $out;
}

sub combine_snvs_and_indels {
    my $self = shift;
    my $cmd = "cat ".File::Spec->join($self->report_directory("snvs"), "report1.fa")." ".
        File::Spec->join($self->report_directory("indels"), "report1.fa"). " > ".$self->combined_fasta_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$self->combined_fasta_file],
    );
}

sub plan_file {
    my ($self, $variant_type) = @_;
    my $base_dir = dirname(dirname(dirname(dirname(__FILE__))));
    return File::Spec->join($base_dir, "plan_files", "flanking_fasta_$variant_type.yaml");
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    my $accessor = "get_detailed_$variant_type"."_vcf";
    return $self->build->$accessor;
}

sub log_directory {
    my ($self, $variant_type) = @_;
    return File::Spec->join($self->output_directory, "logs_$variant_type");
}

sub report_directory {
    my ($self, $variant_type) = @_;
    return File::Spec->join($self->output_directory, "report_$variant_type");
}
1;

