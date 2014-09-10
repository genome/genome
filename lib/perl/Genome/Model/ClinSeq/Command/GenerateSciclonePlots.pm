package Genome::Model::ClinSeq::Command::GenerateSciclonePlots;

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::GenerateSciclonePlots {
    is => 'Command::V2',
    has_input => [
        clinseq_build => {
            is => 'Genome::Model::Build::ClinSeq',
            doc => 'ClinSeq build to make SciClone plots for.',
        },
        outdir => {
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written.',
        },
        maximum_clusters => {
            is => 'Number',
            is_optional => 1,
            doc => 'Maximum number of clusters.',
        },
        minimum_depth => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum depth of variants.',
        },
        microarray_cnv_result => {
            is => 'Boolean',
            doc => 'Link in workflow',
            is_optional => 1
        },
        exome_cnv_result => {
            is => 'Boolean',
            doc => 'Link in workflow',
            is_optional => 1
        },
        wgs_cnv_result => {
            is => 'Boolean',
            doc => 'Link in workflow',
            is_optional => 1
        },
        converge_snv_indel_report_result => {
            is => 'Boolean',
            doc => 'Link in workflow',
            is_optional => 1
        },
    ],
    doc => 'Create clonality plots with SciClone.',
};

sub help_synopsis {
    return <<EOS
        genome model clin-seq generate-sciclone-plots \\
        --outdir=/gscuser/gscuser1/tmp/ \\
        --clinseq-build='a4abcd1313eb4376b59e68a9dd9d5ad2'
EOS
}

sub help_detail {
    return <<EOS
Generate SciClone Plots inside clin-seq builds.
EOS
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    unless (-e $self->outdir && -d $self->outdir) {
        push @errors, UR::Object::Tag->create(
          type => 'error',
          properties => ['outdir'],
          desc => "Outdir: " . $self->outdir .
            " not found or not a directory",
        );
    }
    return @errors;
}

sub parse_variant_file {
    my $self = shift;
    my $clinseq_build = shift;
    my $variant_file = shift;
    my $variant_file_temp = $variant_file . ".tmp";
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input => $variant_file,
    );
    my @headers = qw/chr pos ref_allele var_allele ref_rc var_rc vaf/;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $variant_file_temp,
        separator => "\t",
        headers => \@headers,
        print_headers => 0,
    );
    my $out_data;
    my $tumor_prefix = $self->_get_si_report_tumor_prefix($clinseq_build);
    while (my $data = $reader->next) {
        if ($data->{chromosome_name} =~ /X|Y|MT/) {
            next;
        }
        $out_data->{chr} = $data->{chromosome_name};
        $out_data->{pos} = $data->{start};
        $out_data->{ref_allele} = $data->{reference};
        $out_data->{var_allele} = $data->{variant};
        $out_data->{ref_rc} = $data->{$tumor_prefix . "_ref_count"};
        $out_data->{var_rc} = $data->{$tumor_prefix . "_var_count"};
        $out_data->{vaf} = $data->{$tumor_prefix . "_VAF"};
        $writer->write_one($out_data);
    }
    unlink $variant_file;
    Genome::Sys->move_file($variant_file_temp, $variant_file);
    return $variant_file;
}


sub get_variant_file {
    my $self = shift;
    my $clinseq_build = shift;
    my $outfile = $self->outdir . "/variants.clean.tsv";
    my $snv_indel_report_clean_file =
        $clinseq_build->snv_indel_report_clean_filtered_file;
    if(-e $snv_indel_report_clean_file) {
        Genome::Sys->copy_file($snv_indel_report_clean_file, $outfile);
    } else {
        die $self->error_message("Unable to find variant read-counts file
            for clinseq build " . $clinseq_build->id);
    }
    return $self->parse_variant_file($clinseq_build, $outfile);
}

sub parse_cnv_file {
    my $self = shift;
    my $clinseq_build = shift;
    my $cnv_file = shift;
    my $cnv_file_temp = $cnv_file . ".tmp";
    my @cnvhmm_header = qw/chr start end size adjusted_size tumor_cn
        tumor_adjusted_cn normal_cn normal_adjusted_cn unknown_column status/;
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input => $cnv_file,
        headers => \@cnvhmm_header,
    );
    my @headers = qw/chr start end num_probes copy_number/;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $cnv_file_temp,
        separator => "\t",
        headers => \@headers,
        print_headers => 0,
    );
    my $out_data;
    while (my $data = $reader->next) {
        if ($data->{chr} =~ /CHR|X|Y|MT/) {
            next;
        }
        $out_data->{chr} = $data->{chr};
        $out_data->{start} = $data->{start};
        $out_data->{end} = $data->{end};
        $out_data->{num_probes} = "NA";
        $out_data->{copy_number} = $data->{tumor_adjusted_cn};
        $writer->write_one($out_data);
    }
    unlink $cnv_file;
    Genome::Sys->move_file($cnv_file_temp, $cnv_file);
    return $cnv_file;
}

sub get_cnv_file {
    my $self = shift;
    my $clinseq_build = shift;
    my $outfile = $self->outdir . "/cnvs.tsv";
    my $cnv_file = $clinseq_build->best_cnvhmm_file;
    if(-e $cnv_file) {
        Genome::Sys->copy_file($cnv_file, $outfile);
    } else {
        die $self->error_message("Unable to find CNV file " .
            "for build " . $clinseq_build->id);
    }
    return $self->parse_cnv_file($clinseq_build, $outfile);
}

sub run_sciclone {
    my $self = shift;
    my $variant_f = shift;
    my $cnv_f = shift;
    my $clinseq_build = shift;
    my $outdir = $self->outdir;
    my $clusters_f = $outdir . "/sciclone.clusters.txt";
    my $rscript_f = $outdir . "/sciclone.R";
    my $plot_f = $outdir . "/sciclone.clonality.pdf";
    my $maximum_clusters;
    if ($self->maximum_clusters) {
        $maximum_clusters = $self->maximum_clusters;
    } else {
        $maximum_clusters = 6;
    }
    my $minimum_depth;
    if ($self->minimum_depth) {
        $minimum_depth = $self->minimum_depth;
    } else {
        $minimum_depth = 1;
    }
    my $sample_name = $clinseq_build->subject->name;
    my $sciclone = Genome::Model::Tools::Sciclone->create(
        clusters_file => $clusters_f,
        r_script_file => $rscript_f,
        variant_files => $variant_f,
        copy_number_files => $cnv_f,
        maximum_clusters => $maximum_clusters,
        minimum_depth => $minimum_depth,
        sample_names => $sample_name,
        plot1d_file => $plot_f,
        cn_calls_are_log2 => 1,
        label_highlighted_points => 1,
        plot_only_cn2 => 1,
    );
    $sciclone->execute();
}

#SciClone needs two files
# 1. Variant file of format (chr pos ref_all var_all ref_RC var_RC VAF)
# 2. CNV file of format (chr start stop num_probes copy_number)
sub execute {
    my $self = shift;
    my $clinseq_build = $self->clinseq_build;
    my $variant_f = $self->get_variant_file($clinseq_build);
    my $cnv_f = $self->get_cnv_file($clinseq_build);
    $self->run_sciclone($variant_f, $cnv_f, $clinseq_build);
    return 1;
}

1;
