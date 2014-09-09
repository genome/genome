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
            doc => 'Directory where output files will be written',
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
    $self->parse_variant_file($clinseq_build, $outfile);
}

sub get_cnv_file {
    my $self = shift;
    my $clinseq_build = shift;
    my $outfile = $self->outdir . "/variants.clean.tsv";
    my $cnv_file = $clinseq_build->best_cnv_file;
    if(-e $cnv_file) {
        Genome::Sys->copy_file($snv_indel_report_clean_file, $outfile);
    } else {
        die $self->error_message("Unable to find variant read-counts file 
            for clinseq build " . $clinseq_build->id);
    }
    $self->parse_variant_file($clinseq_build, $outfile);
}

#SciClone needs two files
# 1. Variant file of format (chr pos ref_all var_all ref_RC var_RC VAF)
# 2. CNV file of format (chr start stop num_probes copy_number)
sub execute {
    my $self = shift;
    my $clinseq_build = $self->clinseq_build;
    $self->get_variant_file($clinseq_build);
    $self->get_cnv_file();
    #$self->run_sciclone();
    return 1;
}

1;
