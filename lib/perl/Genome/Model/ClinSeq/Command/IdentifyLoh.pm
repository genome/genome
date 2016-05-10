package Genome::Model::ClinSeq::Command::IdentifyLoh;

use strict;
use warnings;
use File::Spec;
use Genome;

class Genome::Model::ClinSeq::Command::IdentifyLoh {
    is        => ['Command::V2', 'Genome::Model::ClinSeq::Util',],
    has_input => [
        clinseq_build => {
            is  => 'Genome::Model::Build::ClinSeq',
            doc => 'ClinSeq build to identify LOH regions in.
            [Either this or a somatic build is required.]',
            is_optional => 1,
        },
        outdir => {
            is  => 'FilesystemPath',
            doc => 'Directory where output files will be written.',
        },
        somatic_build => {
            is  => 'Genome::Model::Build::SomaticInterface',
            doc => 'SomVar build to identify LOH regions in.
        [Either this or a clinseq build is required.]',
            is_optional => 1,
        },
        minprobes => {
            is          => 'Number',
            is_optional => 1,
            doc         => 'Minimum number of probes to call a loh segment.',
            default     => 10,
        },
        min_segmean => {
            is          => 'Number',
            is_optional => 1,
            doc         => 'Minimum segment mean. This is the percent of LOH SNPs in a segment',
            default     => 0.95,
        },
        bamrc_version => {
            is  => 'String',
            doc => 'version of bam-readcount to use',
        },
        test => {
            is          => 'Boolean',
            doc         => 'set for test-cases',
            is_optional => 1,
            default     => 0,
        },
    ],
    doc => 'Identify regions of LOH using clinseq tumor/normal pairs.',
};

sub help_synopsis {
    return <<EOS
genome model clin-seq identify-loh\\
    --outdir=/tmp/test/ \\
    --clinseq-build='a4abcd1313eb4376b59e68a9dd9d5ad2' \\
    --bamrc-version=0.7
genome model clin-seq identify-loh\\
    --outdir=/tmp/test/ \\
    --somatic-build='9dc0385a1b634c9bb85eb2017b3a0c73' \\
    --bamrc-version=0.7
EOS
}

sub help_detail {
    return <<EOS
Use the results from Varscan in the somatic builds and identify
regions of LOH. Specifically, look at heterozygous single nucleotide variants
in the normal sample and compare with stretches of homozygosity in the tumor.

If a clinseq build is supplied as input, the tool attempts to use the results
from the WGS-somatic build input, if these are unavailable the tool
looks for WEx-somatic build.

Alternatively, instead of a clin-seq build you can supply a
somatic build as an input.

EOS
}

sub __errors__ {
    my $self   = shift;
    my @errors = $self->SUPER::__errors__(@_);
    unless (-e $self->outdir && -d $self->outdir) {
        push @errors,
            UR::Object::Tag->create(
            type       => 'error',
            properties => ['outdir'],
            desc       => "Outdir: " . $self->outdir . " not found or not a directory",
            );
    }
    return @errors;
}

#Get a somatic build from input or from clinseq-build
sub resolve_somatic_build {
    my $self         = shift;
    if ($self->somatic_build) {
        return $self->somatic_build;
    }
    else {
        return $self->clinseq_build->best_somatic_build;
    }
}

#Take Varscan SNVs and split into per-chromosome calls.
sub split_snvs {
    my $self       = shift;
    my $varscan_op = shift;
    my $snv_prefix = shift;
    my $awk_cmd;
    #use subsample of chr1 variants for test
    if ($self->test) {
        $awk_cmd = "awk '\$1==1 { if(NR\%1000 == 0) { print > \"$snv_prefix.\"\$1\".unfiltered\" }}' $varscan_op";
    }
    else {
        $awk_cmd = "awk '{ print > \"$snv_prefix.\"\$1\".unfiltered\" }' $varscan_op";
    }
    Genome::Sys->shellcmd(cmd => $awk_cmd);
    return;
}

#Get Varscan SNVs from somatic-variation
sub get_varscan_snvs {
    my $self         = shift;
    my $somatic_build = shift;
    my $snv_prefix   = shift;
    my $dd           = $somatic_build->data_directory;
    my @varscan_snvs = glob(File::Spec->join($dd, "variants/snv/varscan-somatic-*/snvs.hq.unfiltered"));
    unless (scalar @varscan_snvs) {
        die $self->error_message("Unable to find varscan SNVs in $dd");
    }
    $self->status_message("varscan snvs file found: " . join("\t", @varscan_snvs));
    $self->split_snvs($varscan_snvs[0], $snv_prefix);
}

#Filter out putative false-positive SNVs
sub filter_snvs {
    my $self         = shift;
    my $somatic_build = shift;
    my $snv_prefix   = shift;
    my $refbuild = $somatic_build->reference_sequence_build;
    my $ref_fa   = $refbuild->full_consensus_path('fa');
    my $filter   = Genome::Model::Tools::Varscan::SomaticFilterWorkflow->create(
        outdir        => $self->outdir,
        tumor_bam     => $somatic_build->tumor_bam,
        normal_bam    => $somatic_build->normal_bam,
        prefix        => $snv_prefix,
        reference     => $ref_fa,
        bamrc_version => $self->bamrc_version,
    );
    $filter->execute();
}

#Combine germline, LOH calls and sort them
sub combine_sort_snvs {
    my $self                = shift;
    my $snv_prefix          = shift;
    my $snv_combined        = $snv_prefix . ".filtered";
    my $snv_combined_sorted = $snv_prefix . ".filtered.sorted";
    my @germline            = glob $snv_prefix . "*.formatted.LOH.hc.filtered";
    my @loh                 = glob $snv_prefix . "*.formatted.Germline.hc.filtered";
    unless (scalar @germline and scalar @loh) {
        die $self->error_message(
            "Unable to find germline or
            loh filtered results"
        );
    }
    Genome::Sys->cat(
        input_files => [@germline, @loh],
        output_file => $snv_combined
    );
    my $sort = Genome::Model::Tools::Capture::SortByChrPos->create(
        input_file  => $snv_combined,
        output_file => $snv_combined_sorted,
    );
    $sort->execute;
    return $snv_combined_sorted;
}

#Segment the combined calls
sub segment_loh {
    my $self                 = shift;
    my $snvs_combined_sorted = shift;
    my $loh_basename         = shift;
    my $segmenter            = Genome::Model::Tools::Varscan::LohSegments->create(
        variant_file    => $snvs_combined_sorted,
        output_basename => $loh_basename,
    );
    $segmenter->execute();
}

#Filter out LOH segments without sufficient support
sub filter_loh {
    my $self                  = shift;
    my $loh_basename          = shift;
    my $min_segmean           = $self->min_segmean;
    my $minprobes             = $self->minprobes;
    my $loh_segments          = $loh_basename . ".segments.cbs";
    my $loh_segments_filtered = $loh_segments . ".filtered";
    unless (-e $loh_segments) {
        die $self->error_message("Unable to find loh file $loh_segments");
    }
    my $awk_filter = "awk '\$4 > $minprobes && \$5 >= $min_segmean' $loh_segments " . "> $loh_segments_filtered";
    Genome::Sys->shellcmd(cmd => $awk_filter);
}

#Get rid of intermediate files
sub cleanup {
    my $self     = shift;
    my @patterns = qw(*.Somatic* *.readcounts *.hc *.LOH
        *.lc *.Germline *.removed *.formatted
        *.hc.err *.other *.err *.out *.snp snvs.filtered snvs.*.unfiltered *hc.filtered);
    for my $pattern (@patterns) {
        for my $file (glob File::Spec->join($self->outdir, $pattern)) {
            $self->status_message("unlinking $file");
            unlink $file;
        }
    }
}

#Get, Filter, Combine, Segment, Filter, Cleanup
sub execute {
    my $self         = shift;
    my $somatic_build = $self->resolve_somatic_build;
    my $snv_prefix   = File::Spec->join($self->outdir, "snvs");
    $self->get_varscan_snvs($somatic_build, $snv_prefix);
    $self->filter_snvs($somatic_build, $snv_prefix);
    my $combined_sorted = $self->combine_sort_snvs($snv_prefix);
    my $loh_basename    = File::Spec->join($self->outdir, "loh");
    my $loh_segments    = $self->segment_loh($combined_sorted, $loh_basename);
    $self->filter_loh($loh_basename);
    $self->cleanup;
    return 1;
}

1;

