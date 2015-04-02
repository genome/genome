package Genome::Model::Tools::Vcf::EvaluateVcf;

use strict;
use warnings;

use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;
use Genome;
use Genome::Sys::ShellPipeline;
use IO::File;
use IO::Zlib;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Temp qw(tempdir);
use Path::Class;
use Genome::Model::Tools::Vcf::VcfCompare;


my ($REFERENCE);
my $bgzip_pipe_cmd = "| bgzip -c ";

class Genome::Model::Tools::Vcf::EvaluateVcf {
    is => "Command::V2",
    has_input => [
        bedtools_version => {
            is => "Text",
            doc => "Bedtools version to use",
            default_value => "2.17.0",
        },

        joinx_version => {
            is => "Text",
            doc => "Joinx version to use",
            default_value => "v1.10",
        },

        vcflib_version => {
            is => "Text",
            doc => "vcflib toolset version to use",
            default_value => "1.0",
        },

        vcf => {
            is => "Text",
            doc => "Vcf file to analyze",
        },

        output_directory => {
            is => "Text",
            doc => "Output directory to write to",
            is_output => 1,
        },

        old_sample => {
            is => "Text",
            doc => "Sample name in input vcf file",
        },

        new_sample => {
            is => "Text",
            doc => "Sample to use in output files",
        },

        gold_vcf => {
            is => "Text",
            doc => "Gold standard vcf to compare to",
        },

        gold_sample => {
            is => "Text",
            doc => "Sample name in gold standard vcf file",
        },

        roi => {
            is => "Text",
            doc => "Region of interest bed file",
        },

        true_negative_bed => {
            is => "Text",
            doc => "True negatives bed file",
        },

        pass_only_expression => {
            is => "Text",
            doc => "Expression for vcflib/vcffilter to select only passing (unfiltered) variants",
            default_value => q{-g 'FT = PASS | FT = .'},
        },

        clean_indels => {
            is => "Text",
            doc => "If set, attempt to cleanup indels",
            default_value => 0,
        },
    ]
};

sub execute {
    my $self = shift;

    $self->_setup_prg_tools();

    my $vcf = $self->vcf;
    my $roi = $self->roi;
    my $gold_vcf = $self->gold_vcf;
    my $tn_bed = $self->true_negative_bed;
    my $old_sample = $self->old_sample;
    my $new_sample = $self->new_sample;
    my $gold_sample = $self->gold_sample;
    my $clean_indels = $self->clean_indels;

    my $output_dir = $self->output_directory;

    die "Output dir $output_dir does not exist" unless -d $output_dir;
    my $orig_file = "$output_dir/orig.vcf";
    my $input_file = "$output_dir/clean.vcf";
    my $final_input_file = "$output_dir/final.vcf";
    my $final_gold_file = "$output_dir/gold-final.vcf";
    my $final_tn_file = "$output_dir/tn-final.bed";
    my $variants_dir = "$output_dir/variants";
    my $fp_roi_file = "$output_dir/fp_in_roi.vcf";
    my $compare_file = "$output_dir/compare.txt";

    Genome::Sys->create_directory($variants_dir);
    Genome::Sys->create_symlink($vcf, "$output_dir/orig.vcf");
    Genome::Sys->create_symlink($roi, "$output_dir/roi.bed");
    Genome::Sys->create_symlink($gold_vcf, "$output_dir/gold.vcf");

    $self->restrict($tn_bed, $roi, $final_tn_file);


    # input file processing
    my $tmp = Genome::Sys->create_temp_file_path;
    $self->_clean_vcf($orig_file, $tmp);
    $self->_process_input_file($tmp, $final_input_file);

    # Gold file processing
    #  this should be cleaned here because, presumably, allelic primitives has already been run.
    $self->restrict($gold_vcf, $roi, $tmp);
    $self->_clean_vcf($tmp, $final_gold_file);


    $self->compare_partial(
        $final_input_file,
        $variants_dir,
        $final_gold_file,
        $compare_file,
        $gold_sample,
        $old_sample,
        $new_sample
        );

    my $tn_bed_size = $self->bed_size("$tn_bed.roi.bed.gz");

    my $false_positives_in_roi = $self->number_within_roi(
        $final_input_file,
        $final_tn_file,
        $fp_roi_file,
        );

    my %results = $self->true_positives($final_input_file, $final_gold_file, $compare_file, $new_sample);

    print join("\t",
        # Exact matches
        $results{true_positive_exact},

        # All true cases
        $results{true_positive_exact} + $results{false_negative_exact},

        # Exact sensitivty
        $results{true_positive_exact} / ($results{true_positive_exact} + $results{false_negative_exact}),

        $results{true_positive_partial},
        $results{true_positive_partial} + $results{false_negative_partial},
        $results{true_positive_partial} / ($results{true_positive_partial} + $results{false_negative_partial}),
        $results{false_positive_exact},
        $results{false_positive_partial},
        $tn_bed_size,
        #exact specificity, these are not strictly accurate as the tn_bed may be significantly smaller than the target space ROI
        ($tn_bed_size - $results{false_positive_exact}) / $tn_bed_size,
        ($tn_bed_size - $results{false_positive_partial}) / $tn_bed_size,
        $results{true_positive_exact} / ($results{false_positive_exact} + $results{true_positive_exact}),
        $results{true_positive_partial} / ($results{false_positive_partial} + $results{true_positive_partial}),
        $false_positives_in_roi,
        ($tn_bed_size - $false_positives_in_roi) / $tn_bed_size, #this is more accurate than the other specificity measures, but doesn't take partial into account
    ), "\n";

    return 1;
}

sub joinx {
    my $self = shift;
    my $path = Genome::Model::Tools::Joinx->joinx_path($self->joinx_version);
    return $path;
}

sub bedtools {
    my $self = shift;
    my $base = Genome::Model::Tools::BedTools->path_for_bedtools_version(
        $self->bedtools_version
    );
    my $base = Path::Class::Dir->new($base);
    my $prg = $base->subdir('bin')->file('bedtools');
    return $prg->stringify;
}

sub vcflib_tool {
    my ($self, $tool) = @_;
    my $path = Genome::Model::Tools::Vcflib->vcflib_tool_path(
        $tool, $self->vcflib_version
    );
    return $path;
}

sub reference {
    my $self = shift;
    return "/gscmnt/ams1102/info/model_data/"
           . "2869585698/build106942997/all_sequences.fa";
}

sub _setup_prg_tools {
    my $self = shift;
    $REFERENCE = $self->reference();
}

sub _process_input_file {
    my ($self, $input_file, $output_file) = @_;
    my @cmds = (
        $self->restrict_commands($input_file, $self->roi),
        $self->restrict_to_sample_commands("/dev/stdin", $self->old_sample),
        $self->pass_only_commands("/dev/stdin", $self->pass_only_expression),
        $self->allelic_primitives_commands("/dev/stdin"),
        $self->normalize_vcf_commands("/dev/stdin", $REFERENCE),
        $self->sort_commands("/dev/stdin"),
        $self->restrict_commands("stdin", $self->roi),
        "bgzip -c",
        );

    my $cmd = Genome::Sys::ShellPipeline->create(
        pipe_commands => \@cmds,
        redirects => " > $output_file",
        );
    $cmd->execute;

}

sub _clean_vcf {
    my ($self, $input_file, $output_file) = @_;

    $self->status_message("Cleaning vcf $input_file => $output_file");
    $DB::single = 1;
    my $reader = new Genome::File::Vcf::Reader($input_file);
    my $writer = new Genome::File::Vcf::Writer($output_file, $reader->header);

    $reader->add_filter(_make_bad_indel_filter($self->roi)) if $self->clean_indels;
    $reader->add_filter(_make_per_allele_info_filter($reader->header));

    while (my $entry = $reader->next) {
        $writer->write($entry);
    }
    $writer->close;
}

sub _make_per_allele_info_filter {
    my ($header) = @_;

    my $info = $header->info_types;
    my @remove = grep {$info->{$_}{number} =~ /^(A|R)$/} keys %$info;

    return sub {
        my $entry = shift;
        my $info_hash = $entry->info;
        delete @{$entry}{@remove};
        return $entry;
    };
}

sub _make_bad_indel_filter {
    my ($roi_file) = @_;

    my %bed_ends;

    # This is compressed!
    my $fh = new IO::Zlib;
    $fh->open($roi_file, "rb") or die "Unable to open BED file $roi_file for removal of bad indel lines\n";
    while(my $bedline = $fh->getline) {
        chomp $bedline;
        my ($chr, $start, $stop) = split "\t", $bedline;
        $bed_ends{"$chr\t$stop"} = 1;
    }
    $fh->close;

    return sub {
        my $entry = shift;
        my $x = sprintf "%s\t%s", $entry->{chrom}, $entry->{position};
        return if exists $bed_ends{$x};
        return $entry;
    }
}

sub restrict {
    my ($self, $input_file, $roi_file, $output_file) = @_;
    my $bedtools = $self->bedtools;
    $DB::single = 1;
    my @args = ($bedtools, "intersect", "-header", "-a", $input_file, "-b", $roi_file, "> $output_file");
    Genome::Sys->shellcmd( cmd => join(" ", @args) );

    $DB::single = 1;
    return 1;
}


sub restrict_commands {
    my $self = shift;
    my ($input_file, $roi_file) = @_;
    my $bedtools = $self->bedtools;
    return ("$bedtools intersect -header -a $input_file -b $roi_file");
}

sub restrict_to_sample_commands {
    my $self = shift;
    my ($input_file, $sample) = @_;
    my $vcfkeepsamples = $self->vcflib_tool('vcfkeepsamples');
    return ("$vcfkeepsamples $input_file $sample");
}

sub normalize_vcf_commands {
    my $self = shift;
    my ($input_file, $reference, $output_file) = @_;
    my $joinx = $self->joinx;
    return ("$joinx vcf-normalize-indels -f $reference $input_file");
}

sub sort_commands {
    my $self = shift;
    my ($input_file) = @_;
    my $joinx = $self->joinx;
    return ("$joinx sort $input_file");
}

sub allelic_primitives_commands {
    my $self = shift;
    my ($input_file) = @_;

    my $vcfallelicprimitives = $self->vcflib_tool('vcfallelicprimitives');
    my $vcffixup = $self->vcflib_tool('vcffixup');
    my $vcffilter = $self->vcflib_tool('vcffilter');

    return (
        "$vcfallelicprimitives -t ALLELICPRIMITIVE /dev/stdin",
        "$vcffixup - ",
        "$vcffilter -f 'AC > 0'"
        );
}

sub pass_only_commands {
    my $self = shift;
    my ($input_file, $expression) = @_;
    my $vcffilter = $self->vcflib_tool('vcffilter');
    return ("$vcffilter $expression $input_file");
}


sub restrict_input_file {
    my ($input_file, $roi_file, $output_file, $sample, $filter_exp)  = @_;
    my @cmds = (
        restrict_commands($input_file, $roi_file),
        restrict_to_sample_commands("/dev/stdin", $sample),
        pass_only_commands("/dev/stdin", $filter_exp),
        allelic_primitives_commands("/dev/stdin"),
        normalize_vcf_commands("/dev/stdin", $REFERENCE),
        sort_commands("/dev/stdin"),
        restrict_commands("stdin", $roi_file),
        "bgzip -c",
        );

    my $cmd = Genome::Sys::ShellPipeline->create(
        pipe_commands => \@cmds,
        redirects => " > $output_file",
        );
    $cmd->execute;
}

sub compare_partial {
    my $self = shift;
    my ($input_file, $variant_directory, $gold_file, $output_file, $gold_sample, $eval_sample, $new_sample) = @_;
    my $rename_option = "";
    if($new_sample) {
        if($gold_sample) {
            $rename_option .= " -R $gold_sample=$new_sample";
        }
        if($eval_sample) {
            $rename_option .= " -R $eval_sample=$new_sample";
        }
    }
    my $joinx = $self->joinx;
    _run("$joinx vcf-compare $rename_option -d $variant_directory $input_file $gold_file -s $new_sample > $output_file");
}

sub true_positives {
    my $self = shift;
    my ($input_file, $gold_file, $joinx_output, $new_sample) = @_;
    my $table = Genome::Model::Tools::Vcf::VcfCompare->new($joinx_output);
    #for now only do perfect matches
    return (
        false_positive_exact => $table->unique_count($input_file, "exact_match", $new_sample),
        false_negative_exact => $table->unique_count($gold_file, "exact_match", $new_sample),
        true_positive_exact => $table->joint_count("exact_match", $new_sample),
        false_positive_partial => $table->unique_count($input_file, "partial_match", $new_sample),
        false_negative_partial => $table->unique_count($gold_file, "partial_match", $new_sample),
        true_positive_partial => $table->joint_count("partial_match", $new_sample),
        false_positive_partial_miss => $table->unique_count($input_file, "partial_miss", $new_sample),
        false_negative_partial_miss => $table->unique_count($gold_file, "partial_miss", $new_sample),
        false_positive_complete_miss => $table->unique_count($input_file, "complete_miss", $new_sample),
        false_negative_complete_miss => $table->unique_count($gold_file, "complete_miss", $new_sample),
    );
}


sub number_within_roi {
    my $self = shift;
    my ($input_file, $roi, $output_file) = @_;
    my $bedtools = $self->bedtools;
    _run("zcat $input_file | $bedtools intersect -header -a stdin -b $roi $bgzip_pipe_cmd > $output_file");
    return count($output_file);
}

sub bed_size {
    my $self = shift;
    my $bed = shift;
    my $count = 0;
    my $fh = IO::File->new("zcat $bed |") or die "Unable to open $bed to calculate size\n";
    while(my $line  = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop) = split "\t", $line;
        $count += $stop-$start;
    }
    $fh->close;
    return $count;
}

sub count {
    my $file = shift;
    my @results = `zgrep -v '^#' $file | cut -f1,2,3 | sort -u | wc -l`;
    chomp $results[0];
    return $results[0];
}

sub _run {
    my $cmd = shift;

    return Genome::Sys->shellcmd(cmd => $cmd);
}

1;
