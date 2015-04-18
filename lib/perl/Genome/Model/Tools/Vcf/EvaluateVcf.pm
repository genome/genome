package Genome::Model::Tools::Vcf::EvaluateVcf;

use strict;
use warnings;

use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;
use Genome;
use Genome::Sys::ShellPipeline;
use IO::File;
use IO::Zlib;
use Path::Class;
use Genome::Model::Tools::Vcf::VcfCompare;

class Genome::Model::Tools::Vcf::EvaluateVcf {
    is => "Command::V2",
    has_input => [
        reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'a reference sequence of interest.',
        },

        bedtools_version => {
            is => "Text",
            doc => "Bedtools version to use",
            default_value => "2.17.0",
        },

        joinx_version => {
            is => "Text",
            doc => "Joinx version to use",
            default_value => "1.11",
        },

        vcflib_version => {
            is => "Text",
            doc => "vcflib toolset version to use",
            default_value => "1.0",
        },

        vcf => {
            is => "Text",
            doc => "VCF file to be evaluated",
        },

        output_directory => {
            is => "Text",
            doc => "Output directory to write to",
            is_output => 1,
        },

        old_sample => {
            is => "Text",
            doc => "Sample name in the VCF file to use",
        },

        new_sample => {
            is => "Text",
            doc => "Sample name to call the old sample during comparisons",
        },

        gold_vcf => {
            is => "Text",
            doc => "VCF file containing gold standard variant sites "
                   . "to be used to measure the evaluation VCF"
        },

        gold_sample => {
            is => "Text",
            doc => "Sample name in gold standard vcf file",
        },

        roi => {
            is => "Text",
            doc => "BED file of the target regions to restrict the analysis to",
        },

        true_negative_bed => {
            is => "Text",
            doc => "BED file containing regions where no "
                   . "variant call should be made",
        },

        true_negative_size => {
            is => "Integer",
            doc => "Use this number as the size of the TN BED "
                   . "file rather than calculating on the fly",
            is_optional => '1',
        },

        pass_only_expression => {
            is => "Text",
            doc => "String to pass to vcflib vcffilter in order to "
                   . "select variants from the evaluation VCF",
            default_value => q{-g 'FT = PASS | FT = .'},
        },

        clean_indels => {
            is => "Text",
            doc => "Whether or not to exclude indels that don't "
                   . "actually overlap the ROI",
            default_value => 0,
        },
    ],

    has_transient_optional_output => [
        rawstats => {
            is => "HASH",
            doc => "The raw stats generated during primary execution",
        },
    ],
};

sub execute {
    my $self = shift;

    my $vcf = $self->vcf;
    my $roi = $self->roi;
    my $gold_vcf = $self->gold_vcf;
    my $tn_bed = $self->true_negative_bed;
    my $old_sample = $self->old_sample;
    my $new_sample = $self->new_sample;
    my $gold_sample = $self->gold_sample;
    my $clean_indels = $self->clean_indels;

    my $output_dir = Path::Class::Dir->new($self->output_directory);

    die "Output dir $output_dir does not exist" unless -d $output_dir;
    my $orig_file = $output_dir->file("orig.vcf")->stringify;
    my $input_file = $output_dir->file("clean.vcf")->stringify;
    my $final_input_file = $output_dir->file("final.vcf")->stringify;
    my $final_gold_file = $output_dir->file("gold-final.vcf")->stringify;
    my $final_tn_file = $output_dir->file("tn-final.bed")->stringify;
    my $variants_dir = $output_dir->subdir("variants")->stringify;
    my $fp_roi_file = $output_dir->file("fp_in_roi.vcf")->stringify;
    my $compare_file = $output_dir->file("compare.txt")->stringify;

    Genome::Sys->create_directory($variants_dir);
    Genome::Sys->create_symlink($vcf, $orig_file);
    Genome::Sys->create_symlink($roi, $output_dir->file("roi.bed")->stringify);
    Genome::Sys->create_symlink($gold_vcf, $output_dir->file("gold.vcf")->stringify);

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

    my $tn_bed_size = $self->true_negative_size
      || $self->bed_size(${tn_bed});

    my $false_positives_in_roi = $self->number_within_roi(
        $final_input_file,
        $final_tn_file,
        $fp_roi_file,
        );

    my %results = $self->true_positives(
        $final_input_file,
        $final_gold_file,
        $compare_file,
        $new_sample
    );

    $results{true_negatives} = $tn_bed_size;
    $results{false_positives_in_roi} = $false_positives_in_roi;

    $self->rawstats(\%results);

#    $self->display_all_stats();

    return 1;
}

sub display_all_stats {
    my $self = shift;
    print join("\t",
        $self->stat_true_positive_found_exact(),
        $self->stat_total_true_positive_exact(),
        $self->stat_sensitivity_exact(),
        $self->stat_true_positive_found_partial(),
        $self->stat_total_true_positive_partial(),
        $self->stat_sensitivity_partial(),
        $self->stat_false_positive_exact(),
        $self->stat_false_positive_partial(),
        $self->stat_true_negatives(),
        $self->stat_exact_specificity(),
        $self->stat_partial_specificity(),
        $self->stat_exact_ppv(),
        $self->stat_partial_ppv(),
        $self->stat_vcf_lines_overlapping_tn(),
        $self->stat_lines_specificity_in_tn_only(),
    ), "\n";
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
    $base = Path::Class::Dir->new($base);
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

sub reference_path {
    my $self = shift;
    my $path = $self->reference->full_consensus_path('fa');
    return $path;
}

sub _process_input_file {
    my ($self, $input_file, $output_file) = @_;
    my @cmds = (
        $self->restrict_commands($input_file, $self->roi),
        $self->restrict_to_sample_commands("/dev/stdin", $self->old_sample),
        $self->pass_only_commands("/dev/stdin", $self->pass_only_expression),
        $self->allelic_primitives_commands("/dev/stdin"),
        $self->normalize_vcf_commands("/dev/stdin", $self->reference_path),
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
    my $reader = Genome::File::Vcf::Reader->new($input_file);
    my $writer = Genome::File::Vcf::Writer->new($output_file, $reader->header);

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
    my $fh = IO::Zlib->new;
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
    my $cmd = $self->restrict_commands($input_file, $roi_file);
    $cmd .= " > $output_file";
    Genome::Sys->shellcmd( cmd => $cmd );

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
        "$vcfallelicprimitives -t ALLELICPRIMITIVE $input_file",
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
    my $bgzip_pipe_cmd = "| bgzip -c ";
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

# H E L P E R  F U N C T I O N S ##############################################
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

# S T A T  /  M E T R I C   A C C E S S O R S #################################
sub get_rawstats {
    my $self    = shift;
    my $results = $self->rawstats or die $self->error_message(
        "Please invoke 'execute' to generate rawstats!"
    );
    return $results;
}

sub stat_true_positive_found_exact {
    my $self = shift;
    my $results = $self->get_rawstats();
    my $stat = $results->{true_positive_exact};
    return $stat;
}

sub stat_total_true_positive_exact {
    my $self    = shift;
    my $results = $self->get_rawstats();
    my $stat =
      $results->{true_positive_exact} + $results->{false_negative_exact};
    return $stat;
}

sub stat_sensitivity_exact {
    my $self    = shift;
    my $results = $self->get_rawstats();

    my $total_true = $self->stat_total_true_positive_exact;
    if ($total_true == 0) {
        return 0;
    }

    my $stat = $results->{true_positive_exact} / $total_true;
    return $stat;
}

sub stat_true_positive_found_partial {
    my $self    = shift;
    my $results = $self->get_rawstats();
    my $stat    = $results->{true_positive_partial};
    return $stat;
}

sub stat_total_true_positive_partial {
    my $self    = shift;
    my $results = $self->get_rawstats();
    my $stat    = 
        $results->{true_positive_partial} + $results->{false_negative_partial};
    return $stat;
}

sub stat_sensitivity_partial {
    my $self    = shift;
    my $results = $self->get_rawstats();

    my $total_true = $self->stat_total_true_positive_partial;
    if ($total_true == 0) {
        return 0;
    }

    my $stat = $results->{true_positive_partial} / $total_true;
    return $stat;
}

sub stat_false_positive_exact {
    my $self = shift;
    my $results = $self->get_rawstats();
    my $stat = $results->{false_positive_exact};
    return $stat;
}

sub stat_false_positive_partial {
    my $self = shift;
    my $results = $self->get_rawstats();
    my $stat = $results->{false_positive_partial};
    return $stat;
}

sub stat_true_negatives {
    my $self = shift;
    my $results = $self->get_rawstats();
    my $stat = $results->{true_negatives};
    return $stat;
}

sub stat_exact_specificity {
    my $self = shift;
    my $results = $self->get_rawstats();

    # exact specificity, these are not strictly accurate as the tn_bed
    # may be significantly smaller than the target space ROI
    
    my $tn = $self->stat_true_negatives();
    if ($tn == 0) {
        return 0;
    }

    my $stat =
      ($results->{true_negatives} - $results->{false_positive_exact}) / $tn;
    return $stat;
}

sub stat_partial_specificity {
    my $self    = shift;
    my $results = $self->get_rawstats();

    my $tn = $self->stat_true_negatives();
    if ($tn == 0) {
        return 0;
    }

    my $stat =
      ($results->{true_negatives} - $results->{false_positive_partial}) / $tn;
    return $stat;
}

sub stat_exact_ppv {
    my $self    = shift;

    my $fp_exact = $self->stat_false_positive_exact();
    my $tp_exact = $self->stat_true_positive_found_exact();

    my $denominator = $fp_exact + $tp_exact;
    if ($denominator == 0) {
        return 0;
    }

    my $stat = $tp_exact / $denominator;
    return $stat;
}

sub stat_partial_ppv {
    my $self    = shift;

    my $fp_partial = $self->stat_false_positive_partial();
    my $tp_partial = $self->stat_true_positive_found_partial();

    my $denominator = $fp_partial + $tp_partial;
    if ($denominator == 0) {
        return 0;
    }

    my $stat = $tp_partial / $denominator;
    return $stat;
}

sub stat_vcf_lines_overlapping_tn {
    my $self    = shift;
    my $results = $self->get_rawstats();
    my $stat    = $results->{false_positives_in_roi};
    return $stat;
}

sub stat_lines_specificity_in_tn_only {
    my $self    = shift;
    my $results = $self->get_rawstats();

    # this is more accurate than the other specificity measures, but
    # doesn't take partial into account

    my $tn = $self->stat_true_negatives();
    if ($tn == 0) {
        return 0;
    }

    my $stat = ($tn - $results->{false_positives_in_roi}) / $tn;
    return $stat;
}

sub help_brief {
    return "Compare or Validate VCF files to a Gold Standard";
}

sub help_detail {
    my $doc = q{

    BACKGROUND
    ==========

    Modern genomic pipelines are composed of continually advancing next
    generation sequencing technologies and software tools.  As pipelines
    evolve, it is important, especially in a clinical setting, to
    systematically assess the accuracy and reproducibility of the variant
    calls that are produced as the end product.

    One way to measure the overall accuracy of a genomic pipeline is to
    compare variant calls generated from a special input data set with
    already known variants.  The special data set is called a "gold
    standard", and the process of comparing the pipeline obtained variant
    calls to the known gold standard variant calls is called "validation".

    The validation process categorizes whether a given variant call
    produced by a genomic pipeline is correctly identified or missed.  The
    classifications can be measured overall by the statistical notions of
    sensitivity, specificity, and positive predictive value (PPV).

    The standard format for storing variant calls is a VCF file.
    A procedure to validate VCF files against each other is performed in
    this script.  It roughly follows the NIST recommendations as described
    in the [Genome in a Bottle paper][1].

    METHODOLOGY
    ===========

    This tool compares a gold-standard VCF to an experimental VCF and a BED
    file of true negative positions.

    Genotype comparisons are done using [joinx vcf-compare][2]. This program
    reports an exact and a partial genotype match to the gold-standard VCF.
    Exact matches are straightforward: if both samples report the variant
    then it is an exact match. Partial matches are more complex: joinx
    reports the number of unique, matching alternate alleles at each site
    giving a zygosity-independent measure of concordance. For example, a
    heterozygous variant in the gold sample and a homozygous variant in the
    evaluation sample would yield a single partial match.

    In order to calculate metrics, the following calculations are used and
    are further detailed in Output Formats below:

    Metric                      Calculation
    -------------------------   --------------
    Specificity                 TN / (TN + FP)
    Sensitivity                 TP / (TP + FN)
    Positive Predictive Value   TP / (TP + FP)

    where 
    
        TP = True Positives
        FP = False Positives
        TN = True Negatives
        
    The TP, FP and TN counts are explicitly defined as:

    True Positive
    -------------

        Variant in the Evaluation VCF matches the variant in the Gold Standard
        VCF (either exactly or partially)

    False Positive
    --------------

        Variant present in the Evaluation VCF that is not present in the Gold
        Standard VCF

    True Negative
    -------------

        Position known to not contain a variant in the Gold Standard sample
        that is also not called in the Evaluation VCF

    Workflow
    ========

    The high-level algorithm of this tool is:
    
    1. Remove CAF fields if they exist.
    2. Restrict evaluation VCF to the ROI.
    3. Restrict gold VCF to the ROI.
    4. Restrict true negative BED to the ROI.
    5. Restrict evaluation VCF to the requested sample.
    6. Use the --pass-only-expression to restrict the evaluation VCF to only
       those variants passing the expression.
    7. Break complex indels into simpler ones using vcflib vcfallelicprimitives.
    8. Re-sort the file.
    9. Re-restrict to the ROI.
    10. Compare the resulting VCF to the ROI-restricted gold VCF using
        joinx vcf-compare.
    11. Calculate metrics and print output.

    Results or Output Statistics
    ============================

    The output statistics produced by the execute method are available as 
    accessor methods.  The available statistics are:

    True_Positive_Found_Exact

       Exact Genotype Matches between the Gold and Eval VCF

    Total_True_Positive_Exact

       Total Genotypes in the Gold VCF that were evaluated

    Sensitivity_Exact

       True_Positive_Found_Exact / Total_True_Positive_Exact

    True_Positive_Found_Partial

       At least one alternative allele matches between the Gold and Eval VCF

    Total_True_Positive_Partial

       Total Partial Genotypes in the Gold VCF that were evaluated

    Sensitivity_Partial

       True_Positive_Found_Partial / Total_True_Positive_Partial

    False_Positive_Exact

       Eval VCF sites that didn't match a site in the Gold VCF

    False_Positive_Partial

       Eval VCF alleles that didn't match a site in the Gold VCF

    True_Negatives

       Number of bases in the True Negative BED file

    Exact_Specificity
   
       (True_Negatives - False_Positive_Exact) / True_Negatives

    Partial_Specificity

       (True_Negatives - False_Positive_Partial) / True_Negatives

    Exact_PPV

       True_Positive_Found_Exact / (False_Positive_Exact + True_Positive_Found_Exact)
   
    Partial_PPV 

       True_Positive_Found_Partial / (False_Positive_Partial + True_Positive_Found_Partial)

    VCF_Lines_Overlapping_TN

       VCF lines in the Evaluation VCF found inside of the True Negative
       BED file. More conservative definition of False Positive.

    Lines_Specificity_in_TN_Only

       (True_Negatives - VCF_Lines_Overlapping_TN) / True_Negatives


    References
    ==========

    [1]: http://www.nature.com/nbt/journal/v32/n3/full/nbt.2835.html
    [2]: https://github.com/genome/joinx

    };

    return $doc;
}

1;

__END__
