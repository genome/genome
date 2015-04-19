package Genome::Model::Tools::Vcf::EvaluateVcfs;

use strict;
use warnings;

use Cwd;
use Path::Class;

class Genome::Model::Tools::Vcf::EvaluateVcfs {
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
            default_value => "1.11",
        },

        vcflib_version => {
            is => "Text",
            doc => "vcflib toolset version to use",
            default_value => "1.0",
        },

        config_file => {
            is => 'Text',
            doc => 'config file input',
        },

        output_directory => {
            is => "Text",
            doc => "Output directory to write to",
            is_output => 1,
        },

        gold_snv_vcf => {
            is => "Text",
            doc => "VCF file containing gold standard variant sites "
                   . "to be used to measure the evaluation VCF (for snvs)"
        },

        gold_indel_vcf => {
            is => "Text",
            doc => "VCF file containing gold standard variant sites "
                   . "to be used to measure the evaluation VCF (for indels)"
        },

        gold_sample => {
            is => "Text",
            doc => "Sample name in gold standard vcf file",
            default => undef,
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

        pass_filter_expression => {
            is => "Text",
            doc => "String to pass to vcflib vcffilter in order to "
                   . "select variants from the evaluation VCF",
            default_value => q{-g 'FT = PASS | FT = .'},
        },
    ],

    has_output => [
        rawdata => {
            is => "HASH",
            doc => "The organized stats & metadata generated during "
                   . "primary execution",
        },
    ],
};

sub execute {
    my $self = shift;
    my $configs = $self->parse_config_file();
    $self->collect_statistics($configs);
    $self->rawdata($configs);
    $self->display_summary_stats($configs);
    return 1;
}

sub display_summary_stats {
    my ($self, $configs) = @_;

    print "\n";
    print "==================", "\n";
    print "Valdiation Summary", "\n";
    print "==================", "\n";
    print "\n";

    my @stat_names = $self->stat_types;
    my @headers = ('Name', 'ID', 'VarType', @stat_names);
    print join("\t", @headers), "\n";

    for my $set (@{$configs}) {
        my ($name, $id, $type, $stats) =
          @{$set}{'name', 'id', 'variant_type', 'stats'};
        my @results = @{$stats}{@stat_names};
        print join("\t", $name, $id, $type, @results), "\n";
    }
}

sub collect_statistics {
    my ($self, $inputs) = @_;

    for my $set (@{$inputs}) {
        my $wkspace = $self->create_subdirectory_workspace($set);
        my $stats = $self->_run_evaluate_vcf($wkspace, $set);
        $set->{'stats'} = $stats;
    }
}

sub parse_config_file {
    my $self = shift;
    my $config = $self->fetch_config_file;

    my $fh = $config->openr;

    my @inputs;
    while (my $line = <$fh>) {
        chomp $line;
        my ($name, $id, $path, $variant_type, $sample) = split("\t", $line);
        $self->validate_type($variant_type);
        my $gold_vcf =
          $variant_type eq "snvs" ? $self->gold_snv_vcf : $self->gold_indel_vcf;

        my $clean_indels = $variant_type eq 'indels' ? 1 : 0;

        my $model = Genome::Model->get($id);
        my $build = $model->last_succeeded_build;
        my $reference = $build->reference_sequence_build
          or die $self->error_message(
              "Did not find a reference for build %s", 
              $build->id
          );

        my $i = {
            name => $name,
            id => $id,
            variant_type => $variant_type,
            vcf => Path::Class::File->new($path)->absolute,
            gold_vcf => Path::Class::File->new($gold_vcf)->absolute,
            sample => $sample,
            reference => $reference,
            clean_indels => $clean_indels,
        };
        push(@inputs, $i);
    }

    return \@inputs;
}

sub _run_evaluate_vcf {
    my ($self, $outdir, $inputs) = @_;

    my %params = (
        vcf                  => $inputs->{'vcf'}->stringify,
        old_sample           => $inputs->{'sample'},
        new_sample           => "GOLDSTANDARD_SAMPLE",
        gold_vcf             => $inputs->{'gold_vcf'}->stringify,
        gold_sample          => $self->gold_sample,
        roi                  => Cwd::abs_path($self->roi),
        true_negative_bed    => Cwd::abs_path($self->true_negative_bed),
        output_directory     => $outdir->stringify,
        reference            => $inputs->{'reference'},
        pass_only_expression => $self->pass_filter_expression,
        clean_indels         => $inputs->{'clean_indels'},
    );

    if ($self->true_negative_size) {
        $params{'true_negative_size'} = $self->true_negative_size;
    }

    my $cmd = Genome::Model::Tools::Vcf::EvaluateVcf->create(%params);

    $cmd->execute or die $self->error_message(
        "Trouble running Genome::Model::Tools::Vcf::EvaluateVcf!"
    );

    my $stats = $self->_collect_evaluation_stats($cmd);
    return $stats;
}

sub _collect_evaluation_stats {
    my ($self, $cmd) = @_;

    my %stats;

    $stats{'true_positive_found_exact'} =
      $cmd->stat_true_positive_found_exact();

    $stats{'total_true_positive_exact'} =
      $cmd->stat_total_true_positive_exact();

    $stats{'sensitivity_exact'} = $cmd->stat_sensitivity_exact();

    $stats{'true_positive_found_partial'} =
      $cmd->stat_true_positive_found_partial();

    $stats{'total_true_positive_partial'} =
      $cmd->stat_total_true_positive_partial();

    $stats{'sensitivity_partial'}    = $cmd->stat_sensitivity_partial();
    $stats{'false_positive_exact'}   = $cmd->stat_false_positive_exact();
    $stats{'false_positive_partial'} = $cmd->stat_false_positive_partial();
    $stats{'true_negatives'}         = $cmd->stat_true_negatives();
    $stats{'exact_specificity'}      = $cmd->stat_exact_specificity();
    $stats{'partial_specificity'}    = $cmd->stat_partial_specificity();
    $stats{'exact_ppv'}              = $cmd->stat_exact_ppv();
    $stats{'partial_ppv'}            = $cmd->stat_partial_ppv();

    $stats{'vcf_lines_overlapping_tn'} =
      $cmd->stat_vcf_lines_overlapping_tn();

    $stats{'lines_specificity_in_tn_only'} =
      $cmd->stat_lines_specificity_in_tn_only();

    return \%stats;
}

sub stat_types {
    my $self = shift;

    my @types = (
        'true_positive_found_exact', 
        'total_true_positive_exact', 
        'sensitivity_exact', 
        'true_positive_found_partial', 
        'total_true_positive_partial', 
        'sensitivity_partial', 
        'false_positive_exact', 
        'false_positive_partial', 
        'true_negatives', 
        'exact_specificity', 
        'partial_specificity', 
        'exact_ppv', 
        'partial_ppv', 
        'vcf_lines_overlapping_tn', 
        'lines_specificity_in_tn_only', 
    );

    return @types;
}

sub create_subdirectory_workspace {
    my ($self, $inputs) = @_;
    my $base = $self->fetch_output_directory;
    my ($name, $id, $variant_type) = @{$inputs}{'name', 'id', 'variant_type'};

    my $wkspace = join('_', $name, $id);
    my $subdir = $base->subdir($wkspace)->subdir($variant_type);
    $subdir->mkpath or die $self->error_message(
        "Couldn't make workspace path '$subdir' on the file system!"
    );

    return $subdir;
}

sub validate_type {
    my ($self, $type) = @_;
    if ($type ne 'snvs' && $type ne 'indels') {
        die $self->error_message(
            "Invalid type: '$type'. Only 'snvs' and 'indels' are supported "
            . "as valid types."
        );
    }
    return 1;
}

sub fetch_config_file {
    my $self = shift;
    my $config = Path::Class::File->new($self->config_file());
    unless (-e $config) {
        die $self->error_message("Couldn't find $config on the file system!");
    }

    return $config;
}

sub fetch_output_directory {
    my $self = shift;
    my $dir = Path::Class::Dir->new($self->output_directory());
    unless (-d $dir) {
        $self->status_message(
            "'$dir' doesn't exist on file system."
            . "Making path"
        );

        $dir->mkpath or die $self->error_message(
            "Couldn't make directory $dir on the file system!"
        );
    }

    return $dir;
}

sub help_brief {
    return "Validate multiple VCF files to a specified Gold Standard Set";
}

sub help_detail {
    my $doc = q{
      This tool is a wrapper process to Genome::Model::Tools::Vcf::EvaluateVcf.
      It allows one to compare multiple SNV and Indel VCFs against a Gold Standard
      set.

      A configuration file is used to specify the list of VCFs to compare.  The
      configuration file is a TSV file with each row containing the following ordered
      columns: name, model id, vcf file system path, variant type, and sample name.

      This tools parses each line from the configruation file, appropriately
      passes the parameters into Genome::Model::Tools::Vcf::EvaluateVcf, and records
      the output statistics.  Once all the configuration lines are processed, a
      summary report is printed out to the console.

      The validation output statistics for each input VCF are also available via API
      access for custom reporting.
    };

    return $doc;
}
 
1;

__END__
