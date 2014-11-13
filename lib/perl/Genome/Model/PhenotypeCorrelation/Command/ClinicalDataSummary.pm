package Genome::Model::PhenotypeCorrelation::Command::ClinicalDataSummary;

use Genome;
use Data::Dumper;
use Sort::Naturally qw/nsort/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::ClinicalDataSummary {
    is => "Command::V2",
    doc => "Produce summary tables and plots for clinical data files",
    has_input => [
        glm_model_file => {
            is => "File",
            doc => "The GLM model file that describes clinical traits",
        },

    # clinical data file and samples/nomenclature are mutually exclusive
        clinical_data_file => {
            is => "File",
            doc => "The clinical data file to process. Specify EITHER this OR samples/nomenclature.",
            is_optional => 1,
        },

        samples => {
            is => "Genome::Sample",
            doc => "List of samples (or a population group) by name or id",
            is_many => 1,
            is_optional => 1,
        },
        nomenclature => {
            is => "Genome::Nomenclature",
            doc => "nomenclature used to access clinical data",
            is_optional => 1,
        },
    ######################################################################

        delim => {
            is => "Text",
            doc => "The field delimiter used in the clinical and glm files",
            default => "\t",
        },
        output_delim => {
            is => "Text",
            doc => "Delimiter to use in output files",
            default_value => "\t",
        },
        output_directory => {
            is => "File",
            doc => "Output directory",
            is_output => 1,
        },
        include_boxplots => {
            is => "Boolean",
            doc => "Include box and whisker plots for quantitative traits",
            default_value => 0,
        },
    ],
    has_transient_optional => [
        _clinical_data => {
            is => "Genome::Model::PhenotypeCorrelation::ClinicalData",
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
genome model phenotype-correlation clinical-data-summary \
    --glm-model-file glm-model.txt \
    --clinical-data-file clindat.csv \
    --tables-file tab.txt \
EOS
}

sub help_detail {
    return <<"EOS"
Summarize clinical data with text tables (categorical traits) and
histograms/box plots (quantitative traits).
EOS
}

sub execute {
    my $self = shift;

    Genome::Sys->create_directory($self->output_directory);

    my $tables_file = $self->output_directory . "/categorical.txt";
    my $plot_output_dir = $self->output_directory . "/plots";

    if (defined $self->clinical_data_file) {
        die "Please specify either --clinical-data-file OR --samples and --nomenclature"
            unless !defined $self->samples and !defined $self->nomenclature;
        $self->_clinical_data(Genome::Model::PhenotypeCorrelation::ClinicalData->from_file($self->clinical_data_file));
    } else {
        die "Please specify either --clinical-data-file OR --samples and --nomenclature"
            unless defined $self->samples and defined $self->nomenclature;
        my @samples = $self->samples;
        $self->_clinical_data(Genome::Model::PhenotypeCorrelation::ClinicalData->from_database($self->nomenclature, @samples));
    }

    my $cd = $self->_clinical_data;
    my $glm_config = Genome::Model::PhenotypeCorrelation::GlmConfig->from_file($self->glm_model_file);

    my %attr_types = %{$cd->attribute_types};
    my @categorical_attrs = $glm_config->categorical_attributes;
    my @quantitative_attrs = $glm_config->quantitative_attributes;

    my %tables;
    my %categorical_with_quant_covars;

    if (@categorical_attrs) {
        my $tab_fh = Genome::Sys->open_file_for_overwriting($tables_file);
        for my $attr (@categorical_attrs) {
            my ($table, @quant_covars) = $self->_make_table_data($attr);
            if (@quant_covars) {
                $categorical_with_quant_covars{$attr->{attr_name}} = []
                    unless defined $categorical_with_quant_covars{$attr->{attr_name}};
                push(@{$categorical_with_quant_covars{$attr->{attr_name}}}, @quant_covars);
            }
            if (scalar(%$table)) {
                $self->_write_tables($attr, $table, $tab_fh);
            }

            $tables{$attr->{attr_name}} = { attr => $attr, data => $table };
        }
    } else {
        $self->status_message("No categorical traits found, skipping table generation");
    }

    if (@quantitative_attrs || %categorical_with_quant_covars) {
        $self->_make_plots($plot_output_dir, \@quantitative_attrs, \%categorical_with_quant_covars);
    } else {
        $self->status_message("No quantitative traits found, skipping plot generation");
    }

    return 1;
}

sub _parse_attr {
    my ($self, $line) = @_;
    my ($type, $name, $gene, $covar, $memo) = split($self->delim, $line);
    if ($type ne "B" and $type ne "Q") {
        die "Unknown attribute type $type, expected 'B' or 'Q'";
    }

    my @covariates = grep { $_ ne "NA" } split(' *\+ *', $covar);

    return {
        type => $type,
        name => $name,
        covariates => \@covariates,
    };
}

sub _make_table_data {
    my ($self, $attr) = @_;

    my $cd = $self->_clinical_data;

    my $tables = {};
    my @quant_covars;
    my $attr_values = $cd->attribute_values($attr->{attr_name});
    my %possible_attr_values = map { $_ => 1 } @$attr_values;
    my $self_table = {};
    for my $value (keys %possible_attr_values) {
        $self_table->{$attr->{attr_name}}->{$value} = scalar grep { $value eq $_ } @{$attr_values};
    }
    $tables->{$attr->{attr_name}} = $self_table;

    for my $covar (@{$attr->{covariates}}) {

        my $covar_type = $cd->attribute_types->{$covar}
            or die "Failed to find type information for covariate $covar";

        if (!$covar_type->{categorical}) {
            push(@quant_covars, $covar);
            next;
        }

        my $data = {};
        my $covar_values = $cd->attribute_values($covar);
        for my $i (0..$cd->sample_count - 1) {
            my $covar_value = $covar_values->[$i];
            my $attr_value = $attr_values->[$i];

            $data->{$covar_value} = {} unless defined $data->{$covar_value};
            if (defined $data->{$covar_value}->{$attr_value}) {
                ++$data->{$covar_value}->{$attr_value};
            } else {
                $data->{$covar_value}->{$attr_value} = 1;
            }
        }
        $tables->{$covar} = $data;
    }

    return $tables, @quant_covars;
}

sub _make_plots {
    my ($self, $plot_output_dir, $quant_attrs, $categorical_with_quant_covars) = @_;

    my $cd = $self->_clinical_data;

    my @R_code;
    my @paths;
    for my $attr (@$quant_attrs) {
        my @data = @{$cd->attribute_values($attr)};
        next unless @data;

        my $file_suffix = Genome::Utility::Text::sanitize_string_for_filesystem("$attr.png");
        my $hist_path = join("/", $plot_output_dir, "histogram-$file_suffix");
        my $boxplot_path = join("/", $plot_output_dir, "boxplot-$file_suffix");
        push(@paths, $hist_path, $boxplot_path);
        push(@R_code, $self->_R_code_for_plot("hist", $hist_path, $attr, $attr, \@data));
        push(@R_code, $self->_R_code_for_plot("boxplot", $boxplot_path, $attr, $attr, \@data)); 
    }

    # Process quantitative covariates of categorical attributes
    # We split them into categories for each value the categorical attribute can take on.
    for my $attr (keys %$categorical_with_quant_covars) {
        my @attr_values = @{$cd->attribute_values($attr)};
        my %uniq = map {defined $_ ? ($_ => 1) : ("NA" => 1)} @attr_values;
        my @uniq_attr_values = keys %uniq;

        my @covars = @{$categorical_with_quant_covars->{$attr}};
        for my $covar (@covars) {
            my @covar_values = @{$cd->attribute_values($covar)};
            for my $attr_value (@uniq_attr_values) {
                next if $attr_value eq "NA";
                my @indices = grep { "$attr_values[$_]" eq "$attr_value"} 0..$#attr_values;
                my @data = grep {defined $_} @covar_values[@indices];
                next unless @data;

                my $file_suffix = Genome::Utility::Text::sanitize_string_for_filesystem("$covar-$attr=$attr_value.png");
                my $hist_path = join("/", $plot_output_dir, "histogram-$file_suffix");
                my $boxplot_path = join("/", $plot_output_dir, "boxplot-$file_suffix");
 
                push(@paths, $hist_path, $boxplot_path);

                my $title = "$covar, $attr=$attr_value";
                push(@R_code, $self->_R_code_for_plot("hist", $hist_path, $title, $attr, \@data));
                push(@R_code, $self->_R_code_for_plot("boxplot", $boxplot_path, $title, $attr, \@data));
            }
        }
    }

    Genome::Sys->create_directory($plot_output_dir);
    my ($tmpfh, $tmpfile) = Genome::Sys->create_temp_file;
    $tmpfh->print(join("\n", @R_code));
    $tmpfh->close;
	my $cmd = "R --vanilla --slave \< $tmpfile";

	my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
        output_files => \@paths,
        skip_if_output_is_present => 0,
    );
	unless($return) { 
		$self->error_message("Failed to execute R command $cmd.");
		die $self->error_message;
    }

    $self->status_message("Wrote " . scalar(@paths) . " plots\n");
}

sub _R_code_for_plot {
    my ($self, $type, $path, $title, $attr_name, $data) =@_;

    die "Unknown plot type '$type', expected 'hist' or 'boxplot'"
        unless $type eq 'hist' or $type eq 'boxplot';

    my @valid_obs = grep {$_} @$data;
    my $R_data = "c(" . join(", ", @valid_obs) . ")";
    my $R_code = qq{png("$path"); $type($R_data, main="$title", xlab="$attr_name"); devoff <- dev.off();};
    return $R_code;
}

sub _write_tables {
    my ($self, $attr, $table_data, $fh) = @_;

    my $cd = $self->_clinical_data;

    my $delim = $self->output_delim;
    my $attr_name = $attr->{attr_name};
    my %attr_values = map {$_ => 1} @{$cd->attribute_values($attr_name)};
    my @attr_values = nsort keys %attr_values;

    $fh->print("$attr_name\t" . join("$delim", @attr_values) . "\n");
    for my $covar (keys %$table_data) {
        next if $covar eq $attr_name;
        $fh->print("-"x30 . "\n");
        $fh->print("$covar\n");
        $fh->print("-"x30 . "\n");
        my $d = $table_data->{$covar};
        for my $covar_value (keys %$d) {
            $fh->print("$covar_value");
            for my $tv (@attr_values) {
                my $count = $d->{$covar_value}->{$tv} || 0;
                $fh->print("$delim$count");
            }
            $fh->print("\n");
        }
    }
    $fh->print("-"x30 . "\n");
    $fh->print("Total\t");
    $fh->print(join("\t", map {$table_data->{$attr_name}->{$attr_name}->{$_}} @attr_values) . "\n");
}

1;
