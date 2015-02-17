package Genome::Model::Tools::Joinx::VcfReport;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use Carp qw/confess/;
use POSIX qw( WIFEXITED );

our $MINIMUM_JOINX_VERSION = 1.5;

class Genome::Model::Tools::Joinx::VcfReport {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Vcf File to filter',
            shell_args_position => 1,
        },
    ],
    has_optional_input => [
        per_sample_output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'The per sample report output file',
        },
        per_site_output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'The per site report output file',
        },
        info_fields_from_db => {
            is => 'Text',
            doc => 'Field ids that define whether an allele is novel. Use colons to separate multiple field descriptors.',
            #doing the above because UR autosplits on commas with is_many, but joinx uses commas in its field descriptors
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'zcats the input file into stdin',
            default => 0,
        },
        generate_report_only => {
            is => 'Boolean',
            doc => "don't run joinx, but graph a report based on the passed output file names",
            default => 0,
        },
        generate_report => {
            is => 'Boolean',
            doc => "graph a report based on the output files",
            default => 1,
        },
        report_image_type => {
            is => 'Text',
            doc => 'type of image to graph',
            default => 'pdf',
            valid_values => [ 'pdf', 'png' ],
        },

    ],
};

sub help_brief {
    "Generates some reports on the variants and samples in a VCF file"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx vcf-report --input-file a.vcf --info-fields-from-db dbSNPBuildID
EOS
}

sub execute {
    my $self = shift;

    $self->check_minimum_version($MINIMUM_JOINX_VERSION);

    my ($per_sample_output, $per_site_output) = qw{ per_sample_report.txt per_site_report.txt }; #these are the default output file names in the program

    my $input_file = $self->input_file;

    unless(-s $input_file) {
        die $self->error_message("$input_file does not exist");
    }

    if($self->use_bgzip){
        $input_file = "<(zcat $input_file)";
    }

    unless($self->generate_report_only) {
        my $cmd = $self->joinx_path . " vcf-report" . " --input-file $input_file";
        my $info_fields = " --info-fields-from-db " . join(" --info-fields-from-db ", split /:/, $self->info_fields_from_db) if($self->info_fields_from_db);
        $cmd .= $info_fields if $info_fields;

        if(defined($self->per_site_output_file)) {
            $per_site_output = $self->per_site_output_file;
            $cmd .= " --per-site-file $per_site_output";
        }
        if(defined($self->per_sample_output_file)) {
            $per_sample_output = $self->per_sample_output_file;
            $cmd .= " --per-sample-file $per_sample_output";
        }
        if($self->use_bgzip) {
            $cmd = qq{bash -c "$cmd"};
        }

        my %params = (
            cmd => $cmd,
            output_files => [$per_sample_output, $per_site_output],
            skip_if_output_is_present => 0,
        );
        Genome::Sys->shellcmd(%params);
    }

    if($self->generate_report || $self->generate_report_only) {
        my $rscript = $self->rlib_path . "/VcfReport.R";
        unless(-e $rscript) {
            $self->error_message("Unable to find $rscript for generating report");
            return;
        }
        my $image_type = $self->report_image_type;
        my $cmd = qq{echo 'source("$rscript");write_reports("$per_site_output","$per_sample_output",image_type="$image_type")'};
        $cmd .= "| R --slave";
        print "R:\n$cmd\n";
        WIFEXITED(system $cmd) or confess "Couldn't run: $cmd ($?)";
    }

    return 1;
}

1;
