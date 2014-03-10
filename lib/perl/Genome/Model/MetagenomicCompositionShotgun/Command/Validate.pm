package Genome::Model::MetagenomicCompositionShotgun::Command::Validate;

use strict;
use warnings;
use Genome;
use File::stat;

class Genome::Model::MetagenomicCompositionShotgun::Command::Validate {
    is => 'Genome::Model::MetagenomicCompositionShotgun::Command',
    doc => 'Validate MetagenomicCompositionShotgun build for QC and Metagenomic reports as well as headers.',
    has => [
        build_id => {
            is => 'Int',
        },
        report_dir => {
            is => 'Text',
            is_optional => 1,
        },
        verbose => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
        },
    ],
};

sub execute {
    my ($self) = @_;

    $self->dump_status_messages(1);
    $self->dump_warning_messages(1);
    $self->dump_error_messages(1);

    my $build = Genome::Model::Build->get($self->build_id);
    my $model = $build->model;

    unless ($self->report_dir){
        $self->report_dir($build->data_directory . "/reports");
    }

    my $test_bit = 0b1111; # if all tests pass result will be 1, if any fail it will be greater than 1

    # Validate BAM Headers
    $test_bit = $test_bit ^ $self->header_check();

    # Validate QC Report
    my @qc_files = ('post_trim_stats_report.tsv', 'other_stats_report.txt');
    @qc_files = map { $self->report_dir . "/$_" } @qc_files;
    $test_bit = $test_bit ^ $self->qc_check(@qc_files);

    # Validate Metagenomic Report
    my $ref_cov_report = $self->report_dir . "/metagenomic_refcov_summary.txt";
    $test_bit = $test_bit ^ $self->metagenomic_check($ref_cov_report);

    # Return
    if ($test_bit == 1) {
        $self->debug_message("Passed all checks.");
    }
    else {
        $self->debug_message("Failed to validate!");
    }
    print $self->bit_to_tests($test_bit) if ($self->verbose || $test_bit > 1);

    return $test_bit;
}

sub test_to_bit {
    my $self = shift;
    my $test = shift;

    if ($test eq 'Header') {
        return 0b0010;
    }
    if ($test eq 'QC') {
        return 0b0100;
    }
    if ($test eq 'Metagenomic') {
        return 0b1000;
    }
}

sub bit_to_tests {
    my $self = shift;
    my $bit = shift;

    my $pass = "PASSED:";
    my $fail = "FAILED:";

    ($bit & 0b0010) ? ($fail .= ' Header')             : ($pass .= ' Header'); 
    ($bit & 0b0100) ? ($fail .= ' QC_Report')          : ($pass .= ' QC_Report'); 
    ($bit & 0b1000) ? ($fail .= ' Metagenomic_Report') : ($pass .= ' Metagenomic_Report'); 

    return "$pass\n$fail\n";
}

sub header_check {
    my $self = shift;

    $self->expect64();

    my $flag = 0;

    my $meta_build = Genome::Model::Build->get($self->build_id);
    my @mga_models = $meta_build->model->metagenomic_alignment_models;
    my @data = $mga_models[0]->instrument_data;
    my $data_count = scalar(@data);
    my $combined_bam = $meta_build->_final_metagenomic_bam;
    # TODO: enable once all whole_rmdup_bams are fixed.
    my $msg = "Checking $combined_bam... ";
    my $rg_ids = `samtools view -H $combined_bam | grep \@RG | cut -f 2`; chomp $rg_ids;
    my @rg_ids = split("\n", $rg_ids);
    my $rg_count = 0;
    for my $data (@data) {
        my $data_id = $data->id;
        my @data_found = grep { $_ =~ /$data_id/ } @rg_ids;
        if (@data_found) {
            $rg_count++;
        }
    }
    if ($rg_count == $data_count) {
        $msg .= "PASS (found $rg_count read group(s))";
        $flag = $self->test_to_bit('Header');
        $self->debug_message($msg) if ($self->verbose);
    }
    else {
        $msg .= "FAIL (only found $rg_count read group(s), expected $data_count!)";
        $self->debug_message($msg);
        return $flag;
    }

    return $flag;
}

sub qc_check {
    my $self = shift;
    my @files = @_;
    my $flag = $self->test_to_bit('QC');
    for my $file (@files) {
        my $msg = "Checking for $file... ";
        if (-s $file) {
            $msg .= "PASS";
        }
        else {
            $msg .= "FAIL (empty or does not exist!)";
            $flag = 0;
        }
        $self->debug_message($msg) if ($self->verbose || $msg =~ /FAIL/);
    }
    return $flag;
}

sub metagenomic_check {
    my $self = shift;
    my $file = shift;


    my $flag = 0;
    my $msg = "Checking file ($file)... ";
    if (-s $file) {
        if (stat($file)->mtime >= 1282615095) {
            my $depth = `head -n 1 $file | cut -f 6`;
            chomp $depth;
            if ($depth eq 'Depth') {
                my $count = `cat $file | cut -f 6,7 | grep ^[0-9] | sort -u | wc -l`;
                if ($count > 1) {
                    $msg .= "PASS";
                    $flag = $self->test_to_bit('Metagenomic');
                }
                else {
                    $msg .= "FAIL (possibly corrupt)";
                }
            }
            else {
                $msg .= "FAIL (need to reparse RefCov, column 6 is '$depth' not 'Depth')";
            }
        }
        else {
            $msg .= "FAIL (generated pre-header unfix)";
        }
    }
    else {
        $msg .= "FAIL (empty file)";
    }
    $self->debug_message($msg) if ($self->verbose || $msg =~ /FAIL/);

    return $flag;
}


sub expect64 {
    my $self = shift;
    my $uname = `uname -a`;
    unless ($uname =~ /x86_64/) {
        $self->error_message("Samtools requires a 64-bit operating system.");
        die $self->error_message;
    }
}

1;
