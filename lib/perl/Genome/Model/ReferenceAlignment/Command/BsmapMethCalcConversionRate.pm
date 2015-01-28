package Genome::Model::ReferenceAlignment::Command::BsmapMethCalcConversionRate;

use strict;
use warnings;
use Genome;
use Genome::Model::Tools::Bsmap::MethCalcConversionRate;

class Genome::Model::ReferenceAlignment::Command::BsmapMethCalcConversionRate {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            doc => 'Use genome build ID to calculate methylation conversion',
            shell_args_position => 1,
        },
        output_file => {
            is => 'String',
            is_optional => 1,
            doc => 'Output methylation conversion',
            shell_args_position => 2,
        },
    ],
};

sub help_synopsis {
  return <<EOS
    genome model reference-alignment bsmap-meth-calc-conversion-rate 394d6228a7b5487a9cb0ad0c448b5a44
EOS
}

sub help_brief {
    "calculate methylation conversion using Mitochondrial DNA or spiked lambda"
}

sub help_detail {
    "calculate methylation conversion using Mitochondrial DNA or spiked lambda"
}

sub execute {
    my $self = shift;
    my $output_file = $self->output_file;

    my $cfile;
    if ($output_file) {
        open($cfile, '>', $output_file) or die;
    } else {
        $cfile = \*STDOUT;
    }

    my $build = $self->build;
    my $flagstat = $build->whole_rmdup_bam_flagstat_file;
    my (@field, $total, $duplicates, $mapped, $properly);
    if (-s "$flagstat"){
        print $cfile "\nMethylation alignment status:\n";
        print $cfile $flagstat, "\n";

        my $flagstat_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat);
        print $cfile "Total read\t=\t", $flagstat_data->{total_reads}, "\n";

        print $cfile "Duplicates\t=\t", $flagstat_data->{reads_marked_duplicates}, "\n";
        if ($flagstat_data->{total_reads} != 0) {
            my $dupe_rate = $flagstat_data->{reads_marked_duplicates} / $flagstat_data->{total_reads} * 100;
            print $cfile "Duplicates rate (%)\t=\t", $dupe_rate, "\n";
        }

        print $cfile "Mapped read\t=\t", $flagstat_data->{reads_mapped}, "\n";
        if ($flagstat_data->{total_reads} != 0) {
            print $cfile "Mapped rate (%)\t=\t", $flagstat_data->{reads_mapped_percentage}, "\n";
        }

        print $cfile "Properly paired\t=\t", $flagstat_data->{reads_mapped_in_proper_pairs}, "\n";
        if ($flagstat_data->{total_reads} != 0) {
            print $cfile "Properly paired rate (%)\t=\t", $flagstat_data->{reads_mapped_in_proper_pairs_percentage}, "\n";
        }
    }
    else {
        $self->error_message("can't find flagstat file");
    }

    my $dir = $build->data_directory;
    my %cases = (
        MT => { glob => "$dir/variants/snv/meth-ratio-*/MT/snvs.hq", name => "mtDNA" },
        lambda => { glob => "$dir/variants/snv/meth-ratio-*/gi_9626243_ref_NC_001416.1_/snvs.hq", name => "lambda" },
    );
    for my $chrom (keys %cases) {
        my ($glob, $name) = ($cases{$chrom}{glob}, $cases{$chrom}{name});

        my @file = glob($glob);
        for my $file (@file) {
            if (-s "$file"){
                Genome::Model::Tools::Bsmap::MethCalcConversionRate::bs_rate($file, $chrom, $cfile);
            }
            else {
                $self->error_message("can't find $name snvs file");
            }
        }
    }
    return 1;
}

