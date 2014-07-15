package Genome::Model::Tools::SvSim::CompareBedToGold;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::SvSim::CompareBedToGold {
    is => "Command::V2",
    has_input => [
        bed_file => {
            is => "Text",
            doc => "Bed file containing predicted SV calls",
        },
        gold_file => {
            is => "Text",
            doc => "Bed file containing true SV calls",
        },
        type => {
            doc => "Type of bedfile: bed or bedpe",
            valid_values => ["bed", "bedpe"],
            default_value => "bed",
        }
    ],
    has_output => [
        stats_file => {
            is => "Text",
            doc => "Write stats (TPR, PPV, etc) to file",
            is_optional => 1,
        },
    ],
};

sub intersect_ctx {
    my ($gold, $ctx, $tp_file, $fn_file, $fp_file) = @_;
    my $tmpfile = Genome::Sys->create_temp_file_path;

    my @tp_cmd = (
        "pairToPair",
        "-a", $gold,
        "-b", $ctx,
        " | joinx sort > $tmpfile");
    my $tp_cmd = join(" ", @tp_cmd);
    Genome::Sys->shellcmd(cmd => $tp_cmd);
    dedup_by_column($tmpfile, 7, $tp_file);

    my @fn_cmd = (
        "pairToPair",
        "-type", "notboth",
        "-a", $gold,
        "-b", $ctx,
        " | joinx sort > $fn_file");
    my $fn_cmd = join(" ", @fn_cmd);
    Genome::Sys->shellcmd(cmd => $fn_cmd);

    my @fp_cmd = (
        "pairToPair",
        "-type notboth",
        "-b", $gold,
        "-a", $ctx,
        " | joinx sort > $fp_file");
    my $fp_cmd = join(" ", @fp_cmd);
    Genome::Sys->shellcmd(cmd => $fp_cmd);
}



sub intersect_files {
    my %params = @_;

    my @required = qw(gold_bed calls_bed full_intersect_bed false_neg_bed false_pos_bed);
    for my $req (@required) {
        die "Required param '$req' missing" unless exists $params{$req};
    }

    my @cmd = (
        "joinx", "intersect",
        $params{gold_bed}, $params{calls_bed},
        "--full",
        "--miss-a", $params{false_neg_bed},
        "--miss-b", $params{false_pos_bed},
        " | cut -f4- >",
        $params{full_intersect_bed},
        );

    my $cmd = join(" ", @cmd);
    return Genome::Sys->shellcmd(cmd => $cmd);
}

sub dedup_by_column {
    my ($path, $column_idx, $output_file) = @_;
    my $fh = new IO::File($path, "r");
    my $out_fh = new IO::File($output_file, "w");
    my %seen;
    while (my $line = $fh->getline) {
        chomp $line;
        my @fields = split("\t", $line);
        my $value = $fields[$column_idx];
        if (!exists $seen{$value}) {
            $seen{$value} = 1;
            $out_fh->print("$line\n");
        }

    }
    $out_fh->close();
}

sub filter_mismatched_sizes {
    my %params = @_;
    my @required = qw(max_size_diff_fraction input_bed_path output_bed_path filtered_out_bed_path);
    for my $req (@required) {
        die "Required param '$req' missing" unless exists $params{$req};
    }

    my $max_size_diff_fraction = $params{max_size_diff_fraction};
    my $input_bed_path = $params{input_bed_path};
    my $output_bed_path = $params{output_bed_path};
    my $filtered_out_bed_path = $params{filtered_out_bed_path};

    my $ifh = new IO::File($input_bed_path, "r");
    my $ofh = new IO::File($output_bed_path, "w");
    my $ofh_filt = new IO::File($filtered_out_bed_path, "w");
    while (my $line = $ifh->getline) {
        chomp $line;

        my @fnames = qw(
            chr1 start1 stop1 size1 label
            chr2 start2 stop2 size2
            );

        my @fields = split("\t", $line);
        my %values;
        @values{@fnames} = @fields;

        my $s1 = $values{size1};
        my $s2 = $values{size2};
        my $pass = 0;
        my $diff = abs($s1 - $s2);
        my $maxSize = $s1 > $s2 ? $s1 : $s2;
        my $diffPct = $diff / $maxSize;

        if ($s1 < 100 && $s2 < 100) {
            $pass = 1;
        }
        else {
            $pass = $diffPct <= $max_size_diff_fraction;
        }

        if ($pass) {
            $ofh->print("$line\n");
        }
        else {
            $ofh_filt->print("$line\n");
        }
    }
}

sub dedup_and_trim_bed {
    my ($input_path, $output_path) = @_;
    my $ifh = new IO::File($input_path, "r");
    my $ofh = new IO::File($output_path, "w");
    my %seen;
    while (my $line = $ifh->getline) {
        chomp $line;
        my @f = split(/\t/, $line);
        my $region = join(",", @f[0..2]);
        if (!exists $seen{$region}) {
            # the input is 8 columns, chr1 start1 stop1 size1 <same for ...2>
            # we are "trimming" it to only print the first four
            $ofh->print(join("\t", @f[0..4]) . "\n");
            $seen{$region} = 1;
        }
    }
}

#filtered out bed is 8 columns, we want the last 4
sub add_filtered_out_to_false_positives {
    my ($filtered_out_bed, $false_pos_bed) = @_;
    my $tmp = Genome::Sys->create_temp_file_path;
    Genome::Sys->copy_file($false_pos_bed, $tmp);

    my $ofh = new IO::File($tmp, "a");
    my $ifh = new IO::File($filtered_out_bed, "r");
    my %seen;
    while (my $line = $ifh->getline) {
        chomp $line;
        my @f = split("\t", $line);
        my $out = join("\t", @f[5..8]);
        if (!exists $seen{$out}) {
            $seen{$out} = 1;
            $ofh->print("$out\n");
        }
    }
    my $sort_cmd = "joinx sort -u $tmp -o $false_pos_bed";
    return Genome::Sys->shellcmd(cmd => $sort_cmd);
}

sub calculate_real_false_negatives {
    my ($gold_bed, $true_pos_bed, $false_neg_bed) = @_;
    my @cmd = (
        "joinx intersect",
        $gold_bed,
        $true_pos_bed,
        "-o /dev/null",
        "--miss-a $false_neg_bed");

    my $cmd = join(" ", @cmd);
    return Genome::Sys->shellcmd(cmd => $cmd);
}

sub get_stats {
    my $prefix = shift;

    my $tp = `wc -l ${prefix}-tp.bed* | awk '{print \$1}'`;
    my $fp = `wc -l ${prefix}-fp.bed* | awk '{print \$1}'`;
    my $fn = `wc -l ${prefix}-fn.bed* | awk '{print \$1}'`;
    chomp $tp;
    chomp $fp;
    chomp $fn;

    my %stats;
    $stats{tp_count} = $tp;
    $stats{fp_count} = $fp;
    $stats{fn_count} = $fn;
    if ($tp > 0) {
        $stats{sensitivity} = $tp / ($tp + $fn);
        $stats{ppv} = $tp / ($tp + $fp);
        $stats{fdr} = $fp / ($tp + $fp);
    }
    else {
        $stats{ppv} = 0;
        $stats{sensitivity} = 0;
        if ($fp > 0) {
            $stats{fdr} = 1;
        }
        else {
            $stats{fdr} = "NA";
        }
    }

    $stats{total} = $tp + $fn;

    for my $k (sort keys %stats) {
        print "$k: $stats{$k}\n";
    }
    return %stats;
}

sub do_bedpe {
    my $self = shift;
    my $ctx_bedpe = $self->bed_file;
    my $gold_ctx = $self->gold_file;

    my @basename = $ctx_bedpe =~ /(.*)\.[^.]*$/;
    my $basename = $basename[0];

    my $ctx_bedpe_tp = sprintf("%s-tp.bedpe", $basename);
    my $ctx_bedpe_fn = sprintf("%s-fn.bedpe", $basename);
    my $ctx_bedpe_fp = sprintf("%s-fp.bedpe", $basename);

    intersect_ctx($gold_ctx, $ctx_bedpe,
        $ctx_bedpe_tp,
        $ctx_bedpe_fn,
        $ctx_bedpe_fp);

    return 1;
}

sub execute {
    my $self = shift;

    if ($self->type eq 'bedpe') {
        return $self->do_bedpe;
    }

    my @detectors = qw(bd pd sd all);
    my @sv_types = qw(ins del inv);
    my $max_size_diff_fraction = 0.5;

    my $calls_bed = $self->bed_file;
    my $gold_bed = $self->gold_file;
    my @basename = $calls_bed =~ /(.*)\.[^.]*$/;
    my $basename = $basename[0];

    my $full_intersect_bed = sprintf("%s-tp-full.bed", $basename);
    my $false_pos_bed = sprintf("%s-fp.bed", $basename);
    my $false_neg_bed = sprintf("%s-fn.bed", $basename);

    intersect_files(
        gold_bed => $gold_bed,
        calls_bed => $calls_bed,
        full_intersect_bed => $full_intersect_bed,
        false_pos_bed => $false_pos_bed,
        false_neg_bed => $false_neg_bed,
        );

    my $true_pos_dup_bed = sprintf("%s-tp-dups.bed", $basename);
    my $true_pos_filtered_out_bed = sprintf("%s-tp-filtered-out.bed", $basename);
    my $true_pos_bed = sprintf("%s-tp.bed", $basename);
    filter_mismatched_sizes(
        max_size_diff_fraction => $max_size_diff_fraction,
        input_bed_path => $full_intersect_bed,
        output_bed_path => $true_pos_dup_bed,
        filtered_out_bed_path => $true_pos_filtered_out_bed,
        );

    add_filtered_out_to_false_positives($true_pos_filtered_out_bed, $false_pos_bed);

    dedup_and_trim_bed($true_pos_dup_bed, $true_pos_bed);

    calculate_real_false_negatives($gold_bed, $true_pos_bed, $false_neg_bed);

    my %stats = get_stats($basename);
    my $stats_file = $self->stats_file;
    if ($stats_file) {
        my $fh = Genome::Sys->open_file_for_writing($stats_file);
        $fh->print(join("\n", map {sprintf("%s: %s", $_, $stats{$_})} keys %stats) . "\n");
    }

    return 1;
}

1;
