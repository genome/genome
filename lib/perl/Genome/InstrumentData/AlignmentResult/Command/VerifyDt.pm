package Genome::InstrumentData::AlignmentResult::Command::VerifyDt;

use strict;
use warnings;
use Genome;
use DateTime;

class Genome::InstrumentData::AlignmentResult::Command::VerifyDt {
    is => 'Genome::Command::Base',
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'The instrument data to locate alignment results (and merged BAMs).',
            shell_args_position => 1,
        },
        max => {
            is => 'Integer',
            is_optional => 1,
            doc => 'Maximum number of instrument data to verify.',
        },
        repair => {
            is => 'Boolean',
            default => 0,
            doc => 'Attempt to repair BAMs.',
        }
    ],
    doc => 'Verify that the timestamps in the BAM file are ISO 8601 formatted.',
};

sub help_synopsis {
    my $class = shift;
    return <<EOS;
genome instrument-data alignment-result verify-timestamps 123456,456789
genome instrument-data alignment-result verify-timestamps --instrument-data=123456,456789
EOS
}

sub help_detail {
    my $class = shift;
    return <<'EOS';
Verify that the timestamps in the BAM file are ISO 8601 formatted.
EOS
}

sub execute {
    my $self = shift;

    my @instrument_data = $self->instrument_data;
    if ($self->max) {
        my $max;
        if ($self->max < @instrument_data) {
            $max = $#instrument_data 
        } else {
            $max = $self->max - 1;
        }
        @instrument_data = @instrument_data[0..$max];
    }

    my @instrument_data_ids = map { $_->id } $self->instrument_data;

    if ($self->max) {
        $self->status_message("Reduced number of instrument data to " . @instrument_data_ids . ": " . join(", ", @instrument_data_ids) . ".");
    }

    my @bams;
    for my $instrument_data_id (@instrument_data_ids) {
        my @alignment_results = Genome::InstrumentData::AlignmentResult->get(
            instrument_data_id => $instrument_data_id
        );
        unless (@alignment_results) {
            $self->status_message("No alignment results found for instrument data (ID: $instrument_data_id).");
        }
        for my $alignment_result (@alignment_results) {
            push @bams, $alignment_result->output_dir . "/all_sequences.bam";
            my @builds = map { Genome::Model::Build->get($_->user_id) } $alignment_result->users;
            for my $build (@builds) {
                my @merged_bams = glob($build->data_directory . "/alignments/*_merged_rmdup.bam");
                @merged_bams = grep { $_ =~ /\/\d+_merged_rmdup.bam/ } @merged_bams;
                push @bams, @merged_bams;
            }
        }
    }

    for my $in_bam (@bams) {
        if (-f $in_bam && -s $in_bam) {
            unless(valid_dt_tag($in_bam)) {
                repair_dt($in_bam) if ($self->repair);
            }
        } else {
            print "BAM ($in_bam) does not appear to be a valid file:\n";
            system("ls $in_bam");
        }
    }
    return 1;
}

sub valid_dt_tag {
    my $file = shift;
    chomp(my @dt_lines = qx(samtools view -H $file | grep DT:));
    my $valid = 1;
    for my $dt_line (@dt_lines) {
        my @tags = split("\t", $dt_line);
        my ($dt_tag) = grep { $_ =~ /^DT:/ } @tags;
        if ($dt_tag =~ /^DT:\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}[+-]\d{4}$/
            || $dt_tag =~ /^DT:\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$/) {
            print "Valid: ($dt_tag) $file\n";
        } elsif ($dt_tag) {
            print "Invalid: ($dt_tag) $file\n";
            $valid = 0;
        } else {
            print "Missing DT tag: $file\n";
        }
    }
    return $valid;
}

sub reheader {
    my $header_file = shift || die;
    my $original_bam = shift || die;
    my $fixed_bam = shift || die;

    my $input_header = IO::File->new($header_file);
    my $input_reads  = IO::File->new("samtools view $original_bam |");
    my $output_bam   = IO::File->new("| samtools view -o $fixed_bam -bS -");

    while (my $line = $input_header->getline) {
        $output_bam->print($line);
    }
    while (my $line = $input_reads->getline) {
        $output_bam->print($line);
    }

    return 1;
}

sub repair_dt {
    my $in_bam = shift;
    my $in_bam_lock = "$in_bam.verify-dt.lock";
    if (-e $in_bam_lock) {
        die "Already running a Verify-DT on bam ($in_bam), lock file exists.\n";
    }
    system("touch $in_bam_lock");
    (my $in_sam_h = $in_bam) =~ s/\.bam$/.sam.h/;
    (my $out_bam = $in_bam) =~ s/\.bam$/_fixed.bam/;
    (my $out_sam_h = $in_sam_h) =~ s/\.sam\.h$/_fixed.sam.h/;
    system("samtools view -H $in_bam > $in_sam_h") && die;
    system("cp $in_sam_h $out_sam_h") && die;
    chomp(my @dt_lines = qx(samtools view -H $in_bam | grep DT:));
    print "Repairing headers...\n";
    for my $dt_line (@dt_lines) {
        my @tags = split("\t", $dt_line);
        my ($dt_tag) = grep { $_ =~ /^DT:/ } @tags;
        unless ($dt_tag) {
            print "Could not find DT tag: $dt_line.\n";
            next;
        }
        my ($year, $month, $day, $hour, $min, $sec) = $dt_tag =~ /DT:(\d{4})-(\d{2})-(\d{2})\ (\d{2}):(\d{2}):(\d{2})/;
        unless ($sec) {
            print "Unmatched DT tag: $dt_tag.\n";
            next;
        }
        my $datetime = DateTime->new(
            year      => $year,
            month     => $month,
            day       => $day,
            hour      => $hour,
            minute    => $min,
            second    => $sec,
            time_zone => 'America/Chicago',
        );
        $datetime->set_time_zone('UTC');
        my $new_dt_tag = "DT:${datetime}Z";
        print "\tChanging DT tag from $dt_tag to $new_dt_tag.\n";
        $dt_tag =~ s/\ /\\ /;
        system("sed -i 's/$dt_tag/$new_dt_tag/' $out_sam_h") && die;
    }
    print "Repairing bam...\n";
    #system("samtools reheader $out_sam_h $in_bam > $out_bam") && die;
    reheader($out_sam_h, $in_bam, $out_bam) || die;

    print "Validating new bam...\n";
    my $in_bam_h_md5 = qx(samtools view -H $in_bam | grep -v DT: | md5sum);
    my $out_bam_h_md5 = qx(samtools view -H $out_bam | grep -v DT: | md5sum);
    unless ($in_bam_h_md5 eq $out_bam_h_md5 && length $in_bam_h_md5 > 33) {
        die "ERROR: BAM headers (without DT lines) do not match between $in_bam and $out_bam.\n";
    }
    my $in_bam_md5 = qx(samtools view $in_bam | md5sum);
    my $out_bam_md5 = qx(samtools view $out_bam | md5sum);
    unless ($in_bam_md5 eq $out_bam_md5 && length $in_bam_md5 > 33) {
        die "\tERROR: BAM contents do not match between $in_bam and $out_bam.\n";
    }
    rename($in_bam, "$in_bam.orig") || die;
    rename($out_bam, $in_bam) || die;

    if (-e "$in_bam.md5") {
        unlink("$in_bam.md5") || die;
    }
    print "Regenerating the MD5...\n";
    !system("md5sum $in_bam > $in_bam.md5") || die;

    if (-e "$in_bam.bai") {
        unlink("$in_bam.bai") || die;
    }
    print "Regenerating the BAM index ...\n";
    !system("samtools index $in_bam") || die;

    unlink("$in_bam.orig") || die;
    unlink($out_sam_h) || die;
    unlink($in_sam_h) || die;
    unlink($in_bam_lock) || die;
    print "\n";
}

