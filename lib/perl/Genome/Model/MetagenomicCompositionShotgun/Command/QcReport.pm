package Genome::Model::MetagenomicCompositionShotgun::Command::QcReport;

use strict;
use warnings;
use Genome;
use File::Path;
use File::Find;

$|=1;

class Genome::Model::MetagenomicCompositionShotgun::Command::QcReport{
    is => 'Genome::Model::MetagenomicCompositionShotgun::Command',
    doc => 'Generate QC report for a MetagenomicCompositionShotgun build.',
    has => [
        build_id => {
            is => 'Int',
        },
        overwrite => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
        },
        report_dir => {
            is => 'Text',
            is_optional => 1,
        },
        skip_qc_on_untrimmed_reads => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => "If this flag is specified, QC report will skip metric on the human-free, untrimmed data.",
        },
    ],
};

sub execute {
    my $self = shift;

    # misc. build/model data 
    my $mcs_build = Genome::Model::Build->get($self->build_id);
    my $mcs_model = $mcs_build->model;
    my @meta_models = $mcs_model->_metagenomic_alignment_models;
    my @original_data = $mcs_build->instrument_data;
    my @hcs_data = $meta_models[0]->instrument_data;
    my $dir = $mcs_build->data_directory;
    my ($contamination_bam, $contamination_flagstat, $meta1_bam, $meta1_flagstat, $meta2_bam, $meta2_flagstat) = map{ $dir ."/$_"}(
        "contamination_screen.bam",
        "contamination_screen.bam.flagstat",
        "metagenomic_alignment1.bam",
        "metagenomic_alignment1.bam.flagstat",
        "metagenomic_alignment2.bam",
        "metagenomic_alignment2.bam.flagstat",
    );

    # setup report directory path
    unless ($self->report_dir) {
        $self->report_dir($mcs_build->data_directory . '/reports');
    }
    unless (-d $self->report_dir) {
        mkpath $self->report_dir;
    }
    $self->debug_message('Report Path: ' . $self->report_dir);

    # get existing metrics
    my @mcs_metrics = $mcs_build->metrics;
    my %mcs_metrics;
    for my $mcs_metric (@mcs_metrics) {
        $mcs_metrics{$mcs_metric->name} = $mcs_metric;
    }

    ### POST TRIMMING STATS ###
    $self->post_trimming_stats($contamination_bam, \%mcs_metrics);

    ### PER LANE QC ###
    my @hcs_paired_data_ids = map { $_->id } grep { $_->is_paired_end } @hcs_data;
    my @imported_data_ids = @hcs_paired_data_ids; # replace later
    for my $hcs_data (@hcs_data) {
        unless ($self->skip_qc_on_untrimmed_reads) {
            $self->per_lane_qc($hcs_data, \%mcs_metrics);
        }
    }

    ### OTHER STATS ###
    my $other_stats_output_path = $self->report_dir . '/other_stats_report.txt';
    if (-f $other_stats_output_path) {
        unlink($other_stats_output_path);
    }
    my $other_stats_output = Genome::Sys->open_file_for_writing($other_stats_output_path);
    for my $hcs_data (@hcs_data) {
        $self->other_stats(\%mcs_metrics, $other_stats_output, $hcs_data);
    }

    print "\n\n";
    $self->debug_message('Model: ' . $mcs_model->name);
    $self->debug_message('Build: ' . $mcs_build->id);
    $self->debug_message('Reports: ' . $self->report_dir);
    $self->debug_message('Report Data: ' . $self->report_dir . '/data');
    my $done = system("touch ".$self->report_dir."/FINISHED");
    return 1;
}

sub post_trimming_stats {
    # number of quality trimmed bases per lane
    # average length of quality trimmed bases per lane
    my ($self, $contamination_bam, $mcs_metrics_hash) = @_;
    my %mcs_metrics = %$mcs_metrics_hash;

    my $stats_output_path = $self->report_dir . '/post_trim_stats_report.tsv';
    $self->debug_message("Generating post trimming stats...");
    unlink($stats_output_path) if (-f $stats_output_path);
    my $stats_output = Genome::Sys->open_file_for_writing($stats_output_path);

    $self->debug_message("\tParsing $contamination_bam...");
    my %stats = $self->bam_stats_per_lane($contamination_bam);
    $self->debug_message("\tParsing of $contamination_bam is complete.");

    # print header
    print $stats_output "flow_lane";
    my $flow_lane = (keys(%stats))[0];
    for my $stat (sort(keys %{$stats{$flow_lane}})) {
        print $stats_output "\t$stat";
    }
    print $stats_output "\n";

    # print values
    my $metric_name;
    for my $flow_lane (sort(keys %stats)) {
        print $stats_output $flow_lane;
        for my $stat (sort(keys %{$stats{$flow_lane}})) {
            $metric_name = $flow_lane . "_" . $stat;
            $mcs_metrics{$metric_name}->delete() if($mcs_metrics{$metric_name});
            unless(Genome::Model::Metric->create(build_id => $self->build_id, name => $metric_name, value => $stats{$flow_lane}{$stat})) {
                die $self->error_message("Unable to create build metric (build_id=" . $self->build_id . ", $metric_name)");
            }
            print $stats_output "\t" . $stats{$flow_lane}{$stat};
        }
        print $stats_output "\n";
    }
}

sub per_lane_qc {
    my ($self, $hcs_data, $mcs_metrics_hash) = @_;
    my %mcs_metrics = %$mcs_metrics_hash;
    my $imported_data = $hcs_data; # replace later
    my $temp_dir = Genome::Sys->base_temp_directory;

    my %fastq_files;
    my @imported_fastq;
    my @original_fastq;
    my $data_path = $self->report_dir . '/data';
    mkpath($data_path) unless (-d $data_path);

    my $hcs_data_id = $imported_data->id;
    my $original_data = $self->original_data_from_imported_id($hcs_data_id);

    my $humanfree_bam_path = $self->report_dir . '/data/' . $hcs_data_id . '_humanfree_untrimmed.bam';
    my $original_bam_path = $self->report_dir . '/data/' . $hcs_data_id . '_original_untrimmed.bam';
    my $humanfree_fwd_path = $temp_dir . '/' . $hcs_data_id . '_1_humanfree_untrimmed';
    my $humanfree_rev_path = $temp_dir . '/' . $hcs_data_id . '_2_humanfree_untrimmed';

    ### ###

    $self->debug_message("Extracting FastQ files from original and imported data...");    
    if (! $self->overwrite && -f $humanfree_bam_path && -f $original_bam_path) {
        $self->debug_message("\t$hcs_data_id: Skipping FastQ extraction for " . $hcs_data_id . ", bam files already exists. Use --overwrite to replace...");
    }
    else {
        # untar both imported and original fastq files, only keeping paired files
        if ($imported_data->is_paired_end) {
            my @imported_fastq_filenames = $imported_data->dump_sanger_fastq_files;
            my @original_fastq_filenames;
            @original_fastq_filenames = $original_data->dump_sanger_fastq_files;
            for my $file (@imported_fastq_filenames) {
                my $name = (split('/', $file))[-1];
                $name =~ s/\.txt$//;
                my $output_filename = $temp_dir . '/' . $name . '_imported_trimmed';
                die "failed to rename $file to $output_filename" unless rename($file, $output_filename);
                push @imported_fastq, $output_filename;
                push @{$fastq_files{$hcs_data_id}{imported}}, $output_filename; 
            }
            for my $file (@original_fastq_filenames) {
                my $name = (split('/', $file))[-1];
                $name =~ s/\.txt$//;
                my $output_filename = $temp_dir . '/' . $hcs_data_id . '_' . $name . '_original';
                die "failed to rename $file to $output_filename" unless rename($file, $output_filename);
                push @original_fastq, $output_filename;
                push @{$fastq_files{$hcs_data_id}{original}}, $output_filename;
            }
        }
        else {
            $self->debug_message("\tSkipping unpaired fastq (imported id: " . $imported_data->id . ")...");
            return;
        }
    }

    ### ###
    $self->debug_message("Generating human-free, untrimmed data...");
    if (! $self->overwrite && -f $humanfree_bam_path) {
        $self->debug_message("\t$hcs_data_id: Skipping humanfree creation, humanfree bam file already exists. Use --overwrite to replace...");
    }
    else {
        $self->debug_message("\tGenerating hash of read names for instrument data: $hcs_data_id...");

        my $imported_path = (@{$fastq_files{$hcs_data_id}{imported}})[0];
        my $original_fwd_path = @{$fastq_files{$hcs_data_id}{original}}[0];
        my $original_rev_path = @{$fastq_files{$hcs_data_id}{original}}[1];
        my $imported_file = Genome::Sys->open_file_for_reading($imported_path);
        my $humanfree_fwd_file = Genome::Sys->open_file_for_writing($humanfree_fwd_path);
        my $humanfree_rev_file = Genome::Sys->open_file_for_writing($humanfree_rev_path);

        $self->debug_message("\t\tReading in up to 8M read names...");
        my $reads_left = 1;
        my $readname_re = '[^:]*:(.*)#.*';
        while ($reads_left) {
            my %read_names;
            # read 8M reads at a time to prevent oom
            for (my $count = 0; $count < 8e6; $count++) {
                $self->debug_message("\t\t\tHashed $count reads...") unless ($count % 2e6 || ! $count);
                my $imported_read = read_and_join_lines($imported_file);
                unless ($imported_read) {
                    $self->debug_message("\t\t\tFinished reading $count reads.");
                    $reads_left = 0;
                    last;
                }
                $imported_read =~ /$readname_re/;
                my $imported_readname = $1;
                $read_names{$imported_readname} = 1;
            }

            my $original_fwd_file = Genome::Sys->open_file_for_reading($original_fwd_path);
            my $original_rev_file = Genome::Sys->open_file_for_reading($original_rev_path);

            $self->debug_message("\t\tParsing original forward read file with those hashed reads...");
            while (my $fwd_read = read_and_join_lines($original_fwd_file)) {
                $fwd_read =~ /$readname_re/;
                my $fwd_readname = $1;
                if ($read_names{$fwd_readname}) {
                    print $humanfree_fwd_file $fwd_read ;
                }
            }

            $self->debug_message("\t\tParsing original reverse read file with those hashed reads...");
            while (my $rev_read = read_and_join_lines($original_rev_file)) {
                $rev_read =~ /$readname_re/;
                my $rev_readname = $1;
                if ($read_names{$rev_readname}) {
                    print $humanfree_rev_file $rev_read ;
                }
            }
        }
    }

    ### ###
    $self->debug_message("Creating bams...");

    if (! $self->overwrite && -f $humanfree_bam_path ) {
        $self->debug_message("\t$hcs_data_id: Skipping humanfree bam creation, files already exists. Use --overwrite to replace...");
    }
    else {
        unlink($humanfree_bam_path);

        $self->debug_message("\t$hcs_data_id: Generating humanfree bam file...");

        $self->debug_message("\tVerifying fwd/rev pairs are correct, will swap if not...");

        my $humanfree_fwd_file = Genome::Sys->open_file_for_reading($humanfree_fwd_path);
        my $humanfree_rev_file = Genome::Sys->open_file_for_reading($humanfree_rev_path);
        # If humanfree_untrimmed files are reversed then file contents are probably switched so switch files.
        my $humanfree_fwd_line = $humanfree_fwd_file->getline;
        if ($humanfree_fwd_line =~ /\/2$/) {
            $self->debug_message("\t\t" . (split("/", $humanfree_fwd_path))[-1] . " looks like a reverse file. Swapping...");
            rename($humanfree_fwd_path, $humanfree_fwd_path . '.tmp');
            rename($humanfree_rev_path, $humanfree_fwd_path);
            rename($humanfree_fwd_path . '.tmp', $humanfree_rev_path);
        }

        $self->debug_message("\t" . (split("/", $humanfree_bam_path))[-1]. "...");
        my $humanfree_fastq_to_sam = Genome::Model::Tools::Picard::FastqToSam->create(
            fastq => $humanfree_fwd_path,
            fastq2 => $humanfree_rev_path,
            output => $humanfree_bam_path,
            quality_format => 'Standard',
            platform => 'illumina', # nnutter: Not sure what these are/should be?
            sample_name => $hcs_data->sample_name,
            library_name => $hcs_data->library_name,
            use_version => '1.21',
        );
        unless($humanfree_fastq_to_sam->execute()) {
            die $self->error_message("Failed to convert FastQ to BAM ($humanfree_fwd_path, $humanfree_rev_path).");
        }
    }
    if (! $self->overwrite && -f $original_bam_path) {
        $self->debug_message("\t$hcs_data_id: Skipping original bam creation, file already exists. Use --overwrite to replace...");
    }
    else {
        unlink($original_bam_path);
        $self->debug_message("\t$hcs_data_id: Generating original bam file...");

        $self->debug_message("\tVerifying fwd/rev pairs are correct, will swap if not...");

        # If original files are reversed then they probably just have rev in [0] and fwd in [1] so switch "pointer".
        my $original_fwd_path = @{$fastq_files{$hcs_data_id}{original}}[0];
        my $original_rev_path = @{$fastq_files{$hcs_data_id}{original}}[1];
        my $original_fwd_file = Genome::Sys->open_file_for_reading($original_fwd_path);
        my $original_fwd_line = $original_fwd_file->getline;
        if ($original_fwd_line =~ /\/2$/) {
            $self->debug_message("\t\t" . (split("/", $original_fwd_path))[-1] . " looks like a reverse file. Swapping...");
            my $tmp = @{$fastq_files{$hcs_data_id}{original}}[0];
            @{$fastq_files{$hcs_data_id}{original}}[0] = @{$fastq_files{$hcs_data_id}{original}}[1];
            @{$fastq_files{$hcs_data_id}{original}}[1] = $tmp;
            $original_fwd_path = @{$fastq_files{$hcs_data_id}{original}}[0];
            $original_rev_path = @{$fastq_files{$hcs_data_id}{original}}[1];
        }

        $self->debug_message("\t" . (split("/", $original_bam_path))[-1]. "...");
        my $original_fastq_to_sam = Genome::Model::Tools::Picard::FastqToSam->create(
            fastq => $original_fwd_path,
            fastq2 => $original_rev_path,
            output => $original_bam_path,
            quality_format => 'Standard',
            platform => 'illumina', # nnutter: Not sure what these are/should be?
            sample_name => $hcs_data->sample_name,
            library_name => $hcs_data->library_name,
            use_version => '1.21',
        );
        unless($original_fastq_to_sam->execute()) {
            die $self->error_message("Failed to convert FastQ to BAM ($original_fwd_path, $original_rev_path).");
        }
    }

    ### Genome::Model::Tools::Picard::EstimateLibraryComplexity (human-free) ###
    $self->debug_message("Running Picard EstimateLibraryComplexity report on human-free bam...");
    my $humanfree_report_path = $self->report_dir . '/' . $hcs_data_id . '_humanfree_untrimmed_estimate_library_complexity_report.txt';
    if (! $self->overwrite && -f $humanfree_report_path) {
        $self->debug_message("\t$hcs_data_id: Skipping EstimateLibraryComplexity, files already exists. Use --overwrite to replace...");
    }
    else {
        unlink($humanfree_report_path);

        $self->debug_message("\t$hcs_data_id: Generating EstimateLibraryComplexity report on human-free bam...");

        $self->debug_message("\t\t" . (split("/", $humanfree_report_path))[-1] . "...");
        my $humanfree_picard_elc = Genome::Model::Tools::Picard::EstimateLibraryComplexity->create(
            input_file => [$humanfree_bam_path],
            output_file => $humanfree_report_path,
            use_version => '1.21',
        );
        unless($humanfree_picard_elc->execute()) {
            die $self->error_message("Failed to convert run Picard ELC on humanfree.");
        }

    }

    ### Genome::Model::Tools::Picard::EstimateLibraryComplexity (original) ###
    $self->debug_message("Running Picard EstimateLibraryComplexity report on \"original bam\"...");
    my $original_report_path = $self->report_dir . '/' . $hcs_data_id . '_original_untrimmed_estimate_library_complexity_report.txt';
    if (! $self->overwrite && -f $original_report_path) {
        $self->debug_message("\t$hcs_data_id: Skipping EstimateLibraryComplexity, files already exists. Use --overwrite to replace...");
    }
    else {
        unlink($original_report_path);
        $self->debug_message("\t$hcs_data_id: Generating EstimateLibraryComplexity report on \"original\" bam...");

        $self->debug_message("\t\t" . (split("/", $original_report_path))[-1] . "...");
        my $original_picard_elc = Genome::Model::Tools::Picard::EstimateLibraryComplexity->create(
            input_file => [$original_bam_path],
            output_file => $original_report_path,
            use_version => '1.21',
        );
        unless($original_picard_elc->execute()) {
            die $self->error_message("Failed to convert run Picard ELC on \"original\".");
        }    
    }

    $self->debug_message("Removing unneeded FastQ files...");
    for my $hcs_data_id (keys %fastq_files) {
        last unless(exists($fastq_files{$hcs_data_id}) && exists($fastq_files{$hcs_data_id}{original}) && exists($fastq_files{$hcs_data_id}{imported}));
        unlink(@{$fastq_files{$hcs_data_id}{original}}[0]) if (-f @{$fastq_files{$hcs_data_id}{original}}[0]);
        unlink(@{$fastq_files{$hcs_data_id}{original}}[1]) if (-f @{$fastq_files{$hcs_data_id}{original}}[1]);
        unlink(@{$fastq_files{$hcs_data_id}{imported}}[0]) if (-f @{$fastq_files{$hcs_data_id}{imported}}[0]);
        unlink(@{$fastq_files{$hcs_data_id}{imported}}[1]) if (-f @{$fastq_files{$hcs_data_id}{imported}}[1]);

        unlink($humanfree_fwd_path) if (-f $humanfree_fwd_path);
        unlink($humanfree_rev_path) if (-f $humanfree_rev_path);
    }

    return 1;
}

sub other_stats {
    my ($self, $mcs_metrics_hash, $other_stats_output, $hcs_data) = @_;
    my $hcs_data_id = $hcs_data->id;
    my %mcs_metrics = %$mcs_metrics_hash;

    return unless ($hcs_data->is_paired_end);

    ### OTHER STATS ###
    # Pulled from the $contamination_flagstat including:
    #   the percent mapped
    #   duplicate count
    #   unique, non-human bases
    #   percent duplication of raw data

    # Count bases in humanfree, untrimmed bams
    my %humanfree_base_count;
    $self->debug_message("Counting human-free, untrimmed bases per lane...");

    $self->expect64();
    my $humanfree_bam_path = $self->report_dir . '/data/' . $hcs_data_id . '_humanfree_untrimmed.bam';
    my $bam_fh = IO::File->new("samtools view $humanfree_bam_path |");
    while (<$bam_fh>) {
        my $read = $_;
        my $bases = (split("\t", $read))[9];
        $humanfree_base_count{$hcs_data_id} += length($bases);
    }


    $self->debug_message("Generating other stats...");

    my $orig_data = $self->original_data_from_imported_id($hcs_data_id);
    my $lane = $orig_data->flow_cell_id . "_" . $orig_data->lane;

    my $humanfree_report_path = $self->report_dir . '/' . $hcs_data_id . '_humanfree_untrimmed_estimate_library_complexity_report.txt';
    my $humanfree_report_fh = Genome::Sys->open_file_for_reading($humanfree_report_path);
    while (my $line = $humanfree_report_fh->getline) {
        if ($line =~ /^##\ METRICS/) {
            my $keys = $humanfree_report_fh->getline();
            my $values = $humanfree_report_fh->getline();
            my @keys = split("\t", lc($keys));
            my @values = split("\t", $values);
            my %metrics;
            @metrics{@keys} = @values;
            print $other_stats_output "$lane: Human-free Percent Duplication: " . $metrics{percent_duplication} . "\n";

            my $metric_name = "$lane\_humanfree_percent_duplication";
            $mcs_metrics{$metric_name}->delete() if($mcs_metrics{$metric_name});
            unless(Genome::Model::Metric->create(build_id => $self->build_id, name => $metric_name, value => $metrics{percent_duplication})) {
                die $self->error_message("Unable to create build metric (build_id=" . $self->build_id . ", $metric_name)");
            }

            $metric_name = "$lane\_unique_humanfree_bases";
            $mcs_metrics{$metric_name}->delete() if($mcs_metrics{$metric_name});
            my $unique_bases_count = $humanfree_base_count{$hcs_data_id} * (1 - $metrics{percent_duplication});
            unless(Genome::Model::Metric->create(build_id => $self->build_id, name => $metric_name, value => $unique_bases_count)) {
                die $self->error_message("Unable to create build metric (build_id=" . $self->build_id . ", $metric_name)");
            }
            print $other_stats_output "$lane: Unique, human-free bases: $unique_bases_count\n";
        }
    }

    my $original_report_path = $self->report_dir . '/' . $hcs_data_id . '_original_untrimmed_estimate_library_complexity_report.txt';
    my $original_report_fh = Genome::Sys->open_file_for_reading($original_report_path);
    while (my $line = $original_report_fh->getline) {
        if ($line =~ /^##\ METRICS/) {
            my $keys = $original_report_fh->getline();
            my $values = $original_report_fh->getline();
            my @keys = split("\t", lc($keys));
            my @values = split("\t", $values);
            my %metrics;
            @metrics{@keys} = @values;
            print $other_stats_output "$lane: Original Percent Duplication: " . $metrics{percent_duplication} . "\n";

            my $metric_name = "$lane\_original_percent_duplication";
            $mcs_metrics{$metric_name}->delete() if($mcs_metrics{$metric_name});
            unless(Genome::Model::Metric->create(build_id => $self->build_id, name => $metric_name, value => $metrics{percent_duplication})) {
                die $self->error_message("Unable to create build metric (build_id=" . $self->build_id . ", $metric_name)");
            }
        }
    }

    return 1;
}

sub original_data_from_imported_id {
    my ($self, $id) = @_;
    my $imported_data = Genome::InstrumentData::Imported->get($id);
    (my $alignment_id = $imported_data->original_data_path) =~ s/.*\/([0-9]*)\/.*/$1/;
    return Genome::InstrumentData::AlignmentResult->get($alignment_id)->instrument_data;
}

sub read_and_join_lines {
    my ($fh, $num) = @_;
    $num = 4 unless ($num);

    my @lines;
    for (my $count = 0; $count < 4; $count++) {
        my $line = $fh->getline;
        return undef unless $line;
        push @lines, $line;
    }
    return join('', @lines);
}

sub expect64 {
    my $self = shift;
    my $uname = `uname -a`;
    unless ($uname =~ /x86_64/) {
        die $self->error_message("Samtools requires a 64-bit operating system.");
    }
    return 1;
}

sub bam_stats_per_lane {
    my $self = shift;
    my $bam = shift;

    #$self->debug_message("Opening $bam with samtools...");
    $self->expect64();
    my $bam_fh = IO::File->new("samtools view $bam |");
    unless($bam_fh) {
        die $self->error_message("Failed to open $bam for reading.");
    }

    my %stats;
    my %flow_lane;
    my $total_reads_all_lanes;
    while (<$bam_fh>){
        my @fields = split("\t", $_);

        my $seq = $fields[9];
        (my $id = $fields[11]) =~ s/.*://;
        #print "$total_reads_all_lanes: $id: $seq...";

        # cache the flow_lane -> ID mapping
        unless($flow_lane{$id}) {
            my $data = Genome::InstrumentData->get($id);
            unless($data) {
                die $self->error_message("Unable to find data (imported nor original) by ID $id.");
            }
            $flow_lane{$id} = $data->flow_cell_id . '_' . $data->lane;
        }

        if ($seq eq '*'){
            $self->debug_message("invalid sequence for read in bam: $_");
            next;
        }
        $stats{$flow_lane{$id}}{total_bases} += length($seq);

        # percent mapped and duplication
        my $flag = $fields[1];
        if($flag & 0x0001) {
            $stats{$flow_lane{$id}}{paired_in_sequencing}++;                                    # the read is paired in sequencing, no matter whether it is mapped in a pair
            $stats{$flow_lane{$id}}{properly_paired}++              if(    $flag & 0x0002);     # the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
            $stats{$flow_lane{$id}}{read1}++                        if(    $flag & 0x0040);     # the read is the first read in a pair 1,2
            $stats{$flow_lane{$id}}{read2}++                        if(    $flag & 0x0080);     # the read is the second read in a pair 1,2
        }
        $stats{$flow_lane{$id}}{mapped}++                           unless($flag & 0x0004);     # the query sequence itself is unmapped
        $stats{$flow_lane{$id}}{qc_failure}++                       if(    $flag & 0x0200);     # the read fails platform/vendor quality checks
        $stats{$flow_lane{$id}}{duplicates}++                       if(    $flag & 0x0400);     # the read is either a PCR duplicate or an optical duplicate

        $total_reads_all_lanes++;
        $stats{$flow_lane{$id}}{total_reads}++;
        $self->debug_message("\t\tProcessed " . $total_reads_all_lanes/1e6 . "M reads so far...") unless ($total_reads_all_lanes % 1e6);
        #print " done\n";
    }

    for my $id (keys %flow_lane) {
        $stats{$flow_lane{$id}}{paired_in_sequencing} ||= 0;
        $stats{$flow_lane{$id}}{properly_paired} ||= 0;
        $stats{$flow_lane{$id}}{singletons} ||= 0;
        $stats{$flow_lane{$id}}{with_itself_and_mate_mapped} ||= 0;
        $stats{$flow_lane{$id}}{read1} ||= 0;
        $stats{$flow_lane{$id}}{read2} ||= 0;
        $stats{$flow_lane{$id}}{mapped} ||= 0;
        $stats{$flow_lane{$id}}{qc_failure} ||= 0;
        $stats{$flow_lane{$id}}{duplicates} ||= 0;
        $stats{$flow_lane{$id}}{average_read_length} = $stats{$flow_lane{$id}}{total_bases}/$stats{$flow_lane{$id}}{total_reads};
        $stats{$flow_lane{$id}}{percent_mapped} = $stats{$flow_lane{$id}}{mapped}/$stats{$flow_lane{$id}}{total_reads};
        $stats{$flow_lane{$id}}{percent_duplicates} = $stats{$flow_lane{$id}}{duplicates}/$stats{$flow_lane{$id}}{total_reads};
        $stats{$flow_lane{$id}}{percent_properly_paired} = $stats{$flow_lane{$id}}{properly_paired}/$stats{$flow_lane{$id}}{total_reads};
    }
    return %stats;
}

1;

