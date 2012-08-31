package Genome::Model::Tools::Capture::GermlineModelGroup;
use strict;
use warnings;
use FileHandle;
use Genome;

my %stats = ();
my %already_reviewed = ();
my %wildtype_sites = my %germline_sites = ();

class Genome::Model::Tools::Capture::GermlineModelGroup {
	is => 'Command',                       
	has => [
		group_id		=> { is => 'Text', doc => "ID of model group" , is_optional => 0},
		output_build_dirs	=> { is => 'Text', doc => "If specified, outputs last succeeded build directory for each sample to this file" , is_optional => 1},
		output_coverage_stats	=> { is => 'Text', doc => "Specify a directory to output coverage stats" , is_optional => 1},
		output_flowcell_information	=> { is => 'Text', doc => "Specify a directory to output flowcell information" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Operate on capture somatic model groups"                 
}

sub help_synopsis {
    return help_brief()."\nEXAMPLE: gmt capture somatic-model-group --group-id 3328";
}

sub execute {
    my $self = shift;

    if($self->output_build_dirs) {
        open(BUILDDIRS, ">" . $self->output_build_dirs) or die "Can't open outfile: $!\n";
    }

    if($self->output_coverage_stats) {
        open(OUTFILE, ">" . $self->output_coverage_stats) or die "Can't open output file: $!\n";
        print OUTFILE join("\t",
            "Model_id",
            "Build_id",
            "Subject_name",
            "Build_Dir",
            "Coverage_File",
            "Coverage_Wingspan0_Depth1x",
            "Coverage_Wingspan0_Depth5x",
            "Coverage_Wingspan0_Depth10x",
            "Coverage_Wingspan0_Depth15x",
            "Coverage_Wingspan0_Depth20x",
            "Percent_Target_Space_Covered_1x",
            "Percent_Target_Space_Covered_5x",
            "Percent_Target_Space_Covered_10x",
            "Percent_Target_Space_Covered_15x",
            "Percent_Target_Space_Covered_20x",
            "Mapped_Reads",
            "Percent_Target_Space_Covered_20x_per1Gb",
            "Percent_Duplicates",
            "Mapping_Rate",
            "\n",
        );
    }

    if($self->output_flowcell_information) {
        open(FLOWCELL, ">" . $self->output_flowcell_information) or die "Can't open output file: $!\n";
    }


    for my $model_bridge (Genome::ModelGroup->get($self->group_id)->model_bridges) {
        my $model = $model_bridge->model;
        my $subject_name = $model->subject_name;
        $subject_name = "Model" . $model->id if(!$subject_name);

        if($model->last_succeeded_build_directory) {
            my $build = $model->last_succeeded_build;
            my $last_build_dir = $model->last_succeeded_build_directory;
            if($self->output_build_dirs) {
                my $bam_file = $build->whole_rmdup_bam_file;
                print BUILDDIRS join("\t", $model->id, $subject_name, $build->id, "Succeeded", $last_build_dir, $bam_file) . "\n";
                unless($self->output_coverage_stats) {
                    next;
                }
            }

            if($self->output_flowcell_information) {
#				Genome::InstrumentData->get($RG_id), and then calling flowcell id on that.
#				print $id->index_sequence
                my @instrument_data = $build->instrument_data;
                foreach my $instrument_data (@instrument_data) {
                    my $barcode = $instrument_data->index_sequence;
                    my $flowcell = $instrument_data->flow_cell_id;
                    my $lane = $instrument_data->lane;
                    my $read_length = $instrument_data->read_length;
                    print FLOWCELL join("\t", $subject_name, $model->id, $build->id, $flowcell, $lane, $barcode, $read_length) . "\n";
                }
            }

            if($self->output_coverage_stats) {
                my (@duplicate_pct, $mark_dup_hash_ref);
                my $dup_metric = 'PERCENT_DUPLICATION';
                my @data = $model->instrument_data;
                my %libname;
                foreach my $possible (@data) {
                    my $libname = $possible->library_name;
                    if ($libname =~ m/Pooled/i) {
                        next;
                    }
                    $libname{$libname}++;
                }
                if ($build->rmdup_metrics_file) {
                    foreach my $possible_libname (sort keys %libname) {
                        $mark_dup_hash_ref = $build->mark_duplicates_library_metrics_hash_ref;
                        push(@duplicate_pct, $mark_dup_hash_ref->{$possible_libname}->{$dup_metric});
                    }
                }
                else {
                    print $build->rmdup_metrics_file . " does not exist\n";
                    $duplicate_pct[0] = "N/A";
                }
                my $mark_dup = join(',', @duplicate_pct);

                my ($fs_stats, $mapping_pct, $readcount);
                my $fs_pct = 'reads_mapped_percentage';
                my $fs_readmap = 'reads_mapped';
                my $flagstat_file = $build->whole_rmdup_bam_flagstat_file;
                if (-e $flagstat_file) {
                    $fs_stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
                    $mapping_pct = $fs_stats->{$fs_pct};
                    $readcount = $fs_stats->{$fs_readmap};
                }
                my $read_length = 100; #should calc this but this is good enough for most
                my $gb_sequenced = $readcount * $read_length;

                my $wingspan = 0;
                my $average_depth = 'mean_depth';
                my $target_space_covered = 'pc_target_space_covered';
                my $covered_base_pair = 'covered_base_pair';
                #raw metrics per ROI:
                #Column 13 is the min_depth_filter column.  If you want 20x stats only then grep or parse the file for 20 in column 13.   For a definition of the output format try 'gmt ref-cov standard --help'.
                # summary of the alignment coverage stats:
                my $wingspan_zero_alignment_summary_file = $build->alignment_summary_file($wingspan);
                my $targeted_gbseq_normalized = "NA";
                my $general_stats_file = "NA";
                if (-e $wingspan_zero_alignment_summary_file) {
                    $general_stats_file = $wingspan_zero_alignment_summary_file;
                    my $hash_ref = $build->alignment_summary_hash_ref;
                    my $unique_target_aligned_bp = 'unique_target_aligned_bp';
                    my $unique_on_target_bp = $hash_ref->{$wingspan}->{$unique_target_aligned_bp};
                    if ($unique_on_target_bp > 0 && $gb_sequenced > 0) {
                        $targeted_gbseq_normalized = $unique_on_target_bp / $gb_sequenced;
                    }
                }
                my $wingspan_zero_stats_file = $build->stats_file($wingspan);
                my $specific_stats_file = "NA";
                if ($wingspan_zero_stats_file) {
                    $specific_stats_file = $wingspan_zero_stats_file;
#					# summary of the coverage stats:
#					my $wingspan_zero_summary_stats_file = $build->coverage_stats_summary_file($wingspan);
                    # summary metrics for wingspan 0 and 20x minimum depth:
                    my $hash_ref = $build->coverage_stats_summary_hash_ref;
                    my @minimum_depth = qw(1 5 10 15 20);

                    my @depth_array1;
                    my @depth_array2;
                    my @depth_array3;
                    my $target_space_20x_gbseq_normalized = 0;
                    foreach my $minimum_depth (@minimum_depth) {
                        my $summary_stats_ref = $hash_ref->{$wingspan}->{$minimum_depth}->{$average_depth};
                        my $summary_stats_ref2 = $hash_ref->{$wingspan}->{$minimum_depth}->{$target_space_covered};
                        my $summary_stats_ref3 = $hash_ref->{$wingspan}->{$minimum_depth}->{$covered_base_pair};
                        push(@depth_array1,$summary_stats_ref);
                        push(@depth_array2,$summary_stats_ref2);
                        push(@depth_array3,$summary_stats_ref3);
                        if ($minimum_depth == 20 && $summary_stats_ref2 > 0 && $gb_sequenced > 0) {
                            $target_space_20x_gbseq_normalized = ($summary_stats_ref3 / $gb_sequenced) * 1000000000;
                        }
                    }
                    my $summary_stats_ref = join("\t",@depth_array1);
                    my $summary_stats_ref2 = join("\t",@depth_array2);
                    my $summary_stats_ref3 = join("\t",@depth_array3);
                    print OUTFILE join("\t",
                        $model->id,
                        $build->id,
                        $subject_name,
                        $last_build_dir,
                        $specific_stats_file,
                        $general_stats_file,
                        $summary_stats_ref,
                        $summary_stats_ref2,
                        $gb_sequenced,
                        $targeted_gbseq_normalized,
                        $mark_dup,
                        $mapping_pct,
                        "\n",);
                } else {
                    print OUTFILE join("\t",
                        $model->id,
                        $build->id,
                        $subject_name,
                        $last_build_dir,
                        $specific_stats_file,
                        $general_stats_file,
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        $gb_sequenced,
                        $targeted_gbseq_normalized,
                        $mark_dup,
                        $mapping_pct,
                        "\n",);
                }
            }

#FUTURE: Pull SNP array concordance here
        }

        UR::Context->commit() or die 'commit failed';
        UR::Context->clear_cache(dont_unload => ['Genome::ModelGroup', 'Genome::ModelGroupBridge']);
    }	
    close(OUTFILE);

    return 1;
}

sub byChrPos {
    my ($chr_a, $pos_a) = split(/\t/, $a);
    my ($chr_b, $pos_b) = split(/\t/, $b);

    $chr_a cmp $chr_b
        or
    $pos_a <=> $pos_b;
}

1;
