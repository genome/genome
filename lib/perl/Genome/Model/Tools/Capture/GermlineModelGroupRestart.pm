
package Genome::Model::Tools::Capture::GermlineModelGroupRestart;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ModelGroup - Build Genome Models for Germline Capture Datasets
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	2/09/2011 by W.S.
#	MODIFIED:	2/09/2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

my %already_reviewed = ();
my %wildtype_sites = my %germline_sites = ();

class Genome::Model::Tools::Capture::GermlineModelGroupRestart {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_ids		=> { is => 'Text', doc => "IDs of model group(comma separated)" , is_optional => 0},
		output_file		=> { is => 'Text', doc => "Summarized Output" , is_optional => 0},
		rebuild		=> { is => 'Text', doc => "1 for rebuild, 0 for not and just make summary file" , is_optional => 0, default => 1},
                check_coverage	=> { is => 'Text', doc => "1 for check and dont act on models with low coverage, 0 for not and just check all failed models" , is_optional => 1, default => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Operate on capture somatic model groups"                 
}

sub help_synopsis {
    return <<EOS
Operate on capture somatic model groups
EXAMPLE:	gmt capture somatic-model-group --group-id 3328
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;
    my $group_ids = $self->group_ids;
    my $output_file = $self->output_file;
    my $rebuild_set = $self->rebuild;
    my $check_coverage = $self->check_coverage;
	# Open Output
	unless (open(OUTPUT,">$output_file")) {
	    die "Could not open output file '$output_file' for writing";
	  }

    my (@groups) = split /,/, $group_ids;
	foreach my $group_id (@groups) {
            print "Analyzing Group $group_id\n";
#                print OUTPUT "Analyzing Group $group_id\n";

        ## Get the models in each model group ##

        my $model_group = Genome::ModelGroup->get($group_id);
        my @model_bridges = $model_group->model_bridges;

        foreach my $model_bridge (@model_bridges) {
            my $model = Genome::Model->get($model_bridge->model_id);
	        my $model_id = $model->genome_model_id;
	        my $subject_name = $model->subject_name;
	        $subject_name = "Model" . $model_id if(!$subject_name);

	        if ($subject_name =~ m/Pooled/i) {
		        next;
	        }

	        my $build;
	        my $build_id;
	        my $last_build_dir;
	        my $bam_file = 'n/a';
                    my $build_status;

	        if($model->last_succeeded_build_directory) {
		        $build = $model->last_succeeded_build;
		        $build_id = $build->id;
                $build_status = $build->status;
		        $last_build_dir = $model->last_succeeded_build_directory;
		        $bam_file = $build->whole_rmdup_bam_file;
            }
            else {
		        $build = $model->latest_build;
                if (!defined $build && $rebuild_set) {
                    print OUTPUT join("\t", $model_id, "No_Build_Found", $subject_name, "No_Build_Found", "No_Build_Found") . "\n";
                    my $restart_command = "bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} \'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ -S genome model build start $model_id\'";
                    system($restart_command);
                    my $shortqueue_pending=`bjobs -q $ENV{GENOME_LSF_QUEUE_SHORT} | wc -l`; chomp $shortqueue_pending;
                    if ($shortqueue_pending > 20) {
                        sleep(30);
                    }
                    next;
                }
		        $build_id = $build->id;
                $build_status = $build->status;
            }

            unless ($build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                #unstartable builds restarted without other checks if they have instrument data
                if ($build_status =~ m/Unstartable/i) {
                    my @instrument_data = $build->instrument_data;
                    my $size = @instrument_data;
                    if ($size) {
                        my $new_status = "$build_status-restartable";
                        print OUTPUT join("\t", $model_id, $build_id, $subject_name, $new_status) . "\n";
                        my $restart_command = "bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} \'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ -S genome model build start $model_id\'";
                        unless ($build_status =~ m/abandon/i || $build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                            if ($rebuild_set) {
                                print "Build Id To Abandon: $build_id\n";
                                my @builds = ($build);
                                Genome::Model::Build::Command::Abandon->execute(builds => \@builds);
                            }
                        }
                        unless ($build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                            if ($rebuild_set) {
                                system($restart_command);
                            }
                        }
                    }
                    else {
                        my $new_status = "$build_status-no_inst_data";
                        print OUTPUT join("\t", $model_id, $build_id, $subject_name, $new_status) . "\n";
                    }
                }
                elsif ($check_coverage) {
                    my $read_length = 100; #should calc this but this is good enough for most
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
            #			print Data::Dumper::Dumper($hash_ref);exit;
                        my $unique_target_aligned_bp = 'unique_target_aligned_bp';
                        my $unique_on_target_bp = $hash_ref->{$wingspan}->{$unique_target_aligned_bp};
                    }
                    my $wingspan_zero_stats_file = $build->stats_file($wingspan);
                    my $specific_stats_file = "NA";
                    if ($wingspan_zero_stats_file) {
                        $specific_stats_file = $wingspan_zero_stats_file;
                        my $hash_ref = $build->coverage_stats_summary_hash_ref;
                        my @minimum_depth = qw(1 5 10 15 20);
            #			print Data::Dumper::Dumper($hash_ref);
                        my @deletehash;

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
                        }
                        my $summary_stats_ref = join("\t",@depth_array1);
                        my $summary_stats_ref2 = join("\t",@depth_array2);
                        my $summary_stats_ref3 = join("\t",@depth_array3);
            #			print join("\t", $model_id, $build_id, $subject_name, $build_status, $specific_stats_file, $general_stats_file, $summary_stats_ref, $summary_stats_ref2) . "\n";
                        print OUTPUT join("\t", $model_id, $build_id, $subject_name, $build_status, $bam_file, $summary_stats_ref, $summary_stats_ref2, $summary_stats_ref3) . "\n";
                        if ($depth_array3[0] > 100000) {
                            my $restart_command = "bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} \'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ -S genome model build start $model_id\'";
                            unless ($build_status =~ m/abandon/i || $build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                                    if ($rebuild_set) {
                                        print "Build Id To Abandon: $build_id\n";
                                        my @builds = ($build);
                                        Genome::Model::Build::Command::Abandon->execute(builds => \@builds);
                                    }
                            }
                            unless ($build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                                    if ($rebuild_set) {
                                            system($restart_command);
                                    }
                            }
                        }
                    }
                }
                else {
                print OUTPUT join("\t", $model_id, $build_id, $subject_name, $build_status, $bam_file) . "\n";
                    my $restart_command = "bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} \'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ -S genome model build start $model_id\'";
                    unless ($build_status =~ m/abandon/i || $build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                        if ($rebuild_set) {
                            print "Build Id To Abandon: $build_id\n";
                            my @builds = ($build);
                            Genome::Model::Build::Command::Abandon->execute(builds => \@builds);
#                            system($abandon_command);
                        }
                    }
                    unless ($build_status =~ m/run/i || $build_status =~ m/schedule/i || $build_status =~ m/succeed/i) {
                        if ($rebuild_set) {
                            system($restart_command);
                        }
                    }
                }
            }
            else { #print succeeded, running, or scheduled
                print OUTPUT join("\t", $model_id, $build_id, $subject_name, $build_status, $bam_file) . "\n";
            }
            UR::Context->commit() or die 'commit failed';
            UR::Context->clear_cache(dont_unload => ['Genome::ModelGroup', 'Genome::ModelGroupBridge']);
	        }	
        }
	return 1;
}






sub byChrPos
{
	my ($chr_a, $pos_a) = split(/\t/, $a);
	my ($chr_b, $pos_b) = split(/\t/, $b);
	
	$chr_a cmp $chr_b
	or
	$pos_a <=> $pos_b;
}


1;

