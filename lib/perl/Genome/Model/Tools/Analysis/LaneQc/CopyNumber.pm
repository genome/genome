package Genome::Model::Tools::Analysis::LaneQc::CopyNumber;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Analysis::LaneQc::CopyNumber {
    is => 'Command',
    has => [
        build_id => { 
            type => 'String',
            is_optional => 0,
            doc => "build id of the build to do per-lane copy-number QC on",
        },
        bam2cn_window => { 
            type => 'Number',
            is_optional => 1,
            default => 50000,
            doc => "Window (in bp) for looking at read-depth window-based copy number.",
        },
        output_file_prefix => {
            type => 'String',
            is_optional => 0,
            doc => "Prefix for filename to write per-lane copy-number QC data to. The name of the lane plus \".cnqc\" will be added as a suffix for all outputs if you give a filename, or the name of the lane will be the beginning of the filename if you just give a path, and also \".png\" will be further added for a graphical output as well. Use full path!!",
        },
        lsf => {
            is => 'Boolean',
            default => 1,
            doc => "Send internal tasks to LSF, e.g. R's plot generation and BAM2CN.pl."
        },
    ],
};


sub execute {
    my $self=shift;
    my $build_id = $self->build_id;
    my $outfile_prefix = $self->output_file_prefix;
    my $window = $self->bam2cn_window;

    my $build = Genome::Model::Build->get($build_id);
    unless(defined($build)) {
        $self->error_message("Unable to find build $build_id");
        return;
    }
    my $model = $build->model;
    unless(defined($model)) {
        $self->error_message("Somehow this build does not have a model");
        return;
    }
    printf STDERR "Grabbing information for model %s (build %s)\n", $model->name, $build->build_id;       
    #Grab all alignment events so we can filter out ones that are still running or are abandoned
    # get all align events for the current running build
    my @align_events = Genome::Model::Event->get(event_type => 
        {operator => 'like', value => '%align-reads%'},
        build_id => $build,
        model_id => $model->id,
    );
    printf STDERR "%d lanes in build\n", scalar(@align_events);
    #now just get the Succeeded events to pass along for further processing
    # THIS MAY NOT INCLUDE ANY EVENTS
    my @events = Genome::Model::Event->get(event_type => 
        {operator => 'like', value => '%align-reads%'},
        build_id => $build,
        event_status => 'Succeeded',
        model_id => $model->id,
    );
    # if it does not include any succeeded events - die
    unless (@events) {
        $self->error_message("No alignments have Succeeded on the build");
        return;
    }
    $self->status_message(sprintf "Checking %d lanes", scalar(@events));

    for my $instrument_data ($build->instrument_data) {
        my $lane_name = $instrument_data->__display_name__;
        my @alignments = $build->alignment_results_for_instrument_data($instrument_data);
        unless(@alignments) {
            $self->error_message("No alignment objects for $lane_name");
            return;
        }
        for my $alignment (@alignments) {
            my $alignment_id = $alignment->id;
            my $reference = $alignment->reference_build->full_consensus_path('fa');
            my $instrument_data_id = $alignment->instrument_data_id;    
            my @bams = $alignment->alignment_bam_file_paths;
            unless(@bams) {
                $self->error_message("No alignment bam for $alignment_id");
                return;
            }
            else {
                my $alignment_count = @bams;
                my $path = $bams[0];
                $self->status_message("Found $alignment_count for $alignment_id with default path $path\n");
            }
            my $alignment_file = $bams[0];
            unless(-e $alignment_file) {
                $self->error_message("$alignment_file does not exist");
                return;
            }

            my $user = Genome::Sys->username;

            my $lane_outfile = $outfile_prefix . $alignment_id . ".cnqc";

            my $job1_name = $lane_outfile . "-cn-qc";
            my $job2_name = $job1_name . "-plot";
            my $dependency = "ended($job1_name)";

            my $module_dir = $self->module_dir;

            my $bam2cn_pl = "$module_dir/BAM2CN.pl";
            die $self->error_message("Missing Perl script ($bam2cn_pl)") unless (-s "$bam2cn_pl");
            my $cmd1 = $^X . " $bam2cn_pl -w $window $alignment_file > $lane_outfile";

            my $plot_wholegenome_cn_R = "$module_dir/plot_wholegenome_cn.R";
            die $self->error_message("Missing R script ($plot_wholegenome_cn_R)") unless (-s "$plot_wholegenome_cn_R");
            my $cmd2 = "R --no-save < $plot_wholegenome_cn_R $lane_outfile";

            my $rv;
            if ($self->lsf) {
                $rv = system("bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J $job1_name -R 'select[type==LINUX64]' \"$cmd1\"");
                if ($rv) {
                    $self->error_message("Failed to bsub $job1_name");
                    return;
                }
                $rv = system("bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J $job2_name -w \"$dependency\" \"$cmd2\"");
                if ($rv) {
                    $self->error_message("Failed to bsub $job2_name");
                    return;
                }
            } else {
                $rv = system($cmd1);
                if ($rv) {
                    $self->error_message("$cmd1 failed");
                    return;
                }
                $rv = system($cmd2);
                if ($rv) {
                    $self->error_message("$cmd2 failed");
                    return;
                }
            }

        }
    }

    return 1;

}

sub module_dir {
    my $path = __FILE__;
    my ($dir) = $path =~ /(.*)\//;
    return $dir;
}

sub _dedupe_objects {
    my @objects = @_;
    my %seen = ();
    my @deduped_objects = grep { ! $seen{$_->id} ++ } @objects;
    return @deduped_objects;
}


sub help_brief {
    "runs per-lane copy-number QC on every aligned lane in a build"
}


sub help_detail {
    <<'HELP';
    This script looks through a build, and for each aligned lane, submits a job to create a single-genome copy-number plot. (You do not need to bsub this command - it will bsub other jobs for you.) This plot can be compared to a snp-array copy-number plot from the same sample to determine if there are any lanes which appear to not coincide with the snp-array plot, and thus might not belong to the same sample.
HELP
}

1;
