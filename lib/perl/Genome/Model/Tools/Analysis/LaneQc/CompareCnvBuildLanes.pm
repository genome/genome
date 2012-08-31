package Genome::Model::Tools::Analysis::LaneQc::CompareCnvBuildLanes;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Analysis::LaneQc::CompareCnvBuildLanes {
    is => 'Command',
    has => [
    model_id => { type => 'String', is_optional => 0, doc => "tumor/normal model_id to get the last suceed build to gather per lane bam files",},
    outfile_prefix => { type => 'String', is_optional => 0, doc => "Directory to use for keep outputs",},
    bam2cn_window => { type => 'Number', is_optional => 1, default => 50000, doc => "Window (in bp) for looking at read-depth window-based copy number.",},
    ]
};

sub help_brief {
    "Generates QC CNV plot for every lane in a build for one wgs model"
}

sub help_detail {
    <<'HELP';
This script runs QC check on every lane in a build to generate per lane CNV plot and detect sample swaps. It may also be useful for analysis of quality metrics on a per lane basis.
HELP
}

###################
#
# Code was modified based on Dave's version to fit new pipeline.
#
###################
sub execute {
    my $self=shift;
    my $build_id ="";
    my $outfile_prefix = $self->outfile_prefix;
    my $user = getlogin || getpwuid($<); #get current user name	 
    my $wgs_model_id = $self->model_id;
    my $window = $self->bam2cn_window;
    # step2: To find alignment file of the build or return;
    # grap the last succeed build of the wgs model to find the alignment bam files

    my $wgs_model;
    my $build;

    $wgs_model = Genome::Model-> get(genome_model_id =>$wgs_model_id);
    unless(defined($wgs_model)){
        $self->error_message("Unable to file model $wgs_model_id");
        return;
    }        

    $build= $wgs_model->last_succeeded_build;  
    unless(defined($build)){
        $self->error_message("Unable to find build $build");
        return;
    }
    $build_id=$build->build_id;

    my $model=$wgs_model;

#        print "model:$model\n build:$build\n build_id:$build_id\n model_id:$wgs_model_id\n";

    #Grab all alignment events so we can filter out ones that are still running or are abandoned
    # get all align events for the current running build
=cut    	my @align_events = Genome::Model::Event->get(
#    		event_type => {operator => 'like', value => '%align-reads%'},
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
            $self->error_message(" No alignments have Succeeded on the build ");
            return;
        }
        printf STDERR "Using %d lanes to calculate metrics\n", scalar(@events);
=cut
    my @instrument_data = $build->instrument_data;

# find reference sequences
#	my $reference_file="/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa"; 
    my $reference_file=$model->reference_sequence_build->full_consensus_path('fa') ;
    print "reference: $reference_file User: $user\n";
    # print "Number of idas:$#idas\n";
    for my $idata (@instrument_data) {
        my @alignments = $build->alignment_results_for_instrument_data($idata);
        for my $alignment (@alignments) {
            my $instrument_data = $alignment->instrument_data;
            my $lane=$instrument_data->lane;
            my $flow_cell_id=$instrument_data->flow_cell_id;
            my $lane_name="$flow_cell_id"."_"."$lane";	
            my @bam = $alignment->alignment_bam_file_paths;
            my $alignment_file = $bam[0];
            #	print "$lane_name\t$alignment_file\n";
            if ($alignment_file ne ""){
                $self->error_message("$lane_name : $alignment_file");
                unless(-e $alignment_file) {
                    $self->error_message("$alignment_file does not exist");
                    return;
                }

                my $lane_outfile = $outfile_prefix . "." . $lane_name . ".cnqc";
                my $job1_name = $lane_outfile . "-cn-qc";
                my $job2_name = $job1_name . "-plot";
                my $dependency = "ended($job1_name)";

                my $module_dir = $self->module_dir;

                my $bam2cn_pl = "$module_dir/BAM2CN.pl";
                die $self->error_message("Missing Perl script ($bam2cn_pl)") unless (-s "$bam2cn_pl");
                my $cmd1 = "perl $bam2cn_pl -w $window $alignment_file > $lane_outfile";

                my $plot_wholegenome_cn_R = "$module_dir/plot_wholegenome_cn.R";
                die $self->error_message("Missing R script ($plot_wholegenome_cn_R)") unless (-s "$plot_wholegenome_cn_R");
                my $cmd2 = "R --no-save < $plot_wholegenome_cn_R $lane_outfile";

                print `bsub -N -u $user\@genome.wustl.edu -J $job1_name -R 'select[type==LINUX64]' "$cmd1"`;
                print `bsub -N -u $user\@genome.wustl.edu -J $job2_name -w "$dependency" "$cmd2"`;
            }else{
                $self->error_message("No alignment object for $lane_name");
                return;
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

1;


