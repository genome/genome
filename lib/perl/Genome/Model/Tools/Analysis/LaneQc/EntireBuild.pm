package Genome::Model::Tools::Analysis::LaneQc::EntireBuild;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Analysis::LaneQc::EntireBuild {
    is => 'Command',
    has => [
    build_id => { 
        type => 'String',
        is_optional => 0,
        doc => "build id of the build to do per-lane SNP QC on",
    },
    analysis_dir => {
        type => 'String',
        is_optional => 0,
        doc => "Directory to write per-lane SNP QC to",
    },
    genotype_file => {
        type => 'String',
        is_optional => 0,
        doc => "Genotype file to use as input to gmt analysis lane-qc compare-snps",
    },

    ]
};


sub execute {
    my $self=shift;
    my $build_id = $self->build_id;
    my $genotype_file = $self->genotype_file;

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
    $self->debug_message(sprintf "Checking %d lanes", scalar(@events));
    for my $instrument_data ($build->instrument_data) {
        my $lane_name = $instrument_data->__display_name__;
        my @alignments = $build->alignment_results_for_instrument_data($instrument_data);
        unless(@alignments) {
            $self->error_message("No alignment objects for $lane_name");
            return;
        }
        for my $alignment (@alignments) {
            my $reference = $alignment->reference_build->full_consensus_path('fa');
            my $instrument_data_id = $alignment->instrument_data_id;    
            my @bams = $alignment->alignment_bam_file_paths;
            unless(@bams) {
                $self->error_message("No alignment bam for $lane_name");
                return;
            }
            else {
                my $alignment_count = @bams;
                my $path = $bams[0];
                $self->status_message("Found $alignment_count for $lane_name with default path $path\n");
            }
            my $alignment_file = $bams[0];
            unless(-e $alignment_file) {
                $self->error_message("$alignment_file does not exist");
                return;
            }
            my $dir = $self->analysis_dir;
            my $user = Genome::Sys->username;
            if(-z "$dir/$lane_name.var" || !-e "$dir/$lane_name.var") {
                my $command .= <<"COMMANDS";
samtools pileup -vc -f $reference $alignment_file | perl -pe '\@F = split /\\t/; \\\$_=q{} unless(\\\$F[7] > 2);' > $dir/$lane_name.var
gmt analysis lane-qc compare-snps --genotype-file $genotype_file --variant-file $dir/$lane_name.var > $dir/$lane_name.var.compare_snps
COMMANDS
                print `bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -R 'select[type==LINUX64]' "$command"`;
                #print $command,"\n";
            }
        }

    }
    return 1;

}


1;

sub help_brief {
    "runs per-lane SNP QC on every aligned lane in a build"
}

sub help_detail {
    <<'HELP';
This is a relatively crummy script to aid in running per-lane SNP QC until it is integrated into the pipeline. For every lane aligned in a build, it will produce a .var file of SNP calls and a .compare_snps file containing the output of gmt analysis lane-qc compare-snps. Each file will be named by the convention FlowCellId_Lane
HELP
}
