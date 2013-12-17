package Genome::Model::Tools::Analysis::LaneQc::CompareSnpsBuildLanes;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Analysis::LaneQc::CompareSnpsBuildLanes {
    is => 'Command',
    has => [
    model_id => { type => 'String', is_optional => 0, doc => "tumor/normal model_id to get the last suceed build to gather per lane bam files",},
    analysis_dir => { type => 'String', is_optional => 0, doc => "Directory to use for keep outputs",},
    genotype_file => { type => 'String', is_optional => 1, doc => "Genotype file to use as input to gmt analysis lane-qc compare-snps",},
    sample_name => { type => 'String', is_optional => 0, doc => "Sample name to get imported genotype file, for example H_LC-SJTALL001-G-TB-01-1378 ",},
    ]
};

sub help_brief {
    "Generates QC gold-snp-concordance data on every lane in a build for one wgs model"
}

sub help_detail {
    <<'HELP';
This script runs QC check on every lane in a build to generate gold-snp-concordance and detect sample swaps. It may also be useful for analysis of quality metrics on a per lane basis.
HELP
}

###################
#
# Code was modified based on Dave's version to fit new pipeline.
#
###################
sub execute {
    my $self=shift;
    my $build_id = "";
    my $dir = $self->analysis_dir;
    my $user = getlogin || getpwuid($<); #get current user name
    my $sample_name ="";
    my $genotype_file ="";	 
    my $wgs_model_id = $self->model_id;
    $sample_name = $self->sample_name;
    # step1 : to find genotype file or return;
    $genotype_file = $self->genotype_file;
    if ($sample_name ne "" && $genotype_file eq ""){
        # get owner_id of the microarray_genotype file
        system("genome instrument-data list imported --filter sample_name=$sample_name --noheader | cut -d ' ' -f1 > /tmp/$sample_name");
        open (FH, "/tmp/$sample_name");
        my @owner_id=<FH>;
        my $owner_id=$owner_id[0];
        $owner_id=~s/\s+//;
        print $owner_id;
        close FH;
        # get the path of genotype file
        system("genome disk allocation list  --noheader --filter owner_id=$owner_id | cut -d ' ' -f1 > /tmp/$sample_name.path");
        open (FH2, "/tmp/$sample_name.path"); 
        my @path=<FH2>;
        my $path=$path[0];
        $path=~s/\s+//;
        my $find_genotype_file=$path."/".$sample_name.".genotype";
        print "\n$find_genotype_file\n";
        close FH2;
        system("rm /tmp/$sample_name");
        system("rm /tmp/$sample_name.path");
        # check if genotype_file exists
        if (!(-e $find_genotype_file) && !( -e $genotype_file)){
            $self->error_message("Unable to find genotype file $find_genotype_file and $genotype_file\n please check and supply path to --genotype-file");
            return;
        }		
        if (-e $find_genotype_file ){
            $genotype_file=	$find_genotype_file;
        }
    }

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
    my @align_events = Genome::Model::Event->get(
        event_type => {operator => 'like', value => '%align-reads%'},
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
    my @instrument_data = $build->instrument_data;

# find reference sequences
#	my $reference_file="/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa"; 
    my $reference_file=$model->reference_sequence_build->full_consensus_path('fa') ;
    print "reference: $reference_file User: $user\n";
    for my $data (@instrument_data) {
        my @alignments = $build->alignment_results_for_instrument_data($data);
        for my $alignment (@alignments) {
            my $instrument_data = $alignment->instrument_data;
            my $lane=$instrument_data->lane;
            my $flow_cell_id=$instrument_data->flow_cell_id;
            my $lane_name="$flow_cell_id"."_"."$lane";	
            my @bam = $alignment->alignment_bam_file_paths;
            my $alignment_file = $bam[0];

            if ($alignment_file ne ""){
                $self->error_message("$lane_name : $alignment_file");
                unless(-e $alignment_file) {
                    $self->error_message("$alignment_file does not exist");
                    return;
                }
                my $command .= <<"COMMANDS";
samtools pileup -vc -f $reference_file $alignment_file | perl -pe '\@F = split /\\t/; \\\$_=q{} unless(\\\$F[7] > 2);' > $dir/$lane_name.var
gmt analysis lane-qc compare-snps --genotype-file $genotype_file --variant-file $dir/$lane_name.var > $dir/$lane_name.var.compare_snps
COMMANDS
                print `bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -R 'select[type==LINUX64]' "$command"`;
            }else{
                $self->error_message("No alignment object for $lane_name");
                return;
            }
        }	
    }
    return 1;
}

1;


