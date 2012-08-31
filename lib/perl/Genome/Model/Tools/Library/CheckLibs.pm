package Genome::Model::Tools::Library::CheckLibs;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Library::CheckLibs {
    is => 'Command',
    has => [
    group => { 
        type => 'String',
        is_optional => 1,
        doc => "model group to report outliers for",
    },
    builds => {
        type => 'String',
        is_optional => 1,
        doc => "builds to report outliers for",
    },
    above_cutoff => {
        type => 'float',
        is_optional => 1,
        default => 0.2,
    },
    below_cutoff => {
        type => 'float',
        is_optional => 1,
        default => 0.3,
    },
    output_file => {
        type => 'String',
        is_optional => 1,
        doc => "file to write output to",
    },

    ]
};


sub execute {
    my $self=shift;
    my @builds;
    if($self->builds) {
        @builds = map { Genome::Model::Build->get($_); }  split /\s+/, $self->builds;
    }
    elsif($self->group) {
        my $group = Genome::ModelGroup->get(name => $self->group);
        unless($group) {
            $self->error_message("Unable to find a model group named " . $self->group);
            return;
        }
        for my $model ($group->models) {
            my $build = $model->last_complete_build;
            unless($build) {
                $self->error_message("No complete build for model " . $model->id);
            }
            else {
                push @builds, $build;
            }
        }
    }
    else {
        $self->error_message("You must provide either build id(s) or a model group name to run this script");
        return;
    }

    #set up output file here
    my $output_fh;
    if($self->output_file) {
        $output_fh = IO::File->new($self->output_file,"w");
        unless($output_fh) {
            $self->error_message("Couldn't open " . $self->output_file . " for writing");
            return;
        }
    }
    else {
        $output_fh = IO::File->new_from_fd(1,"w");
        unless($output_fh) {
            $self->error_message("Couldn't open STDOUT for output");
            return;
        }
    }
        
    print $output_fh "Lanes indicating bad libraries\n";
    print $output_fh join("\t",qw(Common_Name Name Library Insert_Size SD_Below SD_Above)), "\n";
    foreach my $build (@builds) {
        my $model = $build->model;
        unless(defined($model)) {
            $self->error_message("Somehow this build does not have a model");
            return;
        }

        #calculate common name like AML11
        my $common_name = $model->subject->source_common_name;

        #this should tell us about whether it's tumor or normal
        my $type = $model->subject->common_name;

        printf STDERR "Grabbing information for model %s (build %s)\n", $model->name, $build->build_id;       
        #Grab all alignment events so we can filter out ones that are still running or are abandoned
        # get all align events for the current running build
        my @align_events = Genome::Model::Event->get(event_type => 
            {operator => 'like', value => '%align-reads%'},
            build_id => $build,
            model_id => $model->id,
        );
        printf STDERR "%d lanes in build\n", scalar(@align_events);
        
        my @inputs = map { $_->instrument_data_input } @align_events;

        my %bad_lanes;
        my %total_lanes;
        
        foreach my $input (@inputs) {
            my $instrument_data = $input->value;
            my $library = $instrument_data->library_name;
            unless(defined($library)) {
                $self->error_message("No library defined for " . $instrument_data->__display_name__);
                next;
            }

            my $ispe = ($instrument_data->is_paired_end && ! defined($input->filter_desc));
            next unless $ispe;

            $total_lanes{$library}++;
            unless(exists($bad_lanes{$library})) {
                $bad_lanes{$library} = 0;
            }

            my $lane_name = $instrument_data->short_name."_".$instrument_data->subset_name;
            my $median_insert_size = $instrument_data->median_insert_size;
            $median_insert_size ||= 0;
            my $sd_above_insert_size = $instrument_data->sd_above_insert_size;
            my $sd_below_insert_size = $instrument_data->sd_below_insert_size;
            if(!$median_insert_size || ($sd_below_insert_size/$median_insert_size) > 0.3 || ($sd_above_insert_size/$median_insert_size) > 0.2) {
                printf $output_fh "%s\t%s\t%s\t%0.2f\t%0.2f\t%0.2f\n","$common_name.$type",$lane_name,$library,$median_insert_size, $median_insert_size ? $sd_below_insert_size/$median_insert_size : 0, $median_insert_size ? $sd_above_insert_size/$median_insert_size : 0;
                $bad_lanes{$library}++;
            }
        }
        foreach my $lib (keys %total_lanes) {
                
            printf $output_fh "%0.2f%% lanes bad for %s %s library %s\n",$bad_lanes{$lib}/$total_lanes{$lib}*100,$common_name,$type,$lib;
        }
    }
    return 1;

}


1;

sub help_brief {
    "Checks to make sure that lanes in each library have insert size distributions within the normal range"
}

sub help_detail {
    <<'HELP';
This script uses the Genome Model API to grab out all alignment events for a model and checks the GERALD insert size metrics for each library. Outliers are flagged and the number of "bad" lanes per library is reported. This program ignores lanes which do not have a library name. Be aware that occasionally, the insert size metrics do not get loaded into the database and thus will have 0's. This is not necessarily an indication that the lane or flowcell is bad. Standard deviations are reported as fractions of the median insert size. Please note that this only queries the GERALD metrics and thus it may report as "bad" lanes which were aligned as fragment.
HELP
}
