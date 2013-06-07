package Genome::Model::Tools::Graph::MutationCountBarPlot;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Graph::MutationCountBarPlot {
    is => 'Command',
    has => [
        model_group_id => {
            type => 'Int',
            doc => "model group id",
            is_optional => 0,
        },
        output_file_path => {
            type => 'String',
            doc => "output plot file",
            is_optional => 0,
        },
        plot_title  => {
            type => 'String',
            doc => "Title of the plot",
            default => "Mutation Counts across Samples Grouped by Tiers",
            example_values => ['Mutation Counts across Samples Grouped by Tiers'],
            is_optional => 1
        },
        res_X => {
            type => 'Int',
            doc => 'plot resolution, X',
            default_value => 1024,
            is_optional => 1,
        },
        res_Y => {
            type => 'Int',
            doc => 'plot resolution, Y',
            default_value => 768,
            is_optional => 1,
        },
        screen_split =>{
            type => 'Int',
            doc => 'screen split',
            default_value => 5,
            is_optional => 1,
        },
    ],
};

sub help_brief {# keep this to just a few words
    "Plot mutation counts grouped by tiers"
}

sub help_synopsis {
<<'HELP';
    Search in a directory for model data builds and count the number of mutations. Plot a bar chart to show the mutation counts grouped by tiers.
HELP
}

sub help_detail {# this is what the user will see with the longer version of help.
<<'HELP';
    The tool takes a user supplied model group id to search for model data builds (last successful builds). It goes through the builds and save the mutation counts into a data table, which is later processed by a R script to plot the bar chart.
HELP
}


################################################################################################################
## MutationCountBarPlot.pm - Plot a group bar chart that shows the mutations of each sample categorized by tiers
## AUTHOR: Yang Li
## EMAIL: yli@genome.wustl.edu
## DATE: June 5 2013
################################################################################################################
sub execute {
    my $self = shift;
    my $model_group_id = $self->model_group_id;
    my $output_file_path = $self->output_file_path;
    my $plot_title = $self->plot_title;
    my $res_X = $self->res_X;
    my $res_Y = $self->res_Y;
    my $screen_split = $self->screen_split;
	my @directories;
	my @names;

	my @model_ids =  map {$_-> id} Genome::Model->get("model_groups.id" => $model_group_id);
	for (my $i=0; $i < @model_ids; $i++){
		my $model_id = $model_ids[$i];
		my $model = Genome::Model->get($model_id);
		my $name = $model -> subject -> name;
		my $build = $model -> last_succeeded_build;
		my $build_dir = $build->data_directory;
		
		
		push (@directories, $build_dir);
		push (@names, $name);
	}
	
	# Create temp files
	my $temp_data_path = Genome::Sys->create_temp_file_path;
	my $temp_R_path = Genome::Sys->create_temp_file_path;
	
	# Create temporary file to store the data
	open (my $temp_data, ">$temp_data_path") or die "Create temporary data file failed. $!.\n";
	
	# Create temporary R script to plot the graph
	open (my $temp_R, ">$temp_R_path") or die "Create temporary R script failed. $!.\n"; 
	
	# Print header for the $temp_data
	my $header = "\tTier1\tTier2\tTier3\tTier4\n"; 
	print $temp_data $header;
	
	# Loop through the input file and get the sample directories stored in a list
	for (my $i=0; $i<@model_ids; $i++){
		my $directory = $directories[$i];
		my $sb .= $names[$i]; # Construct one line of output, starting with the sample name
		
		for (my $i=1;$i<5;$i++){ # Iterate through all four tiers
			my $mutationFile = $directory."/effects/snvs.hq.tier$i.v1.annotated.top";
			if (!(-e $mutationFile)) { # Examine if file exists
				$mutationFile = 0; # If file doesn't exist, use a 0 flag
				$sb .= "\t0"; # Append 0 to the string
			} else {
				my $wc = `wc -l $mutationFile`;
				my @parts = split(' ', $wc);
				my $count = $parts[0];
				$sb .= "\t$count"; # Append a tier count to the string 
			}
		}
		$sb .= "\n";
		print $temp_data $sb; # Write to $temp_data
	}
	
	
	# Split the window into rows 
	my $splits = @model_ids / $screen_split;
	
	# Multiply the vertical resolution for multiple rows
	$res_Y *= $splits;

	# Prepare R script
	print $temp_R "mutations <- read.table('$temp_data_path');\n";
	print $temp_R "png(filename = '$output_file_path', width = $res_X, height = $res_Y, bg = 'white');\n";
	print $temp_R "par(mfrow=c($splits,1));\n";
	for (my $i = 1; $i<=@model_ids; $i=$i+$screen_split) {
		my $upperBound = $i+$screen_split-1;
		print $temp_R "mutations_matrix = t(as.matrix(mutations));\n";
		print $temp_R "barplot(mutations_matrix[,$i:$upperBound], col=rainbow(4), beside=T, main='$plot_title', ylab='Mutation Counts');\n";
	}
	print $temp_R "legend(140, 18000, c('Tier1', 'Tier2', 'Tier3', 'Tier4'), fill=rainbow(4));\n";
	print $temp_R "dev.off();\n";
	
	# Execute R script
	system ("R CMD BATCH $temp_R_path");

	# Close files
	close $temp_data;
	close $temp_R;

    return 1;
}

1;
