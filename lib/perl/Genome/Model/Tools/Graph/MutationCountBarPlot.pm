package Genome::Model::Tools::Graph::MutationCountBarPlot;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Graph::MutationCountBarPlot {
    is => 'Command',
    has => [
        input_file_path => {
            type => 'String',
            doc => "input file path",
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
    The tool takes a user supplied list of input directories to search for model data builds. (Each line contains and only contains a build directory.) It goes through the builds and save the mutation counts into a data table, which is later processed by a R script to plot the bar chart.
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
    my $input_file_path = $self->input_file_path;
    my $output_file_path = $self->output_file_path;
    my $plot_title = $self->plot_title;
    my $res_X = $self->res_X;
    my $res_Y = $self->res_Y;
	
	# Create temp files
	my $temp_data_path = Genome::Sys->create_temp_file_path;
	my $temp_R_path = Genome::Sys->create_temp_file_path;
	
	# Read the input file
	open (my $input_file, "<$input_file_path") or die "Read input directory failed. $!.\n";
	
	# Create temporary file to store the data
	open (my $temp_data, ">$temp_data_path") or die "Create temporary data file failed. $!.\n";
	
	# Create temporary R script to plot the graph
	open (my $temp_R, ">$temp_R_path") or die "Create temporary R script failed. $!.\n"; 
	
	# Print header for the $temp_data
	my $header = "\tTier1\tTier2\tTier3\tTier4\n"; 
	print $temp_data $header;
	
	# Loop through the input file and get the sample directories stored in a list
	while (my $line = <$input_file>){
		chomp $line;
		my $position = rindex($line, "/") + 1;
		my $ending = substr($line, $position);
		my $sb .= $ending; # Construct one line of output, starting with the sample name
		
		for (my $i=1;$i<5;$i++){ # Iterate through all four tiers
			my $tierLine = $line."/effects/snvs.hq.tier$i.v1.annotated.top";
			if (!(-e $tierLine)) { # Examine if file exists
				$tierLine = 0; # If file doesn't exist, use a 0 flag
			}
			my $wc = `wc -l $tierLine`;
			my @parts = split(' ', $wc);
			my $count = $parts[0];
			$sb .= "\t$count"; # Append a tier count to the string 
		}
		$sb .= "\n";
		print $temp_data $sb; # Write to $temp_data
	}
	
	# Prepare R script
	print $temp_R "mutations <- read.table('$temp_data_path');\n";
	print $temp_R "png(filename = '$output_file_path', width = $res_X, height = $res_Y, bg = 'white');\n";
	print $temp_R "barplot(t(as.matrix(mutations)), col=rainbow(4), beside=T, main='$plot_title', ylab='Mutation Counts');\n";
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
