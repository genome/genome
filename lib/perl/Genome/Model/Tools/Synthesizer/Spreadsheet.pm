package Genome::Model::Tools::Synthesizer::Spreadsheet;

use strict;
use warnings;
use Statistics::Descriptive;
use Data::Dumper;
use Genome;

my $DEFAULT_CLUSTERS = '5000';

class Genome::Model::Tools::Synthesizer::Spreadsheet {
	is        => ['Genome::Model::Tools::Synthesizer::Base'],
	has_input => [
		input_stats_file => {
			is  => 'Text',
			doc => 'Input alignment statistics file generated from stats-generator',
		},
		input_intersect_file => {
			is  => 'Text',
			doc => 'Input intersect TSV file generated from annotate-cluster',
		},
		output_spreadsheet => {
			is => 'Text',
			is_output=> 1,
			doc =>'Output speadsheet containing statistics as well as annotation for each cluster in the stats file',
		},
		input_cluster_number => {
            is => 'Text',
            is_optional => 1,
            doc => 'Number of TOP Clusters (sorted by depth) to calculate statistcs',
            default_value => $DEFAULT_CLUSTERS,
	   },

	],
};

sub help_brief {
"Run the Synthesizer Spreadsheet module to consolidate metrics for alignment,coverage as well as annotation in a TSV file";
}

sub help_detail {
"Run the Synthesizer Spreadsheet module to consolidate alignment,coverage as well as annotation in a TSV file. Input files for this module are generated from Annotate-Cluster module as well Stats-Generator module";
}

sub execute {
    my $self = shift;
    my $stats 		= $self->input_stats_file;
    my $annotation 	= $self->input_intersect_file;
    my $output 		= $self->output_spreadsheet;
    my $clusters 	= $self->input_cluster_number;
    
    
    my ($sorted_temp_fh, $sorted_temp_name) = Genome::Sys->create_temp_file();
    
    my $cmd = 'head -'.$clusters.' '.$annotation .'> '.$sorted_temp_name;
    
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_intersect_file],
        output_files => [$sorted_temp_name],
        skip_if_output_is_present => 0,
    );
    
    my $output_fh = Genome::Sys->open_file_for_writing($output);

    	my %id_hash;
    	my %size_hash;
    	print $output_fh join("\t","Cluster",
				"Chr",
				"Start",
				"Stop",
				"Avg Depth",
				"Zenith Depth",
				"Length of Raw Cluster",
				"Positive Strand",
				"Negative Strand",
				"Log Normalization - head bin",
				"Log Normalization -per bin",
				"% Mismatches",
				"ZeroMM",
				"1MM",
				"2MM",
				"3MM",
				"4MM",
				"% 1st Pos MM ",
				"Avg MapQ",
				"Std Dev Map Q",
				"%Zero MapQ",
				"Avg BaseQ",
				"Major Subcluster Loci"
				);
		
	my $annotation_fh = Genome::Sys->open_file_for_reading( $sorted_temp_name);		
		
	while (my $anno_row = $annotation_fh->getline) 
		{
			chomp $anno_row;
			my @array_ids = split(/\t/,$anno_row);
			my $size = scalar(@array_ids); 
    		if ($anno_row =~ /^Cluster/)
    		{
    			print $output_fh "\t". join("\t",@array_ids[4 .. $size -1])."\n";
    			
    		}
    		else
    		{
    			my $anno_cluster = $array_ids[0];
    		    my $anno_chr	 = $array_ids[1];
    		   	my $anno_start 	 = $array_ids[2];
    		    my $anno_stop 	 = $array_ids[3];
    		    
				my $stats_fh = Genome::Sys->open_file_for_reading( $stats);
	    	    while (my $row = $stats_fh->getline) 
	    	    {	
    				chomp $row;		    
    				my ($name,$chr,$start,$stop) = split ("\t", $row);	
    		   		
    		   		if ($chr eq $anno_chr && $start eq $anno_start && $stop eq $anno_stop )
    		    	{
  						print $output_fh $row."\t". join("\t",@array_ids[4 .. $size -1])."\n";
    		    	}
    			}
    		}
    		
  		}
    	    	
return 1;
    
}
1;	

