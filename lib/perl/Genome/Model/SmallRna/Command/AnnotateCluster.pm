package Genome::Model::SmallRna::Command::AnnotateCluster;

## add curated cluster
## separate annotation in diff columns

use strict;
use warnings;
use Data::Dumper;
use Genome;

class Genome::Model::SmallRna::Command::AnnotateCluster {
    is => ['Genome::Model::SmallRna::Command::Base'],
    has_input => [
        cluster_bed_file => {
        	is => 'Text',
            doc => 'Input top N clusters BED file from Stats-Generator',
        },
        annotation_bed_file => {
            is => 'String',
            doc => 'Input BED file containing annotation,For more than one, supply a comma delimited list',
        },
        annotation_name => {
            is => 'String',
            doc => 'Comma delimited list of the Annotation Tracks. Should be in the same order as the list of annotation bed files.',
        },
        output_tsv_file => {
            is => 'Text',
            is_output=> 1,
            doc => 'Raw Output file from Intersectbed',
        },
    ],
};


sub execute {
    my $self = shift;
    my $cluster_bed = $self->cluster_bed_file;
    my @track_names = split(',',$self->annotation_name);	
	my @annotation  = split(',',$self->annotation_bed_file);
	my $output_fh = Genome::Sys->open_file_for_writing($self->output_tsv_file);	
	
    unless (-s $cluster_bed) {
        $self->warning_message("Input cluster_bed_file: $cluster_bed is not valid.");
        return 1;
    }
	########### INTERSECTING CLUSTER BED FILE WITH ANNOTATION #############
	
	my @bed_names_array;
	for my $annotation(@annotation)
	{
		
		my (undef, $bed_temp_name) = Genome::Sys->create_temp_file();
   		unless (Genome::Model::Tools::BedTools::Intersect ->execute
  			 	(
   					input_file_a 		=> $cluster_bed,
   					input_file_a_format => 'bed',
   					input_file_b 		=> $annotation,
   					intersection_type	=> 'overlap_both',
   					output_file 		=> $bed_temp_name,
   					
   				)
   			)
   			{die;}			
   		push (@bed_names_array,$bed_temp_name);
    
	}
    ############# PRINTING HEADER TO OUTPUT FILE #############
	print $output_fh join("\t","Cluster#", "Chr", "Start", "Stop");
    for my $track_name (@track_names)
   	{
   		print $output_fh "\t".$track_name ."\t". "Overlap"."\t". "Size of Species";
   	}
   	
	my %hash_of_annotation;
	
	##### ITERATING OVER EACH ANNOTATION INTERSECT BED FILE########
	for my $intersect_bed_file(@bed_names_array)			
	{
		#print "file"."\t".$intersect_bed_file."\n";
		my $annotation_fh = Genome::Sys->open_file_for_reading($intersect_bed_file);  
		while (my $anno_row = $annotation_fh->getline)
				    {
				        chomp $anno_row;  
				        my (@match, @size, @overlap);
    					my @array_ids = split (/\t/,$anno_row); 
    		    		my $species_size = ($array_ids[6] - $array_ids[5]) + 1; 
    		    		my $anno_chr	 = $array_ids[0];
    		   	 		my $anno_start 	 = $array_ids[1];
    		    		my $anno_stop 	 = $array_ids[2];
    		    		my $overlap 	 = $array_ids[8];
    		    		my $species_name = $array_ids[7];
    		    		my $main_cluster_name = $array_ids[3];
    		    	    
    		    	    
    		    	    push (@{$hash_of_annotation{$main_cluster_name}{$intersect_bed_file}{'match'}}, $species_name )	;
    		    	    push (@{$hash_of_annotation{$main_cluster_name}{$intersect_bed_file}{'size'}}, $species_size )	;
    		    	    push (@{$hash_of_annotation{$main_cluster_name}{$intersect_bed_file}{'overlap'}}, $overlap )	;
    		    	    	    	
    		    		#push (@size, $species_size);   	
    		    		#push (@overlap,$overlap )	;    		
    		    	}
	}
	
	
	############# OPENING MAIN CLUSTERS FILE AND WRITING OUTPUT #############
	my $cluster_bed_fh = Genome::Sys->open_file_for_reading($cluster_bed);			
	while (my $row = $cluster_bed_fh->getline) 
	{
    	    chomp $row;
    	    my ( $chr, $start, $stop, $name ) = split( "\t", $row );
    	    print $output_fh "\n".$name."\t".$chr."\t".$start."\t".$stop;
			
			for my $bed_file(@bed_names_array)			
			{
				print $output_fh "\t";		
				if (exists($hash_of_annotation{$name}) )		
				{
					print $output_fh join(',',@{$hash_of_annotation{$name}{$bed_file}{'match'}})."\t".
					join(',',@{$hash_of_annotation{$name}{$bed_file}{'overlap'}})."\t".
					join(',',@{$hash_of_annotation{$name}{$bed_file}{'size'}});
				
				}
			}
	}


return 1;   
}

1;


__END__

#3	52302348	52302376	CLUSTER-1	3	52301034	52304727	WDR82:ENST00000296490:intron:2:rev	28
#3	52302348	52302376	CLUSTER-1	3	52301034	52304727	WDR82:NM_025222.3:intron:2:rev	28
