package Genome::Model::Tools::HmpShotgun::SpeciesBlastx;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::SpeciesBlastx {
    is  => ['Command'],
    has => [
        model_id => {
            is  => 'String',
            is_input => '1',
            doc => 'The model id to process.',
        },
        working_directory => {
            is => 'String',
            is_input => '1',
            doc => 'The working directory where results will be deposited.',
        },
        query_file => {
            is => 'String',
            is_input => '1',
            doc => 'The reads to query.',
        },
        blastx_database => {
        	is => 'String',
            is_input => '1',
            doc => 'The database to query against.',
        }, 
        blastx_jobs => {
        	is => 'Integer',
            is_input => '1',
            doc => 'The number of jobs to spread the blastx query across.',
        },
        final_file => {
            is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The model id to process.',
        },
        delete_intermediates => {
            is => 'Integer',
            is_input =>1,
            is_optional =>1,
            default=>0,
        },
    ],
};


sub help_brief {
    'Use Blastx to align against 800 proteomes';
}

sub help_detail {
    return <<EOS
    Use Blastx to align against 800 proteomes
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

	#gmt wu-blast blastx parallel 
	
	#--database=/gscmnt/temp206/info/seqana/species_independant/sabubuck/KEGG/KEGG-52/genesV52.KO.faa 
	#--query-file=/gscmnt/temp206/info/seqana/species_independant/jpeck/testdata/s_2_1.small.fasta 
	#--params="wordmask=seg W=4 T=12 B=1 V=1" 
	#--jobs=100 --output-file=blast_off.t12.out    
	
	#This one works!
	#specify an output directory not an output file!
	#gmt wu-blast blastx parallel 
	#--database=/gscmnt/temp206/info/seqana/species_independant/sabubuck/KEGG/KEGG-52/genesV52.KO.faa 
	#--query-file=/gscmnt/temp206/info/seqana/species_independant/jpeck/testdata/s_2_1.tiny.fasta 
	#--jobs=2 
	#--params="wordmask=seg W=4 T=12 B=1 V=1" 
	#--output-directory=/gscmnt/temp206/info/seqana/species_independant/jpeck/tmp/blastx_output

    my $now = UR::Time->now;
    $self->status_message(">>>Starting SpeciesBlastX execute() at $now"); 
    $self->status_message("Bailing for alignment testing...");
    $self->final_file("species_blastx_final_file_path");
    $self->status_message("<<<Completed SpeciesBlastX at ".UR::Time->now);
    return 1;
    
    
    #my $fail_file = $self->working_directory."/blastx/FAIL";
    #if (-e $fail_file) {
    #	die "Found a FAIL file.  Quitting.";
    #} else {
    #	Genome::Sys->shellcmd(cmd=>"touch $fail_file" );
    #}
    
    #my $database="/gscmnt/temp206/info/seqana/species_independant/sabubuck/KEGG/KEGG-52/genesV52.KO.faa";
    my $database= $self->blastx_database;
    my $query_file = $self->query_file;
    
    #for testing.
    #$query_file = "/gscmnt/temp206/info/seqana/species_independant/jpeck/testdata/s_2_1.tiny.fasta";
    
    my $params_string = "wordmask=seg W=4 T=12 B=1 V=1";
    my $jobs = $self->blastx_jobs;
    
    my $output_directory = $self->working_directory."/blastx/";
    my $query_file_name = File::Basename::basename($query_file);
    my $output_file = $output_directory."/".$query_file_name.".blast";
    
    my  @expected_output_files = ($output_file);
    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
    if (defined($rv_check)) {
	    if ($rv_check == 1) {
	    	#shortcut this step, all the required files exist.  Quit.
	    	$self->status_message("Skipping this step.  If you would like to regenerate these files, remove them and rerun.");
	   	    $self->status_message("<<<Completed SpeciesBlastX at ".UR::Time->now);
	   	    return 1;
	    }
	}
    
    my $blaster = Genome::Model::Tools::WuBlast::Blastx::Parallel->create(database=>$database,
    																	  query_file=>$query_file,
    																	  jobs=>$jobs,
    																	  output_directory=>$output_directory,
    																	  params=>$params_string,);
    																	  
    my $rv_blast = $blaster->execute;
    if ($rv_blast != 1) {	
    	$self->error_message("Blastx error.  Return value: $rv_blast");
    	die "SpeciesBlastx could not complete.";
    } 							
    					  
    Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    
    $self->final_file("species_blastx_final_file_path");
    $self->status_message("<<<Completed species_blastx execute() at ".UR::Time->now); 
    return 1;

}
1;
