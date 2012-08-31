
package Genome::Model::Tools::BioDbFasta::Build;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Build.pm - 		Build a Bio::DB::Fasta index for a directory containing FASTA files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/15/2008 by D.K.
#	MODIFIED:	10/16/2008 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use Bio::DB::Fasta;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::BioDbFasta::Build {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
        dir      => { is => 'Text',       doc => "directory containing FASTA files" },
        rebuild      => { is => 'Boolean',    doc => "if set to 1, rebuild the directory index", is_optional => 1 },
    ], 
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Index a directory for Bio::DB::Fasta access"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
This command indexes a given directory for fast random sequence access with Bio::DB::Fasta
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
EXAMPLE 1: Index all FASTA files in the directory mypath/fasta_dir/
 gmt bio-db-fasta build --dir mypath/fasta_dir/

EXAMPLE 2: Force reindexing of the same directory
 gmt bio-db-fasta build --dir mypath/fasta_dir/ --rebuild 1
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;
	
	if(!(-d $self->dir))
	{
		print "Input directory does not exist. Exiting...\n";
		return(0);
	}
	
	if($self->rebuild)
	{
		print "Rebuilding Bio::DB::Fasta index for " . $self->dir . "\n";	
		my $seq_db = Bio::DB::Fasta->new($self->dir);	
	}
	else
	{
		print "Building Bio::DB::Fasta index for " . $self->dir . "\n";
		my $seq_db = Bio::DB::Fasta->new($self->dir);	
	}
	

    
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;

