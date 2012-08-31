
package Genome::Model::Tools::BioDbFasta::Subsequence;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Subsequence.pm - 	Retrieve subsequence from a fasta file in a a Bio::DB::Fasta indexed directory
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
use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::BioDbFasta::Subsequence {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
        dir      => { is => 'Text',       doc => "Indexed directory of sequence files" },
        name      => { is => 'Text',       doc => "Name of sequence to extract" },        
        start      => { is => 'Text',       doc => "Sequence start position" },
        stop      => { is => 'Text',       doc => "Sequence stop position" },        
    ], 
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Extract subsequence by name, start, and stop"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
This command extracts a sub-sequence using the Bio::DB::Fasta index of a directory containing FASTA files
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
The name provided should match a record name from one of the FASTA files in the directory provided.
For example, consider the following FASTA entry:
 >read1
 AAAAAAAGGGGGGGCCTTTTTTTTTTTAAAAAAAACCCCCCCCC
 
The following command would extract a sub-sequence of "read1":
 gmt bio-db-fasta subsequence --dir mypath/fasta_dir/ --name read1 --start 6 --stop 10
 AAGGG

Note: To extract sequence from the human Hs36 refseq, use "hs36" as your directory:
 gmt bio-db-fasta subsequence --dir hs36 --name chr2 --start 1704463 --stop 1704468
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $dir = $self->dir;

	## If requested, point to the human refseq ##
	
	if($dir eq "hs36")
	{
		$dir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
	}

	## Verify that directory exists and is indexed ##

	if(-d $dir)
	{
		if(!(-e $dir . '/directory.index'))		
		{
                    $self->error_message('Input directory index not yet built (run gmt bio-db-fasta build). Exiting...');
                    return;
		}
	}
	else
	{
            $self->error_message('Input directory does not exist. Exiting...');
            return(0);
	}

	## Access the index and retrieve the sequence ##
	my $seqdb = Bio::DB::Fasta->new($dir); 
        unless(defined($seqdb))
        {
            $self->error_message("Failed to create new Bio::DB::Fasta object for dir '$dir'");
            return;
        }
	my $spiece = $seqdb->seq($self->name, $self->start => $self->stop);

	print "$spiece\n" if($spiece);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

