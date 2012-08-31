
package Genome::Model::Tools::BioDbFasta::Subquality;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Subquality.pm - 	Retrieve quality values from a fasta file in a a Bio::DB::Fasta indexed directory
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
use Carp;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::BioDbFasta::Subquality {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
        dir      => { is => 'Text',       doc => "Indexed directory of quality files" },
        name      => { is => 'Text',       doc => "Name of sequence to extract" },        
        start      => { is => 'Text',       doc => "Sequence start position" },
        stop      => { is => 'Text',       doc => "Sequence stop position" },        
    ], 
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Extract quality values by sequence name, start, and stop"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
This command extracts one or more base qualities using the Bio::DB::Fasta index of a directory containing QUALITY files
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
The name provided should match a record name from one of the FASTA files in the directory provided.
For example, consider the following QUAL entry:
 >read1
 15 15 20 28 28 29 31 31 28 28 24 22 19 14 7

The following command would extract a sub-sequence of "read1":
 gmt bio-db-fasta subquality --dir mypath/qual_dir/ --name read1 --start 6 --stop 10
 29 31 31 28 28
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $dir = $self->dir;

	## Verify that directory exists and is indexed ##

	if(-d $dir)
	{
		if(!(-e $dir . '/directory.index'))		
		{
			print "Input directory index not yet built (run gmt bio-db-fasta build). Exiting...\n";
			return(0);		
		}
	}
	else
	{
		print "Input directory does not exist. Exiting...\n";
		return(0);			
	}

	## Access the index and retrieve the sequence ##
	my $seqdb = Bio::DB::Fasta->new($dir) or 
            croak "Bio::DB::Fasta call failed!?"; 
	my $qpiece = $seqdb->seq($self->name, $self->start => $self->stop);

	## Convert to numeric ##
	
	my $qpiece_num = "";
	
	for(my $baseCounter = 0; $baseCounter < length($qpiece); $baseCounter++)
	{
		my $Qchar = substr($qpiece, $baseCounter, 1);
		my $Qnum = ord($Qchar) - 33;
		$qpiece_num .= " " if($qpiece_num);
		$qpiece_num .= $Qnum;
	}

	## Print it ##
	
	#print "$qpiece_num\n" if($qpiece && length($qpiece) > 0);
	print "$qpiece_num\n" if($qpiece);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

