package VarScan::SequenceFile;

use warnings;
use strict;

=head1 NAME

VarScan::SequenceFile - Functions for working with FASTA and QUALITY files

=head1 VERSION

    Version 1.03

=cut

our $VERSION = '1.03';

=head1 SYNOPSIS

    This module contains functions for working with sequence and quality files in FASTA format.

=head1 FUNCTIONS


=head2 build_bio_db_fasta - a module to index a directory for use with Bio::DB::Fasta

=cut

sub build_bio_db_fasta
{
    my $dir = shift(@_);
    
    if(-d $dir)
    {
        print "Creating Bio::DB::Fasta object for $dir...\n";
	my $seq_db = Bio::DB::Fasta->new($dir) or die "!$\n";
        return($seq_db);
    }
    
    return(0);
}


=head2 build_bio_db_qual - a module to index a directory for use with Bio::DB::Qual

=cut

sub build_bio_db_qual
{
    my $dir = shift(@_);
    
    if(-d $dir)
    {
        print "Creating Bio::DB::Qual object for $dir...\n";

	my $seq_db = Bio::DB::Qual->new($dir) or die "!$\n";

        return($seq_db);
    }
    
    return(0);
}

=head2 quality_to_fasta - Converts a quality file with Phred-like scores to single-character FASTA format

=cut

sub quality_to_fasta
{
        (my $infile, my $outfile) = @_;
        
	my $input = new FileHandle ($infile);
	
	open(OUTFILE, ">" . $outfile) or die "Can't open outfile: $!\n";
	
	my $lineCounter = 0;
	my $record_name = my $record_seq = "";
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	

			if(substr($line, 0, 1) eq ">")
			{
				if($record_name && $record_seq)
				{
#					printQualString($record_name, $record_seq);
					my $qual_string = getQualString($record_seq);
					print OUTFILE ">$record_name\n$qual_string\n";				
				}				
				
				$record_name = $record_seq = "";
				
				$record_name = substr($line, 1, 999);	
			}
			else
			{
				$record_seq .= " " if($record_seq);				
				$record_seq .= $line;
			}

	}

	if($record_name && $record_seq)
	{
#		printQualString($record_name, $record_seq);
		my $qual_string = getQualString($record_seq);
		print OUTFILE ">$record_name\n$qual_string\n";
	}
	
	close(OUTFILE);
}



#############################################################
=head2 getQualString - convert numeric to ASCII qual values
=cut
#
#############################################################

sub getQualString
{
    my $record_seq = shift;
	
	my @qualValues = split(/\s+/, $record_seq);
	my $numValues = @qualValues;
	
	my $record_squal = "";
	
	for(my $vCounter = 0; $vCounter < $numValues; $vCounter++)
	{
		my $Q = $qualValues[$vCounter];
		$Q-- if($Q == 29); ## Fix the problematic value corresponding to '>' ##
		my $Qchar = chr(($Q<=93? $Q : 93) + 33);
		$record_squal .= "$Qchar";
	}
	
	return($record_squal);
}






=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Center at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use.

=cut

1; # End of VarScan::SequenceFile
