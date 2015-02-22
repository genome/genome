package Genome::Model::Tools::Sam::ResurrectDuplicates;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use Getopt::Long;

class Genome::Model::Tools::Sam::ResurrectDuplicates
{
    is => 'Genome::Model::Tools::Sam',
    has_param => [
            lsf_resource =>          
            {
                                    default_value => "-M 15000000 -R 'rusage[mem=15000]'",
            },
    ],
    has_input => [
            raw_dedup_sam1=> {
                                    doc         => 'first input file',
                                    is          => 'String',
                               },
            raw_dedup_sam2=> {
                                    doc         => 'second input file',
                                    is          => 'String',
                               },
            orig_sam1=> {
                                    doc         => 'first original file',
                                    is          => 'String',
                                   },
            orig_sam2=> {
                                    doc         => 'second original file',
                                    is          => 'String',
                                   },
            output_file1 =>      {
                                    doc         => 'first output file',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                },
            output_file2 =>      {
                                    doc         => 'second output file',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                },
         ],
};

sub help_brief 
{
    return <<"EOS"
NOTE: The original .sam files MUST EACH CONTAIN the full range of reads that will be requested when looking at the union of the two read sets.  If they do not this script will die with a message, but the 
EOS
}

sub help_detail
{
    return <<"EOS"
    This script will take 2 .sam files that are the product of 
    the deduplication steps, and 2 original .sam files (that would
    have all the reads), and it will build a merged, unique list of
    all reads in both the raw, dedup sam files.  Then it will pull
    both end versions of all those root names and build new, resurrected
    dedup sam files.
EOS
}

sub help_synopsis 
{
    return <<"EOS"
resurrect_duplicates.pl -raw_dedup_sam1 <.sam file resulting from deduplication of end1> -raw_dedup_sam2 <.sam file resulting from deduplication of end2> -orig_sam1 <.sam file that was the input of dedup of end1> -orig_sam2 <.sam file that was the input of dedup of end2> -output <output root name>
EOS
}

sub create 
{
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute 
{
    my $self = shift;

    my ($raw_dedup_sam1,       $raw_dedup_sam2,         $orig_sam1_file,    $orig_sam2_file) = 
       ($self->raw_dedup_sam1, $self->raw_dedup_sam2,   $self->orig_sam1,   $self->orig_sam2); 

    my ($output_file1, $output_file2) = ($self->output_file1 ?
                                         $self->output_file1 :
                                         $orig_sam1_file. 'resurrected_end1.sam', 
                                         $self->output_file2 ?
                                         $self->output_file2 :
                                         $orig_sam2_file. 'resurrected_end2.sam'); 

    print "Resurrect begun with $orig_sam1_file and $orig_sam2_file\n";

    #Build master hash holding the rootname of all wanted reads
    #This will be a unique listing of all reads for BOTH ends of
    #the raw, deduplicated .sam files.
    my %master_hash;
    my $dedup1 = new IO::File $raw_dedup_sam1;
    while (<$dedup1>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        $master_hash{$line[0]} = 1;
    }
    $dedup1->close;
print "attempting $raw_dedup_sam2\n";
    my $dedup2 = new IO::File $raw_dedup_sam2;
    while (<$dedup2>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        $master_hash{$line[0]} = 1;
    }
    $dedup2->close;

    #Now that I have the hash of all wanted reads that should be represented
    #in both end-specific, resurrected deduplication .sam files, I will
    #stream through each of the original, end-specific .sam files and pull
    #these guys out into the end-specific output files
    my $orig_sam1 = new IO::File $orig_sam1_file;
    my $end1_output = new IO::File ">>$output_file1";
    while (<$orig_sam1>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        if (defined $master_hash{$line[0]}) 
        {
	    my $header = $line[0] . "/1";
       	    print $end1_output ">$header\n";
	    print $end1_output "$line[9]\n";
	    $master_hash{$line[0]}++; #This is to make sure every read in the master hash is found
        }
    }
    $orig_sam1->close;
    $end1_output->close;
    my $orig_sam2 = new IO::File $orig_sam2_file;
    my $end2_output = new IO::File ">>$output_file2";
    while (<$orig_sam2>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        if (defined $master_hash{$line[0]}) 
        {
	    my $header = $line[0] . "/2";
	    print $end2_output ">$header\n";
	    print $end2_output "$line[9]\n";
	    $master_hash{$line[0]}++; #This is to make sure every read in the master hash is found
        }
    }
    $orig_sam2->close;
    $end2_output->close;

    #Final sanity check to make sure every read in the master hash was found for both ends
    foreach my $name (keys(%master_hash)) 
    {
        unless ($master_hash{$name} == 3) 
        {
	    warn "problem detected...removing the resurrected sam files created by this script\n";
	    unlink("$output_file1");
	    unlink("$output_file2");
	    die "did not find =>$name<= in at least one of the original sam files (this should never happen)\n";
        }
    }

    print "Resurrect completed with $orig_sam1_file and $orig_sam2_file\n";
    return 1;
}

1;
