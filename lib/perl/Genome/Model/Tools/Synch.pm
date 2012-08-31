package Genome::Model::Tools::Synch;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Bio::SeqIO;

class Genome::Model::Tools::Synch
{
    is => 'Command',
    has_input => [
            name1=> {
                                    doc => 'first sam file',
                                    is => 'String',
                                },
            name2=> {
                                    doc => 'second sam file',
                                    is => 'String',
                                },
            fastq1=> {
                                    doc => 'first fastq file',
                                    is => 'String',
                                },
            fastq2=> {
                                    doc => 'second fastq file',
                                    is => 'String',
                                },
            output=> {
                                    doc => 'file to write to',
                                    is => 'String',
                                    is_output => 1,
                                    is_optional => 1,
                                },
            prefix1=> {
                                    doc => 'file to write to',
                                    is => 'String',
                                    is_output => 1,
                                    is_optional => 1,
                                },
            prefix2=> {
                                    doc => 'file to write to',
                                    is => 'String',
                                    is_output => 1,
                                    is_optional => 1,
                                },
            output_only1 => {
                                doc         => 'first single-end output',
                                is          => 'String',
                                is_output   => 1,
                                    is_optional => 1,
                            },
            output_only2 => {
                                doc         => 'second single-end output',
                                is          => 'String',
                                is_output   => 1,
                                    is_optional => 1,
                            },
            output_both1 => {
                                doc         => 'first paired-end output',
                                is          => 'String',
                                is_output   => 1,
                                    is_optional => 1,
                            },
            output_both2 => {
                                doc         => 'second paired-end output',
                                is          => 'String',
                                is_output   => 1,
                                    is_optional => 1,
                            },
         ],
};

sub help_brief 
{
    "remove reads from file containing N";
}

sub help_detail
{   
    "Removes reads that have internal N's, or more than cutoff amount of N's on ends.  Cutoff is 6 for 75-mer, 9 for 100-mer";
}

sub help_synopsis 
{
    return <<"EOS"
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
    my %ref_hash;
    my $output = $self->output;

    ####print "Reading the first fasta file into the reference hash......\n";

    my $name1  = new IO::File $self->name1;
    while (<$name1>) 
    {
	chomp;
	next if ( $_ =~ /^\s$/ );
	my @line = split(/\s+/);
	$ref_hash{ $line[0] } = 1;
    }
    $name1->close;

    ####my $size = keys(%ref_hash);
    ####print "Size of first fasta file is $size\n";

    ####print "Reading the second fasta file into the reference hash......\n";
    my $name2_counter=0;

    my $name2  = new IO::File $self->name2;
    while (<$name2>) 
    {
        chomp;
        next if ( $_ =~ /^\s$/ );
        my @line = split(/\s+/);
        $name2_counter++;
        if ( defined $ref_hash{ $line[0] } ) 
        {
	    $ref_hash{ $line[0] } = 3;
        } 
        else 
        {
	    $ref_hash{ $line[0] } = 2;
        }
    }
    $name2->close;

    ####print "Size of second fasta file is $name2_counter\n";
    ####$size = keys(%ref_hash);
    ####print "Size of ref hash is $size\n";

    #At this point I want to stream through each fastq file in linear fashion,
    #and shuffle each read into its appropriate location.  If we are looking at
    #pair end /1, then a $ref_hash{<name>} value of 1 means this is an
    #END_1 singlet, a 3 means its a PAIR, and a 2 means that its an END_2
    #singlet and must be DISCARDED from end 1 (i.e. not entered into any
    #output file).
    #
    #If we are looking at pair end /2, then a $ref_hash{<name>} value of 1 
    #means this is an END_1 singlet and must be DISCARDED from end 2(i.e. not
    #entered into any output file), a 3 means its a PAIR, and a 2 means that 
    #its an END_2 singlet.

    #___END_1
    my $fastq1 = new IO::File $self->fastq1;
    my $output_only_file1_name = $self->prefix1 . $output . ".only_end1.fastq";
    my $output_both_file1_name = $self->prefix1 . $output . ".paired_end1.fastq";
    my $output_only_file1 = new IO::File ">>$output_only_file1_name";
    my $output_both_file1 = new IO::File ">>$output_both_file1_name";
    my $name;
    while (<$fastq1>) 
    { #I will pull in 4 lines at a time to get the full fastq record
       my $header = $_;
        chomp($header);
        my $sequence = <$fastq1>; #peel another line out of $fastq1
        chomp($sequence);
        my $qual_header = <$fastq1>; #peel another line out of $fastq1
        chomp($qual_header);
        my $qual_sequence = <$fastq1>; #peel another line out of $fastq1
        chomp($qual_sequence);
        $name = $header;
        $name =~ s/^\@//;
        $name =~ s/\s+$//;
        $name =~ s/\/[12]$//;
        if (defined $ref_hash{$name}) 
        {
	    if ($ref_hash{$name} == 1) 
            {
	        print $output_only_file1 "$header\n$sequence\n$qual_header\n$qual_sequence\n";
            } 
            elsif ($ref_hash{$name} == 3) 
            {
	        print $output_both_file1 "$header\n$sequence\n$qual_header\n$qual_sequence\n";
	    } 
            elsif ($ref_hash{$name} == 2) 
            {
	        my $var = 1; #This guy is present in the END_2 list, but NOT in the END_1 list, so skip this read in the fastq file here
	    } 
            else 
            {
	        die "unexpected state in ref_hash ($ref_hash{$name}) for END_1 read =>$header<=\n";
    	    }
        } #if the fastq root name is not in the $ref_hash, do nothing with the record
    }
    $output_only_file1->close;
    $output_both_file1->close;
    $fastq1->close;

    #___END_2
    my $fastq2 = new IO::File $self->fastq2;
    my $output_only_file2_name = $self->prefix2 . $output . ".only_end2.fastq";
    my $output_both_file2_name = $self->prefix2 . $output . ".paired_end2.fastq";
    my $output_only_file2 = new IO::File ">>$output_only_file2_name";
    my $output_both_file2 = new IO::File ">>$output_both_file2_name";
    while (<$fastq2>) 
    { #I will pull in 4 lines at a time to get the full fastq record
        my $header = $_;
        chomp($header);
        my $sequence = <$fastq2>; #peel another line out of $fastq2
        chomp($sequence);
        my $qual_header = <$fastq2>; #peel another line out of $fastq2
        chomp($qual_header);
        my $qual_sequence = <$fastq2>; #peel another line out of $fastq2
        chomp($qual_sequence);
        $name = $header;
        $name =~ s/^\@//;
        $name =~ s/\s+$//;
        $name =~ s/\/[12]$//;
        if (defined $ref_hash{$name}) 
        {
	    if ($ref_hash{$name} == 2) 
            {
	        print $output_only_file2 "$header\n$sequence\n$qual_header\n$qual_sequence\n";
	    } 
            elsif ($ref_hash{$name} == 3) 
            {
	        print $output_both_file2 "$header\n$sequence\n$qual_header\n$qual_sequence\n";
	    } 
            elsif ($ref_hash{$name} == 1) 
            {
	        my $var = 1; #This guy is present in the END_1 list, but NOT in the END_2 list, so skip this read in the fastq file here
	    } 
            else 
            {
	        die "unexpected state in ref_hash ($ref_hash{$name}) for END_2 read =>$header<=\n";
	    }
        } #if the fastq root name is not in the $ref_hash, do nothing with the record
    }
    $output_only_file2->close;
    $output_both_file2->close;
    $fastq2->close;

    $self->output_only1($output_only_file1_name);
    $self->output_both1($output_both_file1_name);
    $self->output_only2($output_only_file2_name);
    $self->output_both2($output_both_file2_name);

    return 1;
}

1;
