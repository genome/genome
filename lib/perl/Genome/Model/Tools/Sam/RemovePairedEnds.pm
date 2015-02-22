package Genome::Model::Tools::Sam::RemovePairedEnds;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use Getopt::Long;

class Genome::Model::Tools::Sam::RemovePairedEnds
{
    is => 'Genome::Model::Tools::Sam',
    has_param => [
            lsf_resource =>          
            {
                value => "-R 'span[hosts=1] rusage[mem=10000]' -M 10000000",
            },
    ],
    has_input => [
            sam1            => {
                                    doc         => 'first input file to remove N\'s from',
                                    is          => 'String',
                               },
            sam2            => {
                                    doc         => 'second input file to remove N\'s from',
                                    is          => 'String',
                               },
            paired_end_removed_file1 => {
                                    doc         => 'file to write to',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                   },
            paired_end_removed_file2 => {
                                    doc         => 'file to write to',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                   },
            paired_end_removed_count1 => {
                                    doc         => '# of reads removed from 1st file',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                   },
            paired_end_removed_count2 => {
                                    doc         => '# of reads removed from 2nd file',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                   },
         ],
};

sub help_brief 
{
    return <<"EOS"
    paired end removal - given 2 files, remove reads from each file which are not also present in the other
    NOTE: The input .sam files should be the UNALIGNED output of the bwa contamination screen step
EOS
}

sub help_synopsis 
{
    return <<"EOS"
remove_unpaired_contaminated_reads.pl -sam1 <UNALIGNED .sam file for END_1> -sam2 <UNALIGNED .sam file for END_2> -output <output root name>
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
    my ($sam_file1, $sam_file2) = ($self->sam1, $self->sam2);
    my ($perc1, $perc2) = (0,0);

    #Build %wanted_read hash
    my %wanted_read;
    my $sam1_starting_count;
    my $sam1 = new IO::File $sam_file1;
    while (<$sam1>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        $wanted_read{$line[0]} = 1;
        $sam1_starting_count++;
    }
    $sam1->close;

    #Stream through sam2 file and build sam2 output where all unpaired, contaminated reads are removed
    my $sam2_output_name = $self->paired_end_removed_file2 ? $self->paired_end_removed_file2 : $sam_file2 . ".PE_contam_free.END2.sam";
    my $sam2_output = new IO::File ">>$sam2_output_name";
    my $sam2 = new IO::File $sam_file2;
    my $sam2_final_count;
    while (<$sam2>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        if (defined $wanted_read{$line[0]}) 
        {
	    $wanted_read{$line[0]}++;
	    print $sam2_output "$line\n";
	    $sam2_final_count++;
        }
        else
        {
            $perc2++;
        }
    }
    $sam2->close;
    $sam2_output->close;

    #Finally stream through sam1 file again and build sam1 output where all unpaired, contaminted reads are removed
    my $sam1_output_name = $self->paired_end_removed_file1 ? $self->paired_end_removed_file1 : $sam_file1 . ".PE_contam_free.END1.sam";
    my $sam1_output = new IO::File ">>$sam1_output_name";
    $sam1 = new IO::File $sam_file1;
    my $sam1_final_count;
    while (<$sam1>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        my @line = split(/\s+/,$line);
        if ($wanted_read{$line[0]} == 2) {
	    print $sam1_output "$line\n";
	    $sam1_final_count++;
        }
        else
        {
            $perc1++;
        }
    }
    $sam1->close;
    $sam1_output->close;

    $self->paired_end_removed_file1($sam1_output_name);
    $self->paired_end_removed_file2($sam2_output_name);

    $self->paired_end_removed_count1($perc1);
    $self->paired_end_removed_count2($perc2);

    return 1;
}

1;
