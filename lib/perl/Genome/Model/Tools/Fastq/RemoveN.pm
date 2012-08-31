package Genome::Model::Tools::Fastq::RemoveN;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Bio::SeqIO;

class Genome::Model::Tools::Fastq::RemoveN
{
    is => 'Genome::Model::Tools::Fastq',
    has_input => [
        n_removed_file => {
            doc => 'file to write to',
            is => 'Text',
            is_output => 1,
            is_optional => 1,
        }, 

        n_removal_threshold =>   {
            doc => 'minimum # of N\'s to screen on.  Set to 0 to disable',
            is => 'Number',
            is_optional => 1,
        },
        non_n_base_threshold =>   {
            doc => 'minimum # of consecutive non N\'s to screen on. Set to 0 to disable',
            is => 'Number',
            is_optional => 1,

        },
        save_screened_reads => 
        {
            doc => 'save screened reads in separate file',
            is => 'Boolean',
            is_optional => 1,
            default => 0,
        },
    ],
    has_output => [
        passed_read_count => {
            is=>'Number',
            doc => 'number of reads passed screening',
            is_optional => 1
        },
        failed_read_count => {
            is=>'Number',
            doc => 'number of reads failed screening',
            is_optional => 1
        },
    ]
};

sub help_brief 
{
    "Remove reads from file containing N or remove reads without expected number of non-N bases";
}

sub help_detail
{   
    "If N removal cut=off is set, removes reads that have internal N's, or more than cutoff amount of N's on ends.  By default, removes for a single N.  Set cutoff to 0 to disable. If the non-N threshold is set, it removed reads that have lesser than the desired number of non-N bases.  "; 
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
    my $fastq_file = $self->fastq_file;
    my $n_removed_file = ($self->n_removed_file ? $self->n_removed_file : $fastq_file . "n_removed");
    my $save_screened_reads = $self->save_screened_reads;

    if (!$self->n_removal_threshold && !$self->non_n_base_threshold) {
        $self->error_message("Need one threshold set.");
        return;
    }

    my $input_fh = IO::File->new($fastq_file);
    unless ($input_fh) {
        $self->error_message("Failed to open input file " . $fastq_file . ": $!");
        return;
    }

    my $output_fh = IO::File->new('>'.$n_removed_file);
    unless ($output_fh) {
        $self->error_message("Failed to open output file " . $n_removed_file . ": $!");
        return;
    }

    my $passed_reads = 0;
    my $failed_reads = 0;

    while (my $header = $input_fh->getline) 
    {
        my $seq = $input_fh->getline;
        my $sep = $input_fh->getline;
        my $qual = $input_fh->getline;
        my $count = 0;
        my $cutoff;
        if ($self->n_removal_threshold){
            $cutoff=$self->n_removal_threshold;
            $seq=~s/(N)/$count++;$1/eg; # get N-count
            if($cutoff > 0 and $count >= $cutoff) {
                $failed_reads++;
            }else {
                $passed_reads++;   
                $output_fh->print("$header$seq$sep$qual");
            }
        }   
        elsif ($self->non_n_base_threshold){
            $cutoff=$self->non_n_base_threshold;
            $seq=~s/([AGCTagct])/$count++;$1/eg; # get non=-N-count
            if ($cutoff > 0 and $count < $cutoff) {
                $failed_reads++;
            } else {
                $passed_reads++;   
                $output_fh->print("$header$seq$sep$qual");
            }
        }
    }
    $input_fh->close;
    $output_fh->close;

    $self->passed_read_count($passed_reads);
    $self->failed_read_count($failed_reads);

    return 1;
}

1;
