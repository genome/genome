package Genome::Model::Tools::Fasta::RemoveN;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Bio::SeqIO;

class Genome::Model::Tools::Fasta::RemoveN
{
    is => 'Genome::Model::Tools::Fasta',
    has_input => [
            n_removed_file => {
                                    doc => 'file to write to',
                                    is => 'Text',
                                    is_output => 1,
                                    is_optional => 1,
                                },
            cutoff =>   {
                                    doc => 'minimum # of N\'s to screen on.  Set to 0 to disable',
                                    is => 'Number',
                                    is_optional => 1,
                                    default => 10, 
                        },
            save_screened_reads => 
                                {
                                    doc => 'save screened reads in separate file',
                                    is => 'Boolean',
                                    is_optional => 1,
                                    default => 1,
                                },
            screened_count => 
                                {
                                    doc => '# of reads screened',
                                    is => 'Integer',
                                    is_optional => 1,
                                    is_output => 1,
                                },
         ],
};

sub help_brief 
{
    "remove reads from file containing N";
}

sub help_detail
{   
    "Removes reads that have internal N's, or more than cutoff amount of N's on ends.  By default, removes for a single N.  Set cutoff to 0 to disable";
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
    my $fasta_file = $self->fasta_file;
    my $n_removed_file = ($self->n_removed_file ? $self->n_removed_file : $fasta_file . ".n_removed");
    my $screened_file = ($n_removed_file=~m/(.*.)(\..*)/ and "$1.SCREENED.$2");
    my $cutoff = $self->cutoff;
    my @screened;

    my $input_fh = IO::File->new($fasta_file);
    unless ($input_fh) 
    {
        $self->error_message("Failed to open input file " . $fasta_file . ": $!");
        return;
    }

    my $output_fh = IO::File->new('>'.$n_removed_file);
    unless ($output_fh) 
    {
        $self->error_message("Failed to open output file " . $n_removed_file . ": $!");
        return;
    }

    my ($header, $seq, $count, $screened_count) = ("", "", 0, 0);
    while (my $line = $input_fh->getline) 
    {
        if ($line=~/^>.*/) #found a header
        {
            $header = $line;
        }
        else    
        {
            $seq .= $line;
        }
        
        while ($line = $input_fh->getline) #accumulate lines for read, until next header encountered
        {
            if ($line=~/^>.*/) #found a new header - read has been accumulated 
            {
                last;
            }
            else
            {
                $seq .= $line;
            }
        }

        $seq=~s/(N)/$count++;$1/eg; #get N-count
        if ($cutoff > 0 and $count >= $cutoff and ++$screened_count) # check if cutoff disabled, then compare count
        {
            push (@screened, "$header$seq");
        }
        else
        {
            $output_fh->print("$header$seq");
        }
        
        #reset
        $count = 0;
        $seq = '';
        $header = $line;
     }

    $input_fh->close;
    $output_fh->close;

    if ($self->save_screened_reads)
    {
        my $screened_name = "$n_removed_file.SCREENED";
        my $screened_fh = IO::File->new('>'.$screened_name);

        unless ($screened_fh) 
        {
            $self->error_message("Failed to open output file " . $screened_name . ": $!");
            return;
        }
        $screened_fh->print(@screened);
        $screened_fh->close;
    }

    $self->screened_count($screened_count);

    return 1;
}

1;
