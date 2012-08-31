package Genome::Model::Tools::Bwa::UnalignedDataToFastq;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bwa::UnalignedDataToFastq {
    is => 'Command',
    has => [
        in              => { is => 'Text', 
                            doc => "Pathname to a file generated from samtools sam(p|s)e option", 
                            shell_args_position => 1 },
    ],
    has_optional => [
        fastq           => { is => 'Text', 
                            doc => 'the output pathname for "forward" reads (or all reads on a fragment run)', 
                            shell_args_position => 2 },

        reverse_fastq   => { is => 'Text', 
                            doc => 'the output pathname for "reverse" for paired-end data', 
                            shell_args_position => 3 },
    ],
    doc => "Create a new fastq-format file containing reads that aligned poorly in the prior align-reads step"
};

sub help_synopsis {
    return <<"EOS"
    gmt bwa unaligned-data-to-fastq in.bwa-unaligned out.fwd.fastq out.rev.fastq
    
    gmt bwa unaligned-data-to-fastq -i /path/to/inputpathname.unaligned -f /path/to/forward.fastq -r /path/to/reverse.fastq
EOS
}

sub help_detail {                           
    return <<EOS 
As part of the aligmnent process, the -u option to maq will create a file containing 
reads that had low alignment quailty, or did not align at all.  This command
will take that file and create a fastq-formatted file containing those reads.

For paired-end data, this tool will make 2 files, and expect 2 outputs to be specified.

It will correct the error common in some maq unaligned files that read #2 is misnamed with the name of its mate.
EOS
}

sub execute {
    my $self = shift;
    
$DB::single = $DB::stopper;

    my $unaligned_file = $self->in();
    my $unaligned = IO::File->new($unaligned_file);
    unless ($unaligned) {
        $self->error_message("Unable to open $unaligned_file for reading: $!");
        return;
    }

    my $unaligned_fastq_file1 = $self->fastq();
    my $fastq1 = IO::File->new(">$unaligned_fastq_file1");
    unless ($fastq1) {
        $self->error_message("Unable to open $unaligned_fastq_file1 for writing: $!");
        return;
    }
  
    my $unaligned_fastq_file2 = $self->reverse_fastq();
    my $fastq2;

    if ($unaligned_fastq_file2) {
        $fastq2 = IO::File->new(">$unaligned_fastq_file2");
        unless ($fastq2) {
            $self->error_message("Unable to open $unaligned_fastq_file2 for writing: $!");
            return;
        }
    }

    my ($read_name,$alignment_quality,$sequence,$read_quality);
    my $last_read_name;
    my $warned;
    if ($fastq2) {
        while(<$unaligned>) {
            
            # bwa seems to lop off the strand identifier.  assuming the first input seen is the forward strand
            # and the second listing is the reverse.
            
            chomp;
            $fastq1->print($self->_chunk_sam_into_fastq($_,1));            
            my $rev = <$unaligned>;
            chomp $rev;
            $fastq2->print($self->_chunk_sam_into_fastq($rev,2));    
        }
        
    }
    else {
        while(<$unaligned>) {
            chomp;
            $fastq1->print($self->_chunk_sam_into_fastq($_));
        }
    }

    $unaligned->close();
    $fastq1->close();
    $fastq2->close() if $fastq2;
    return 1;
}

sub _chunk_sam_into_fastq {
    my ($self, $in, $strand) = @_;
        
    my @components = split /\s+/, $in;
    
    my ($read_name,$sequence,$read_quality) = ($components[0], $components[9], $components[10]);
    if ($strand) {
        $read_name .= "/" . $strand;
    }
    
    return "\@$read_name\n$sequence\n\+$read_name\n$read_quality\n";
}

1;

