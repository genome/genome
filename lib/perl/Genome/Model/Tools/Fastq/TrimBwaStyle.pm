package Genome::Model::Tools::Fastq::TrimBwaStyle;

#Adapting from TrimBwaStyle.pl

use strict;
use warnings;

use Genome;

use File::Basename;
require Genome::Model::Tools::Sx::FastqReader;
require Genome::Model::Tools::Sx::FastqWriter;

class Genome::Model::Tools::Fastq::TrimBwaStyle {
    is  => 'Command',
    has_input => [
        fastq_file  => {
            is  => 'Text',
            doc => 'the input fastq file path',
            is_optional => 1,
        }, 
        out_file    => {
            is  => 'Text',
            doc => 'the file path of the output file, default is xxx.trimmed.fastq in fastq_file dir',
            is_optional => 1,
        },
        trim_qual_level => {
            is  => 'Integer',
            doc => 'trim quality level',
            default => 10,
            is_optional => 1,
        },
        qual_type   => {
            is  => 'Text',
            doc => 'The fastq quality type, must be either sanger(Qphred+33) or illumina(Qphred+64)',
            valid_values  => ['sanger', 'illumina'],
            default_value => 'sanger',
            is_optional   => 1,
        },
        _qual_str => {
            is => 'Text',
        },
        _qual_thresh => {
            is => 'Number',
        },
        report_file => {
            is  => 'Text',
            doc => 'the file path of the trim report file, default is trim.report in the same dir as out_file',
            is_optional => 1,
        },
    ],
};

sub help_synopsis {
    return <<EOS
gmt fastq trim-bwa-style --fastq-file=lane1.fastq --out-file=lane1.trimmed.fastq
EOS
}

sub help_detail {
    return <<EOS 
Trims fastq reads with BWA trimming style aka bwa aln -q
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

   if ( $self->qual_type eq 'sanger' ) {
        $self->_qual_str('#');
        $self->_qual_thresh(33);
    }
    else {
        $self->_qual_str('B');
        $self->_qual_thresh(64);
    }
    
    return $self;
}

sub execute {
    my $self = shift;
        
    my $reader = $self->_open_reader
        or return;
    my $writer = $self->_open_writer
        or return;

    my $out_file = $self->out_file;
    my $out_dir = dirname $out_file;
    $self->out_file($out_file);
    my $report = $self->report_file || $out_dir . '/trim.report';
    if (-e $report) {
        $self->warning_message("Reprot file: $report existing. Overwrite it");
        unlink $report;
    }
    my $report_fh = Genome::Sys->open_file_for_writing($report);
    unless ($report_fh) {
        $self->error_message("Failed to open report file " . $report . ": $!");
        return;
    }
    binmode $report_fh, ":utf8";
    
    my $ori_ct     = 0;
    my $trim_ct    = 0;
    my $rd_ori_ct  = 0;
    my $rd_trim_ct = 0;

    while ( my $seq = $reader->read ) {
        my $seq_length = length $seq->{seq};
        $ori_ct += $seq_length;
        $rd_ori_ct++;

        my $trimmed_length = $self->_trim($seq);

        if ($trimmed_length) {
            $trim_ct += $trimmed_length;
            $rd_trim_ct++;
            $report_fh->print($seq->{id}."\tT\t".$trimmed_length."\n"); #In report T for trimmed
        }

        $writer->write($seq);
    }

    my $new_ct  = $ori_ct - $trim_ct;
    my $percent = 100*$new_ct/$ori_ct;

    $report_fh->print("\nNumberOfOriginalBases  NumberOfTrimmedBases   NumberOfResultingBases  Percentage  NumberOfTrimmedReads\n");
    $report_fh->printf("%21s%22s%24s%11.1f%%%21s\n", $ori_ct, $trim_ct, $new_ct, $percent, $rd_trim_ct);
    $report_fh->close;

    return 1;
}

sub _open_reader {
    my $self = shift;

    my $fastq_file = $self->fastq_file;
    my $reader;
    eval{
        $reader = Genome::Model::Tools::Sx::FastqReader->create(
            file => $fastq_file,
        );
    };
    unless ( $reader ) {
        $self->error_message("Can't create fastq reader for file ($fastq_file): $@");
        return;
    }

    return $reader;
}

sub _open_writer {
    my $self = shift;

    my $out_file = $self->out_file;
    unless ( $out_file ) {
        my ($base_name, $base_dir) = fileparse($self->fastq_file);
        $base_name =~ s/\.fastq$//;
        $out_file = $base_dir."/$base_name.trimmed.fastq";
        $self->out_file($out_file);
        if (-e $out_file) {
            $self->warning_message("Out file: $out_file existing. Overwrite it");
            unlink $out_file;
        }
    }

    my $writer;
    eval{
        $writer = Genome::Model::Tools::Sx::FastqWriter->create(
            file => $out_file,
        );
    };
    unless ( $writer ) {
        $self->error_message("Can't create fastq writer for file ($out_file): $@");
        return;
    }

    return $writer;
}

sub trim {
    my ($self, $seqs) = @_;

    unless ( $seqs and ref $seqs eq 'ARRAY' and @$seqs ) {
        Carp::confess(
            $self->error_message("Expected array ref of sequences, but got: ".Dumper($seqs))
        );
    }

    for my $seq ( @$seqs ) {
        $self->_trim($seq);
    }

    return $seqs;
}

sub _trim {
    my ($self, $seq) = @_;

    my $seq_length = length $seq->{seq};

    my $trimmed_length;
    my ($pos, $maxPos, $area, $maxArea) = ($seq_length, $seq_length, 0, 0);

    while ($pos > 0 and $area >= 0) {
        $area += $self->trim_qual_level - (ord(substr($seq->{qual}, $pos-1, 1)) - $self->_qual_thresh);
        if ($area > $maxArea) {
            $maxArea = $area;
            $maxPos = $pos;
        }
        $pos--;
    }

    if ($pos == 0) { 
        # scanned whole read and didn't integrate to zero?  replace with "empty" read ...
        $seq->{seq}  = 'N';
        $seq->{qual} = $self->_qual_str;
        $maxPos = 1;  # reset to leave 1 base/qual as N/# there
    }
    else {  # integrated to zero?  trim before position where area reached a maximum (~where string of qualities were still below 20 ...)
        $seq->{seq}  = substr($seq->{seq},  0, $maxPos);
        $seq->{qual} = substr($seq->{qual}, 0, $maxPos);
    }

    return $seq_length - $maxPos; # trimmed length - may be 0
}

1;

#$HeadURL$
#$Id$
