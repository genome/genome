package Genome::Model::Tools::Sx::Base;

use strict;
use warnings;

use Genome;
use File::Copy::Recursive;

class Genome::Model::Tools::Sx::Base {
    is  => 'Command',
    has => [
        input => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => <<DOC
Input reader configurations. Give 'key=value' pairs, separated by a colon (:). Readers may have additonal options.

DO NOT USE this option when piping from sx commands.

Standard options:
 file => The file to read. The use of the preceding 'file=' is optional.
          It is assumed that the bare option is the file. Use '-' to read from STDIN. The file can be gzipped.
 type => The type of input. Not required if type can be determined from the file.
          Required when reading from STDIN. Valid types: sanger, bam, sam, sff, illumina, phred (fasta), fasta,
          efasta (enhanced asta).
 cnt => The number of sequences to read from the input. If the input is paired, use 2.

Additional options, by type:
 phred
  qual_file => read qualities from this file.

DOC
        }, 
        input_metrics => {
            is => 'Text',
            is_optional => 1,
            doc => 'Capture sequence metrics for the intput to this file. Current metrics include: count, bases',
        },
        output => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => <<DOC
Output writer configurations. Give 'key=value' pairs, separated by a colon (:). Writers may have additonal options.

DO NOT USE this option when piping from sx commands.

Standard options:
 file => The file to write. The use of the preceding 'file=' is optional.
          It is assumed that the bare option is the file. Use '-' to write to STDOUT. Use '.gz' to write as gzipped.
 type => The type of output. Not required if type can be determined from the file.
          Type when writing to STDOUT defaults to sanger. Valid types: sanger, illumina, phred (fasta), fasta, bed.
 mode => The mode to open the output in. Use 'w' to write to a new file (file cannot exist) and 'a' to append to a file. Default is 'w'.
 name => The name of the writer.  If using commands that attach a writer name to a sequence,
          they will be written to the specified writer.
          
          Names pair, fwd, rev and sing are reserved. They have special behavior, and don't
           require the underlying command to tag writer names to sequences.
          Examples:
          name=pair:FILE
           write only pairs
          name=pair:FILE name=sing:FILE2
           write pairs to one and singletons to another
          name=fwd:FILE name=rev:FILE2
           write first sequence to fwd, second to rev, discard singletons
          name=fwd:FILE name=rev:FILE2 name=sing:FILE3
           write first sequence to fwd, second to rev, singletons to sing
          name=sing:FILE
           write singletons to sing, discardc pairs


Additional options, by type:
 phred
  qual_file => write qualities to this file

DOC
        },
        output_metrics => {
            is => 'Text',
            is_optional => 1,
            doc => 'Capture sequence metrics for the output to this file. Current metrics include: count, bases',
        },
        _input => { is_optional => 1, },
        _output => { is_optional => 1, },
    ],
};

sub help_brief {
    return 'Transform sequences';
}

sub help_synopsis {
    return <<HELP;
    Transform sequences. See sub-commands for a additional functionality.

    Types Handled
    * sanger => fastq w/ snager quality values
    * illumina => fastq w/ illumina quality values
    * phred => fasta/quality
    * efasta => enhanced fasta w/ ambiguous bases (input only, no qual)
    * bam/sam
    * sff
    * gzipped fastq, fasta

    Things The Base Command Can Do
    * collate two inputs into one (sanger, illumina only)
    * decollate one input into two (sanger, illumina only)
    * convert type
    * remove quality fastq headers
    
HELP
}

sub help_detail {
    return <<HELP;

 / INPUT & OUPTUT FORMAT /

 * Short, where type can be determined from the file:
  gmt sx --input file.fastq --output file.fasta # qual file not written

 * Long w/ file and type:
  gmt sx --input file=file.fastq:type=illumina --output file=file.fasta:qual_file=file.qual:type:phred

 / CONVERT TYPE /

  * illumina fastq to sanger
   gmt sx --input illumina.fastq:type=illumina --output sanger.fastq

  * sanger fastq to phred fasta
   gmt sx --input file.fastq --output file.fasta

  * sanger fastq to phred fasta w/ quals
   gmt sx --input file.fastq --output file.fasta:qual_file=file.qual

 / COLLATE /

  * from paired fastqs (type-in resolved to sanger, type-out defaults to sanger)
   gmt sx --input fwd.fastq rev.fastq --output collated.fastq

  * to paired STDOUT (type-in resolved to sanger, type-out defaults to sanger)
   gmt sx --input fwd.fastq rev.fastq --output -

 / DECOLLATE /

  * from illumina to individal fastqs, discard singletons
   gmt sx --input collated.fastq --output fwd.fastq:name=fwd rev.fastq,name=rev

  * from paired illumina STDIN (type-in req'd = illumina, type-out defaults to sanger)
   gmt sx --input -:name=pair:type=illumnia --output fwd.fastq:name=fwd rev.fastq:name=rev

 / USE IN PIPE / (cmd represents an sx sub command)

  **from singleton fastq file
   gmt sx cmd1 --cmd1-options --input sanger.fastq | gmt sx cmd2 --cmd2-options --output sanger.fastq

  ** from paired STDIN to paired fastq and singleton (assuming the sub commands filter singletons)
   cat collated_fastq | gmt sx cmd1 --cmd1-options --input - --paired-input | gmt sx cmd2 --cmd2-options --output pairs.fastq

HELP
}

# FIXME rm - should have named these this way!
sub _reader { return $_[0]->_input; }
sub _writer { return $_[0]->_output; }

sub _init {
    my $self = shift;

    my $input = $self->_init_input;
    return if not $input;

    my $output = $self->_init_ouptut;
    return if not $output;

    return 1;
}

sub _init_input {
    my $self = shift;

    my @input = $self->input;
    @input = (qw/ stdinref /) if not @input;
    my %reader_params = (
        config => \@input,
    );
    if ( $self->input_metrics ) { 
        my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->create();
        return if not $metrics;
        $reader_params{metrics} = $metrics;
    }
    my $reader = Genome::Model::Tools::Sx::Reader->create(%reader_params);
    return if not $reader;
    return $self->_input($reader);
}

sub _init_ouptut {
    my $self = shift;

    my @output = $self->output;
    @output = (qw/ stdoutref /) if not @output;
    my %writer_params = (
        config => \@output,
    );
    if ( $self->output_metrics ) {
        my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->create();
        return if not $metrics;
        $writer_params{metrics} = $metrics;
    }
    my $writer = Genome::Model::Tools::Sx::Writer->create(%writer_params);
    return if not $writer;
    return $self->_output($writer);
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    $self->_add_result_observer;

    return $self;
}

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_input;
    my $writer = $self->_output;

    my $evaluator = $self->_create_evaluator;
    return if not $evaluator;

    while ( my $seqs = $reader->read ) {
        next if not $evaluator->($seqs);
        $writer->write($seqs);
    }

    return 1;
}

sub _create_evaluator { return sub{ 1;}; }

sub _add_result_observer { # to write metrics file
    my $self = shift;

    my $result_observer = $self->add_observer(
        aspect => 'result',
        callback => sub {
            #print Dumper(\@_);
            my ($self, $method_name, $prior_value, $new_value) = @_;
            # skip if new result is not successful
            if ( not $new_value ) {
                return 1;
            }

            for my $io (qw/ input output /) {
                my $method = '_'.$io;
                my $obj = $self->$method;
                return if not $obj;

                my $metrics_method = $io.'_metrics';
                my $metrics_file = $self->$metrics_method;
                next if not $metrics_file;

                my $metrics = $obj->metrics;
                if ( not $metrics ) { # very bad
                    Carp::confess("Requested to write metrics for $io, but none were found");
                }

                unlink $metrics_file if -e $metrics_file;
                my $write = $metrics->to_file($metrics_file);
                return if not $write;
            }

            return 1;
        }
    );

    if ( not defined $result_observer ) {
        Carp::confess("Cannot create result observer");
    }

    return 1;
}

1;

