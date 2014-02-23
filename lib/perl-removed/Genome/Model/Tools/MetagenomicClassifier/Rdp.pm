package Genome::Model::Tools::MetagenomicClassifier::Rdp;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MetagenomicClassifier::Rdp {
    is => 'Command',
    has => [ 
        input_file => {
            type => 'String',
            doc => "Path to fasta file"
        },
        output_file => { 
            type => 'String',
            doc => "Path to output file."
        },
        training_set => {
            type => 'String',
            valid_values => [qw/ 4 6 9 broad /],
            doc => 'Name of training set.',
        },
        version => {
            type => 'String',
            valid_values => [qw/ 2x1 2x2 2x3 2x5 /],
            doc => 'Version of rdp to run.',
        },
        format => {
            is => 'Text',
            is_optional => 1,
            valid_values => [qw/ hmp_fix_ranks hmp_all_ranks/],
            default_value => 'hmp_fix_ranks',
            doc => <<DOC,
The format of the output.
  hmp_fix_ranks => name;complemented('-' or ' ');taxon:confidence;[taxon:confidence;]
    prints only root, domain, phylum, class, order, family, genus from classification
  hmp_all_ranks => name;complemented('-' or ' ');taxon:confidence;[taxon:confidence;]
    prints ALL taxa in classification
DOC
        },
        metrics => { 
            is => 'Boolean',
            is_optional => 1, 
            doc => "Write metrics to file. File name is classification file w/ '.metrics' extension",
        },
    ],
};

sub execute {
    my $self = shift;
    
    my $reader = eval {
        Genome::Model::Tools::Sx::PhredReader->create(
            file => $self->input_file,
        );
    };
    return if not $reader;

    my $writer = Genome::Model::Tools::MetagenomicClassifier::ClassificationWriter->create(
        file => $self->output_file,
        format => $self->format,
    );
    return if not $writer;

    $self->_metrics_file; # set

    my $batch_cnt = 0;
    my $batch_sz = 50000;
    while ( 1 ) {
        my @seqs;
        while ( @seqs < $batch_sz ) {
            my $seq = $reader->read;
            last if not $seq;
            push @seqs, $seq;
        }
        last if not @seqs;
        my $batch_pos = $batch_cnt * $batch_sz;
        $batch_cnt++;
        $self->debug_message('Batch '.$batch_cnt.' from '.($batch_pos + 1).' to '.($batch_pos + @seqs));
        my $pid = UR::Context::Process->fork();
        if ( not defined $pid ) {
            die 'Cannot fork!';
        }
        elsif ( $pid == 0 ) {
            #child
            my $classifier_class = 'Genome::Model::Tools::MetagenomicClassifier::Rdp::Version'.$self->version;
            my $classifier = $classifier_class->create(
                training_set => $self->training_set,
            );
            die $self->error_message('Failed to create classifier') if not $classifier;
            my $metrics = $self->_load_metrics;
            die $self->error_message('Failed to load metrics') if not $metrics;
            for my $seq ( @seqs ) {
                #print $seq->{id}."\n";
                $metrics->{total}++;
                my $classification = $classifier->classify($seq);
                if ( not $classification ) {
                    $metrics->{error}++;
                    next;
                }
                $writer->write($classification); # TODO check?
                $metrics->{success}++;
            }
            my $write_metrics = $self->_write_metrics($metrics);
            die $self->error_message('Failed to write metrics') if not $write_metrics;
            exit 0;
        }
        else {
            # parent
            waitpid($pid, 0);
            my $exit_code = $? >> 8;
            if(  $exit_code != 0 ) {
                $self->error_message('Error in running child RDP process.');
                return;
            }
        }
    }

    return 1;
}

sub _metrics_file {
    my $self = shift;

    if ( not $self->{_metrics_file} ) {
        if ( $self->metrics ) {
            $self->{_metrics_file} = $self->output_file.'.metrics';
            unlink $self->{_metrics_file};
        }
        else { # use tmp file that is cleaned up
            my $tmpdir = File::Temp::tempdir(CLEANUP => 1); 
            $self->{_metrics_file} = $tmpdir.'/metrics';
        }
    }

    return $self->{_metrics_file};
}

sub _load_metrics {
    my $self = shift;

    my $metrics_file = $self->_metrics_file;
    my %metrics;
    @metrics{qw/ total success error /} = (qw/ 0 0 0 /);
    return \%metrics if not -e $metrics_file;

    my $metrics_fh = eval{ Genome::Sys->open_file_for_reading($metrics_file); };
    if ( not $metrics_fh ) {
        $self->error_message("Failed to open metrics file ($metrics_file) for reading: $@");
        return;
    }

    while ( my $line = $metrics_fh->getline ) {
        chomp $line;
        my ($key, $value) = split('=', $line);
        $metrics{$key} = $value;
    }

    return \%metrics;
}


sub _write_metrics {
    my ($self, $metrics) = @_;

    my $metrics_file = $self->_metrics_file;
    unlink $metrics_file if -e $metrics_file;

    my $metrics_fh = eval{ Genome::Sys->open_file_for_writing($metrics_file); };
    if ( not $metrics_fh ) {
        $self->error_message("Failed to open metrics file ($metrics_file) for writing: $@");
        return;
    }

    for my $metric ( keys %$metrics ) {
        $metrics_fh->print($metric.'='.$metrics->{$metric}."\n");
    }

    $metrics_fh->close;

    return 1;
}

#< HELP >#
sub help_brief {
    "Classify sequences with rdp",
}

sub help_detail {
    return <<HELP;
   This tool will take a fasta file and output RDP classifications. An attempt will be made to classify each sequence. If it cannot be classified, an error will be displayed, and the classifier will contiue on to the next sequence.

HELP
}

1;

