package Genome::Model::Tools::Sx::Writer;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Writer {
    has => [
        config => { is => 'Text', is_many => 1, },
        metrics => { is => 'Genome::Model::Tools::Sx::Metrics', is_optional => 1, },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my @config = grep { defined } $self->config;
    if ( not @config ) { 
        $self->error_message('No config given to write');
        return;
    }

    if ( grep { $_ =~ /stdoutref/ } @config ) {
        if ( @config > 1 ) {
            $self->error_message('Cannot write stdout refs and to other writers');
            return;
        }
        my $writer = Genome::Model::Tools::Sx::StdoutRefWriter->create;
        return if not $writer;
        $self->{_writer} = $writer;
        $self->{_strategy} = 'write_stdoutref';
        return $self;
    }

    my @writers;
    for my $config ( $self->config ) {
        $self->debug_message('Parsing writer config: '.$config);
        my ($writer_class, $params) = $self->parse_writer_config($config);
        return if not $writer_class;

        $self->debug_message('Config: ');
        $self->debug_message('writer => '.$writer_class);
        for my $key ( keys %$params ) {
            $self->debug_message($key.' => '.$params->{$key});
        }

        my $writer = $writer_class->create(%$params);
        if ( not $writer ) {
            $self->error_message('Failed to create '.$writer_class);
            return;
        }
        push @writers, $writer;
    }

    my $strategy = $self->_resolve_strategy(@writers);
    return if not $strategy;
    $self->debug_message('Write strategy: '.$self->{_strategy});

    return $self;
}

sub parse_writer_config {
    my ($self, $config) = @_;

    Carp::confess('No config to parse') if not $config;

    my %params = Genome::Model::Tools::Sx::Functions->config_to_hash($config);
    if ( not %params ) {
        $self->error_message("Failed to parse config! $config");
        return;
    }

    if ( not $params{file} ) {
        $self->error_message('Failed to get "file" from config: '.$config);
        return;
    }

    my $type = delete $params{type};
    if ( not $type ) {
        $type = $self->_type_for_file($params{file});
        return if not $type;
    }

    my $writer_class = $self->_writer_class_for_type($type);
    return if not $writer_class;

    return ($writer_class, \%params);
}

sub _type_for_file {
    my ($self, $file) = @_;

    Carp::confess('No file to get type') if not $file;

    if ( $file eq '-' ) {
        return 'sanger';
    }

    $file =~ s/\.gz$//;
    my ($ext) = $file =~ /\.(\w+)$/;
    if ( not $ext ) {
        $self->error_message('Failed to get extension for file: '.$file);
        return;
    }

    my %exts_and_types = (
        bed => 'bed',
        fasta => 'phred',
        fna => 'phred',
        fa => 'phred',
        fastq => 'sanger',
        fq => 'sanger',
    );
    if ( $exts_and_types{$ext} ) {
        return $exts_and_types{$ext}
    }

    $self->error_message('Failed to determine type for file: '.$file);
    return;
}

sub _writer_class_for_type {
    my ($self, $type) = @_;

    Carp::confess('No type to get writer class') if not $type;

    my %types_and_classes = (
        bed => 'BedWriter',
        phred => 'PhredWriter',
        sanger => 'FastqWriter',
        illumina => 'IlluminaFastqWriter',
    );

    if ( exists $types_and_classes{$type} ) {
        my $writer_class = 'Genome::Model::Tools::Sx::'.$types_and_classes{$type};
        return $writer_class;
    }

    $self->error_message('Invalid type: '.$type);
    return;
}

sub _resolve_strategy {
    my ($self, @writers) = @_;

    Carp::confess('No writer to resolve strategy!') if not @writers;

    my @names = grep { defined } map { $_->name } @writers;
    # no names - default
    if ( not @names ) { # Write allseqs to each writer
        $self->{_writers} = \@writers;
        $self->{_strategy} = 'write_to_all';
        return 1;
    }

    if ( @names != @writers ) {
        $self->error_message('Can not mix writers with names and ones without');
        return;
    }

    # pair fwd rev sing do not have to have writer names in the seqs
    my $pfrs = $self->_resolve_strategy_for_pair_fwd_rev_and_singleton(@writers); # -1 undef 1
    return if not $pfrs;
    return 1 if $pfrs == 1;

    # writers are named, write only seqs to them w/ that name - seqs must have that writer name
    $self->{_strategy} = 'write_to_named_writers';
    for my $writer ( @writers ) {
        $self->{ $writer->name } = $writer
    }

    return 1;
}

sub _resolve_strategy_for_pair_fwd_rev_and_singleton {
    my ($self, @writers) = @_;

    my @pair = grep { $_->name eq 'pair' } @writers;
    my @fwd = grep { $_->name eq 'fwd' } @writers;
    my @rev = grep { $_->name eq 'rev' } @writers;
    my @sing = grep { $_->name eq 'sing' } @writers;

    if ( not ( @pair or @fwd or @rev or @sing ) ) {
        return -1 
    }

    if ( ( @pair or @fwd or @rev or @sing ) and @pair + @fwd + @rev + @sing != @writers ) {
        $self->error_message('Can not mix pair/singleton writers and other named writers');
        return;
    }
    if ( @pair and ( @fwd or @rev ) ) {
        $self->error_message('Can not mix pair/fwd/rev writers');
        return;
    }
    if ( ( @fwd and not @rev ) or ( @rev and not @fwd ) ) {
        $self->error_message('Have one fwd/rev writer, but missing the other');
        return;
    }

    if ( @pair ) {
        if ( @sing ) {
            $self->{_writers} = \@pair;
            $self->{_sing} = \@sing;
            $self->{_strategy} = 'write_pairs_and_singletons_separately';
        }
        else {
            $self->{_writers} = \@pair;
            $self->{_strategy} = 'write_to_all';
        }
        return 1;
    }
    elsif ( @fwd and @rev ) {
        $self->{_fwd} = \@fwd;
        $self->{_rev} = \@rev;
        if ( @sing ) {
            $self->{_sing} = \@sing;
            $self->{_strategy} = 'write_forward_reverse_and_singletons_separately';
        }
        else {
            $self->{_strategy} = 'write_forward_and_reverse_separately';
        }
        return 1;
    }
    elsif ( @sing ) {
        $self->{_sing} = \@sing;
        $self->{_strategy} = 'write_singletons_only';
        return 1;
    }

    return -1;
}

sub write {
    my ($self, $seqs) = @_;
    Carp::confess('No sequences to write!') if not $seqs or not @$seqs;
    my $strategy = $self->{_strategy};
    return $self->$strategy($seqs);
}

sub write_stdoutref {
    $_[0]->metrics->add_sequences($_[1]) if $_[0]->metrics;
    return $_[0]->{_writer}->write($_[1]);
}

sub write_to_all {
    my ($self, $seqs) = @_;

    for my $writer ( @{$self->{_writers}} ) {
        for my $seq ( @$seqs ) {
            $writer->write($seq) or return;
        }
    }

    $self->metrics->add_sequences($seqs) if $self->metrics;

    return 1;
}

sub write_pairs_and_singletons_separately {
    my ($self, $seqs) = @_;

    if ( @$seqs == 1 ) {
        return $self->write_singletons_only($seqs);
    }
    elsif ( @$seqs == 2 ) {
        return $self->write_to_all($seqs);
    }

    $self->error_message('Too many sequences to write: '.scalar(@$seqs));
    return;
}

sub write_singletons_only {
    my ($self, $seqs) = @_;

    return 1 if @$seqs != 1;

    for my $writer ( @{$self->{_sing}} ) {
        $writer->write($seqs->[0]) or return;
    }

    $self->metrics->add_sequences($seqs) if $self->metrics;

    return 1;
}

sub write_forward_reverse_and_singletons_separately {
    my ($self, $seqs) = @_;

    if ( @$seqs == 1 ) {
        return $self->write_singletons_only($seqs);
    }
    elsif ( @$seqs == 2 ) {
        return $self->write_forward_and_reverse_separately($seqs);
    }

    $self->error_message('Too many sequences to write: '.scalar(@$seqs));
    return;
}

sub write_forward_and_reverse_separately {
    my ($self, $seqs) = @_;

    return 1 if @$seqs != 2;

    for my $writer ( @{$self->{_fwd}} ) {
        $writer->write($seqs->[0]) or return;
    }

    for my $writer ( @{$self->{_rev}} ) {
        $writer->write($seqs->[1]) or return;
    }

    $self->metrics->add_sequences($seqs) if $self->metrics;

    return 1;
}

sub write_to_named_writers {
    my ($self, $seqs) = @_;

    for my $seq ( @$seqs ) {
        if ( not $seq->{writer_name} ) { # OK, write to discard or nothing
            if ( $self->{discard}  ) {
                $self->{discard}->write($seq);
            }
            next;
        }
        if ( not $self->{ $seq->{writer_name} } ) { # is there a writer with this name?
            $self->error_message('Attempting to write sequences to named writers, but there is not a writer for name: '.$seq->{writer_name});
            return;
        }
        $self->{ $seq->{writer_name} }->write($seq) or return;
        $self->metrics->add_sequence($seq) if $self->metrics;
    }

    return 1;
}

1;

