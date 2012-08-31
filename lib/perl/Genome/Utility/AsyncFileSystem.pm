package Genome::Utility::AsyncFileSystem;

use strict;
use warnings;

use Genome;
use AnyEvent::Util;
use Exporter qw( import );

our @EXPORT_OK = qw(on_each_line);

class Genome::Utility::AsyncFileSystem {};

# execute a shell command in a standard way instead of using system()\
# verifies inputs and ouputs, and does detailed logging...
# this version returns a condition variable, to wait for its results
# you need to call ->recv on it which will return 1 just like GUFS->shellcmd.
sub shellcmd {
    my ( $self, %params ) = @_;
    my $cmd                    = delete $params{cmd};
    my $output_files           = delete $params{output_files};
    my $input_files            = delete $params{input_files};
    my $output_directories     = delete $params{output_directories};
    my $input_directories      = delete $params{input_directories};
    my $allow_failed_exit_code = delete $params{allow_failed_exit_code};
    my $allow_zero_size_output_files =
      delete $params{allow_zero_size_output_files};
    my $skip_if_output_is_present = delete $params{skip_if_output_is_present};
    $skip_if_output_is_present = 1 if not defined $skip_if_output_is_present;

    my %ae_params = ();
    my @ae_param_list =
      ( '<', '>', 'on_prepare', 'close_all', 'close_all_fds_except', '$$' );

    for (@ae_param_list) {
        $ae_params{$_} = delete $params{$_} if exists $params{$_};
    }

    for my $p ( keys %params ) {
        if ( $p =~ /\d>/ || $p =~ /\d</ ) {
            $ae_params{$p} = delete $params{$p};
        }
    }

    if (%params) {
        my @crap = %params;
        Carp::confess("Unknown params passed to shellcmd: @crap");
    }

    if ( $output_files and @$output_files ) {
        my @found_outputs = grep { -e $_ } grep { not -p $_ } @$output_files;
        if (    $skip_if_output_is_present
            and @$output_files == @found_outputs )
        {
            $self->status_message( "SKIP RUN (output is present):     $cmd\n\t"
                  . join( "\n\t", @found_outputs ) );
            my $cv = AnyEvent->condvar;
            $cv->send(1);
            return $cv;
        }
    }

    if ($input_files and @$input_files) {
        my @missing_inputs = grep { not -s $_ } grep { not -p $_ } @$input_files;
        if (@missing_inputs) {
            Carp::croak("CANNOT RUN (missing input files):     $cmd\n\t"
                         . join("\n\t", map { -e $_ ? "(empty) $_" : $_ } @missing_inputs));
        }
    }

    if ($input_directories and @$input_directories) {
        my @missing_inputs = grep { not -d $_ } @$input_directories;
        if (@missing_inputs) {
            Carp::croak("CANNOT RUN (missing input directories):     $cmd\n\t"
                        . join("\n\t", @missing_inputs));
        }
    }

    $self->status_message("RUN: $cmd");

    my $exit_code;
    my $wrapper_cv = AnyEvent->condvar;
    $wrapper_cv->begin(
        sub {
            my ($cv) = shift;

            eval {
                if ( $exit_code == -1 )
                {
                    Carp::croak("ERROR RUNNING COMMAND. Failed to execute: $cmd");
                } elsif ( $exit_code & 127 ) {
                    my $signal = $exit_code & 127;
                    my $withcore = ( $exit_code & 128 ) ? 'with' : 'without';
                    Carp::croak("COMMAND KILLED. Signal $signal, $withcore coredump: $cmd");
                } elsif ( $exit_code >> 8 != 0 ) {
                    $exit_code = $exit_code >> 8;
                    $DB::single = $DB::stopper;
                    if ($allow_failed_exit_code) {
                        Carp::carp("TOLERATING Exit code $exit_code, msg $! from: $cmd");
                    } else {
                        Carp::croak("ERROR RUNNING COMMAND.  Exit code $exit_code, msg $! from: $cmd");
                    }
                }

                my @missing_output_files;
                if ( $output_files and @$output_files ) {
                    @missing_output_files =
                      grep { not -s $_ } grep { not -p $_ } @$output_files;
                }
                if (@missing_output_files) {
                    if (    $allow_zero_size_output_files
                        and @$output_files == @missing_output_files )
                    {
                        for my $output_file (@$output_files) {
                            Carp::carp("ALLOWING zero size output file '$output_file' for command: $cmd");
                            my $fh = $self->open_file_for_writing($output_file);
                            unless ($fh) {
                                Carp::croak("failed to open $output_file for writing to replace missing output file: $!");
                            }
                            $fh->close;
                        }
                        @missing_output_files = ();
                    }
                }

                my @missing_output_directories;
                if ( $output_directories and @$output_directories ) {
                    @missing_output_directories =
                      grep { not -s $_ }
                      grep { not -p $_ } @$output_directories;
                }

                if ( @missing_output_files or @missing_output_directories ) {
                    for (@$output_files) { unlink $_ or Carp::croak("Can't unlink $_: $!"); }
                    Carp::croak("MISSING OUTPUTS! "
                                . join(', ', @missing_output_files)
                                . " "
                                . join(', ', @missing_output_directories));
                }
            };
            if ($@) {
                $cv->croak($@);
            } else {
                $cv->send(1);
            }
        }
    );

    my $acv = $self->_async_cmd( $cmd, \%ae_params );

    $acv->cb(
        sub {
            $exit_code = shift->recv;
            $wrapper_cv->end;
        }
    );

    return $wrapper_cv;
}

# this function returns an AnyEvent::CondVar
# to block until the cmd returns call ->recv
# on the condvar.  recv returns the exit code
# in the same form as $?
sub _async_cmd {
    my ( $self, $cmd, $options ) = @_;

    $options ||= {};
    $options->{close_all} = 1
      unless exists $options->{close_all}
          || exists $options->{close_all_fds_except};

    $options->{'<'} = \*STDIN
      unless exists $options->{'<'};
    $options->{'>'} = \*STDOUT
      unless exists $options->{'>'};
    $options->{'2>'} = \*STDERR
      unless exists $options->{'2>'};

    return run_cmd $cmd, %$options;
}

## Utility functions for writing i/o callbacks
#  These should be in EXPORT_OK
sub on_each_line (&) {
    my $code = shift;

    my $buf = '';
    return sub {
        my $data = shift;

        $buf .= $data if (defined $data);
        while (1) {
            my $pos = index( $buf, "\n" );
            last if ( $pos < 0 );

            my $line = substr( $buf, 0, $pos + 1, '' );
            $code->($line);
        }

        unless (defined $data) {
            $code->($buf) if ($buf ne '');
            undef $buf;
            $code->();
        }
      }
}

1;
