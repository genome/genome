use strict;
use warnings;

package Genome::Utility::Test;
use base 'Test::Builder::Module';

use Exporter 'import';

our @EXPORT_OK = qw(compare_ok run_ok capture_ok abort strip_ansi
    command_execute_ok command_execute_fail_ok is_equal_set);

use Carp qw(croak);
use IPC::System::Simple qw(capture);
use Sub::Install;
use Test::More;
use File::Spec qw();
use Set::Scalar;


my %ERRORS = (
    'REPLACE_ARRAY_REF' => q('replace' value should be an ARRAY ref),
);
sub ERRORS {
    my ($class, $key) = @_;
    return $ERRORS{$key};
}

sub abort {
    diag "  Aborted.";
    die;
}

sub _compare_ok_parse_args {
    # First two args are always the files.
    # Then accept one scalar argument, the name.
    # And one HASH ref, the options.
    # Then validate the option names.
    my $file_1 = shift;
    my $file_2 = shift;

    my $name = (@_ % 2 == 1) ? shift : undef;
    my %o = @_;

    if ($name && $o{name}) {
        die 'duplicate name argument not expected';
    } elsif ($name && !$o{name}) {
        $o{name} = $name;
    }

    my %vo; # validated_o
    $vo{name}    = delete $o{name};
    my $filters = delete($o{filters}) || [];
    my $replace = delete($o{replace}) || [];
    $vo{xform} = [];
    $vo{diag}    = delete $o{diag} // 1;
    my @k = keys %o;
    if (@k) {
        croak 'unexpected options passed to compare_ok: ' . join(', ', @k);
    }

    my $replace_ref = ref($replace);
    if (defined $replace && (!$replace_ref || $replace_ref ne 'ARRAY')) {
        die sprintf(q(%s: %s\n), $ERRORS{REPLACE_ARRAY_REF}, $replace);

    } else {
        # inspect each thing in replace
        for(my $i = 0; $i < @$replace; $i++) {
            my $r = $replace->[$i];

            my $coderef;
            if (ref($r) eq 'CODE') {
                # coderef - use it directly
                $coderef = $r;
            } elsif (ref($r) eq 'ARRAY' and (@$r >= 1) and (@$r <= 2)) {
                # list of one or two somethings

                my($regex, $replacement);
                if (ref($r->[0]) eq 'Regexp') {
                    # refered form
                    $regex = $r->[0];
                } elsif (! ref($r->[0])) {
                    # string - convert to regex
                    $regex = $r->[0];
                    $regex = qr($regex);
                } else {
                     Carp::croak("Unexpected options passed to compare_ok's 'replace' item $i. "
                                    . "Expected the first element to be a string or Regex, "
                                    . "but got ".ref($r->[0]));
                }
                # now convert to a coderef
                $replacement = defined($r->[1]) ? $r->[1] : '';  # Default replacement is empty string
                $coderef = sub {
                    my $orig = shift;
                    (my $changed = $orig) =~ s/$regex/$replacement/g;
                    return $changed;
                };

            } else {
                Carp::croak("Unexpected options passed to compare_ok's 'replace' item $i. "
                            . "Expected either a coderef or an arrayref of 1 or 2 items, "
                            . "but got ".ref($r));
            }
            push @{$vo{xform}}, $coderef;
        }
    }

    my $filters_ref = ref($filters);
    if (defined $filters && (!$filters_ref || $filters_ref ne 'ARRAY')) {
        # it's a simple string
        $filters = [ $filters ];
    }
    for (my $i = 0; $i < @$filters; $i++) {
        my $filter = $filters->[$i];

        my $coderef;
        if (ref($filter) eq 'CODE') {
            $coderef = $filter;
        } elsif (!ref($filter) or ref($filter) eq 'Regexp') {
            # simple string or regex
            $coderef = sub {
                my $orig = shift;
                (my $changed = $orig) =~ s/$filter//;
                return $changed;
            };
        } else {
            Carp::croak("Unexpected options passed to compare_ok's 'filter' item $i. "
                        . "Expected a coderef, Regexp or string, but got ".ref($filter));
        }
        push @{$vo{xform}}, $coderef;
    }


    return ($file_1, $file_2, %vo);
}

sub _compare_ok_iterator_for_file {
    my $file = shift;

    my $fh = IO::File->new($file, 'r') || die "Can't open $file for reading: $!";
    return sub {
        my $arg = shift;

        if ($arg eq 'input_line_number') {
            # Hack to get the line number
            return $fh->input_line_number();
        }

        # $arg must be a listref of transform functions
        NEXT_LINE:
        while(my $line = $fh->getline()) {
            foreach my $xform ( @$arg ) {
                $line = $xform->($line);
            }
            return $line if(defined($line) and length($line));   # go again if the xforms eliminated this line
        }
        return;  # must be EOF
    };
}

sub compare_ok {
    my ($file_1, $file_2, %o) = _compare_ok_parse_args(@_);

    my $tb = __PACKAGE__->builder;

    $tb->no_diag(!$o{diag});

    my(@iters,@filename);
    foreach my $file ( $file_1, $file_2 ) {
        push @filename, $file;
        push @iters, _compare_ok_iterator_for_file($file);
    }

    my @xforms;
    if ($o{xform}) {
        # replace means apply this xform to all input files
        # FIXME add options to apply a transform to either file1 or file2
        for (my $i = 0; $i < @iters; $i++) {
            $xforms[$i] ||= [];
            push @{ $xforms[$i] }, @{$o{xform}};
        }
    }

    my @lines = (1); # a dummy value to satisfy the grep below the first time through the loop
    my $result = 1;
    GET_LINE_FROM_FILES:
    while (grep { defined($_) && length($_) } @lines) {

        # fill in the next line for each
        for (my $i = 0; $i < @iters; $i++) {
            $lines[$i] = $iters[$i]->( $xforms[$i] );
        }

        # now compare them
        COMPARISON:
        for (my $i = 1; $i < @lines; $i++) {
            {   no warnings 'uninitialized';
                next COMPARISON if $lines[0] eq $lines[$i];
            }

            # different
            my($line1,$line2) = @lines[0,$i];
            chomp($line1, $line2);

            $tb->diag(sprintf("First diff:\n--- %s line %d\n+++ %s line %d\n-%s\n+%s\n",
                    $filename[0], $iters[0]->('input_line_number'),
                    $filename[$i], $iters[$i]->('input_line_number'),
                    $line1, $line2));

            $result = 0;
            last GET_LINE_FROM_FILES;
        }
    }

    return $tb->ok($result, $o{name});
}

sub capture_ok {
    my ($command, $test_name) = @_;

    my $tb = __PACKAGE__->builder;

    my @command = ref $command ? @$command : $command;
    $test_name //= @command > 1 ? join(' ', @command) : $command[0];

    my @output = eval { capture(@command) };
    my $error = $@;
    my $exit_zero = ($? == 0);
    $tb->ok($exit_zero, $test_name) or diag $error, @output;

    if (wantarray) {
        return ($exit_zero, @output);
    } else {
        return $exit_zero;
    }
}

sub run_ok {
    my ($command, $test_name) = @_;

    my $tb = __PACKAGE__->builder;

    my @command = ref $command ? @$command : $command;
    $test_name //= @command > 1 ? join(' ', @command) : $command[0];

    my $exit_zero = (system(@command) == 0);
    $tb->ok($exit_zero, $test_name);

    return $exit_zero;
}

sub is_equal_set {
    my ($ol, $el, $name) = @_;
    my $tb = __PACKAGE__->builder;
    my $os = Set::Scalar->new(@$ol);
    my $es = Set::Scalar->new(@$el);
    my $ok = $tb->ok($os->is_equal($es), $name);
    unless ($ok) {
        my $extra = $os - $es;
        unless ($extra->is_empty) {
            diag 'Extra: ' . join(', ', $extra->members);
        }

        my $missing = $es - $os;
        unless ($missing->is_empty) {
            diag 'Missing: ' . join(', ', $missing->members);
        }
    }
}

sub data_dir_ok {
    my $data_dir = data_dir(@_);
    my $tb = __PACKAGE__->builder;
    $tb->ok(-d $data_dir, "data_dir exists: $data_dir");
    return $data_dir;
}

sub data_dir {
    my ($class, $package, $test_version) = @_;

    # "validate" package
    $package->class;

    (my $dirname = $package) =~ s/::/-/g;
    my @parts = ($ENV{GENOME_TEST_INPUTS}, $dirname);
    if ($test_version) {
        push @parts, $test_version;
    }
    my $dirpath = File::Spec->join(@parts);
    return $dirpath;
}

sub strip_ansi {
    my $string = shift;
    $string =~ s/\e\[\d+(?>(;\d+)*)m//g;
    return $string;
}

my @command_message_types = qw(status warning error debug usage);
sub _command_execute_ok_parse_args {
    if (@_ < 1 or @_ > 3) {
        Carp::croak("Expected 1, 2 or 3 args to command_execute_ok(), but got ".scalar(@_));
    }

    my $command = shift;
    unless (ref($command) && $command->isa('Command')) {
        Carp::croak('Arg 1 to command_execute_ok() must be an instance of a Command object');
    }

    my($message_config, $ok_msg);
    if (! @_) {
        # no more args
        $message_config = {};
        $ok_msg = "execute ".ref($command);

    } elsif (@_ == 1) {
        $ok_msg = shift;
        $message_config = {};

    } else { # @_ == 2
        $message_config = shift;
        $ok_msg = shift;
    }

    # validate $message_config contents
    foreach my $type ( @command_message_types ) {
        # ! exists means don't filter these messages at all
        next unless exists $message_config->{$type.'_messages'};
        # ! defined means don't print them to the terminal, and don't check them afterward
        next unless defined $message_config->{$type.'_messages'};

        if (ref($message_config->{$type.'_messages'}) ne 'ARRAY') {
            Carp::croak("Expected an arrayref for ${type}_messages, but got "
                            . ref($message_config->{$type.'_messages'}));
        }
        foreach my $check_msg ( @{ $message_config->{$type.'_messages'} }) {
            if (ref($check_msg) && ref($check_msg) ne 'Regexp') {
                Carp::croak("Values for ${type}_messages must be strings or Regexps");
            }
        }
    }
    return ($command, $message_config, $ok_msg);
}

my %format_numeral = qw( 1 st 2 nd 3 rd 4 th 5 th 6 th 7 th 8 th 9 th 0 th );

_command_execute_ok_builder('command_execute_ok', 1);
_command_execute_ok_builder('command_execute_fail_ok', 0);
sub _command_execute_ok_builder {
    my($subname, $should_work) = @_;

    my $sub = sub {
        my($command, $message_config, $message) = _command_execute_ok_parse_args(@_);

        foreach my $type (@command_message_types) {
            if (exists $message_config->{$type.'_messages'}) {
                my($dump, $queue) = map { "${_}_${type}_messages" } qw(dump queue);
                $command->$dump(0);
                $command->$queue(1) if (defined $message_config->{$type.'_messages'});
            }
        }

        my $tb = __PACKAGE__->builder;
        my $result = $command->execute();
        { no warnings 'uninitialized';
            # result might be undef
            ($result ^ $should_work) && return $tb->ok(0, "$message (execute() returned false)");
        }

        foreach my $type (@command_message_types) {
            my $method = $type.'_messages';
            if (defined $message_config->{$method}) {
                my @expected_messages = @{ $message_config->{$method} };
                my @got_messages = $command->$method();
                for (my $i = 0; @expected_messages || @got_messages; $i++) {
                    my($expected) = map { defined $_ ? $_ : '' } shift @expected_messages;
                    my($got) = map { defined $_ ? $_ : '' } shift @got_messages;

                    if (ref($expected) and $got !~ m/$expected/) {
                        my $rv = $tb->ok(0, $message);
                        $i++;
                        $tb->diag("For the $i" .$format_numeral{substr($i, -1)}
                                    . ' ' . substr($method, 0, -1) # remove the 's'
                                    . ", '$got' didn't match $expected");
                        return $rv;

                    } elsif (!ref($expected) and $got ne $expected) {
                        my $rv = $tb->ok(0, $message);
                        $i++;
                        $tb->diag("For the $i" .$format_numeral{substr($i, -1)}
                                    . ' ' . substr($method, 0, -1) # remove the 's'
                                    . ", got '$got' but expected '$expected'");
                        return $rv;
                    }
                }
            }
        }
        return $tb->ok(1, $message);
    };
    Sub::Install::install_sub({
        as => $subname,
        code => $sub,
    });
}


1;

__END__

=pod

=head1 NAME

Genome::Utility::Test

=head1 SYNOPSIS

    use Genome::Utiltiy::Test qw(compare_ok);

    compare_ok($file_1, $file_1);

=head1 METHODS

=over

=item strip_ansi($string)

Returns the given string after removing ANSI escape sequences.

=item compare_ok($file1, $file2, %options)

Compare two files line-by-line

  compare_ok($file_1, $file_2)

With no options, it directly compares two files.  At the first difference, it
will stop and print a diagnostic message about the difference (if diag => 1
is passed as an option).

Options include

=over

=item diag => 1 || 0

Do or do not print a diag() message if the files are different

=item filters => [ string1, string2, ..., regex1, regex2, ...]

A list of string or regular expressions to remove from both files before
comparing them.  Useful to strip out timestamps or usernames that may be
different in normal operation.

=item replace => [ [ string1 => replace_string1 ], [regex1 => replace_string2] ]

Like filter, but allows the matching string or regex to be replaced with
another string before the line comparison is done.

=item name

Specify a name for a test but requires 'name' key.

=item diag

Disable diag output when a diff is encountered. Added this in case people want to surpress any output.

=item test

Disable test usage, just return status. Added this so I could test compare_ok.

=back

=item command_execute_ok

  command_execute_ok($cmd)
  command_execute_ok($cmd, $test_name);
  command_execute_ok($cmd, $expected_messages, $test_name)

Executes an instance of a Command class, and optionally checks that it
produced the expected error_messages, status_messages, etc.  $expected_messages
is a hashref keyed by the message types (error_messages, status_messaegs,
debug_messages, warning_messages and usage_messages) with values that are
arrayrefs of strings or regexes, for example

  command_execute_ok($cmd,
                  { error_messages => ['error 1', qr(broken) ],
                    status_messages => ['in progress'],
                    warning_messages => [],
                    debug_messages => undef },
                  'Try command');

An empty list means the test expects the command to generate none of that type
of message.  undef means that type of message will not be printed to the
terminal during the execution, but the message contents will not be checked.

=item command_execute_fail_ok

Like command_execute_ok(), but this test fails if the command's execution returns true.

=back

=head1 SEE ALSO

Test::More

