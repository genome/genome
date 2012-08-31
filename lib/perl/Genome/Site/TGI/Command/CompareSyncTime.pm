package Genome::Site::TGI::Command::CompareSyncTime;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Site::TGI::Command::CompareSyncTime {
    is => 'Command::V2',
    has_input_optional => [
        log_dir => { is => 'String', default_value => '/gsc/var/log/genome/postgres/',
                     doc => 'Directory the timing log files live in' },
        date => { is => 'String', doc => 'Date to inspect in the formay yyyy-mm-dd'},
        commit => { is => 'String', doc => 'Only show information for this one commit hash' },
    ],
};

sub execute {
    my $self = shift;

$DB::single=1;
    my($pg_files, $ora_files) = $self->find_files();

    if (@$pg_files != @$ora_files) {
        $self->warning_message("The number of timing files do not match.  There are "
                               . scalar(@$pg_files) . " PostgreSQL files and "
                               . scalar(@$ora_files). " Oracle files");
    }

    @$pg_files = sort @$pg_files;
    @$ora_files = sort @$ora_files;

    if ($self->commit) {
        $self->find_matching_commit($pg_files, $ora_files);
        return 1;
    }

    $self->print_header();

    foreach my $pg_file ( @$pg_files ) {
        my $expected_ora_file = $pg_file;
        $expected_ora_file =~ s/pg/oracle/;

        unless (grep { $_ eq $expected_ora_file } @$ora_files) {
            $self->warning_message("PostgreSQL timing file $pg_file has no matching Oracle file $expected_ora_file");
            next;
        }

        my $pg_reader  = $self->_generate_reader_for_file($pg_file);
        my $ora_reader = $self->_generate_reader_for_file($expected_ora_file);

        my $next_missing = [undef, undef, undef, undef];
        # Print all the items in the Pg file
        while (my $next_pg = $pg_reader->()) {
            my $pg_commit = $next_pg->[0];
            my $next_ora = $ora_reader->($pg_commit);
            unless ($next_ora) {
                $next_ora = $next_missing;
            }
            $self->print_result_line($next_pg, $next_ora);
        }

        # The oracle reader may still have items left
        while (my $next_ora = $ora_reader->()) {
            $self->print_result_line($next_missing, $next_ora);
        }
    }
    print "\n";
}

sub _generate_reader_for_file {
    my($self, $file) = @_;

    my $fh = IO::File->new($file);
    return unless $fh;

    my @lines;
    my %commits;
    
    my $read_a_line = sub {
        my $line = <$fh>;
        chomp($line);
        return unless length($line);

        my @fields = split("\t", $line);
        my $commit = $fields[0];
        push @lines, \@fields;
        $commits{$commit} = \@fields;
    };

    my $return_next_line = sub {
        my $fields = shift @lines;
        my $commit = $fields->[0];
        delete $commits{$commit};
        return $fields;
    };

    my $return_commit = sub {
        my $commit = shift;
        my $fields = delete $commits{$commit};
        return unless $fields;

        for (my $i = 0; $i < @lines; $i++) {
            if ($lines[$i]->[0] eq $commit) {
                splice(@lines, $i, 1);
                return $fields;
            }
        }
        return;
    };

    return sub {
        my $commit = shift;
        unless (@lines) {
            # The buffer is empty
            $read_a_line->();
        }

        return $commit ? $return_commit->($commit) : $return_next_line->();
    };
}


sub find_matching_commit {
    my($self, $pg_files, $ora_files) = @_;

    my %matching_line;
    my %file_list = ( 'pg' => $pg_files, 'oracle' => $ora_files);

    my $commit = $self->commit;
    my $commit_re = qr($commit);

    SEARCH_FILE_LIST:
    foreach my $db ( 'pg','oracle' ) {
        foreach my $file ( @{ $file_list{$db} } ) {
            my $fh = IO::File->new($file);
            next unless $fh;

            while (my $line = <$fh>) {
                if ($line =~ $commit_re) {
                    my @data = split("\t",$line);
                    $matching_line{$db} = \@data;
                    next SEARCH_FILE_LIST;
                }
            }
        }
    }

    foreach my $db ( 'pg','oracle' ) {
        unless ($matching_line{$db}) {
            $self->warning_message("There was no line matching commit $commit in the $db file")
        }
    }

    $self->print_header;
    $self->print_result_line($matching_line{'pg'}, $matching_line{'oracle'});
    print "\n";
}

sub print_header {
    my $header_format = "%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%s\n";
    printf($header_format, 'PgTime', 'Oracle', 'diff', 'commit', 'exec', 'cmdline');
    printf($header_format, '------', '------', '------', '------', '------', '------');
}

sub print_result_line {
    my($self, $pg_data, $ora_data) = @_;

    # the data/file format columns are: commit-hash exec-hash time cmdline
    my @output = ( $pg_data->[2] || 'nan',
                   $ora_data->[2] || 'nan',
                   $pg_data->[2] - $ora_data->[2],
                   $pg_data->[0] || $ora_data->[0],
                   $pg_data->[1] || $ora_data->[1],
                   $pg_data->[3] || $ora_data->[3],
                );

    my $data_format = "%6.5g\t%6.5g\t%6.5g\t%6.6s\t%6.6s\t%s\n";
    printf($data_format, @output);
}


# Return pathnames for files we'll search through, given the date arg
sub find_files {
    my $self = shift;

    my($year,$month,$day) = split('-', $self->date);
    $year = '*' unless defined $year;
    $month = '*' unless defined $month;
    $day = '*' unless defined $day;

    my $base_glob_pattern = sprintf('%s/%s/%s-%s', $self->log_dir, $year, $month, $day);

    my @pg_files  = glob($base_glob_pattern . '-pg.timing');
    my @ora_files = glob($base_glob_pattern . '-oracle.timing');

    return (\@pg_files, \@ora_files);
}



1;
