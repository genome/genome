package Genome::SoftwareResult::Command::DiskUsage;

use strict;
use warnings;

use feature qw(say);

use Genome;

class Genome::SoftwareResult::Command::DiskUsage {
    is => 'Command::V2',
    has => [
        format => {
            is => 'Text',
            valid_values => ['tsv', 'fixed'],
            default_value => 'fixed',
            doc => 'Delimit the report columns with tabs, or output them in a fixed-width format',
        },
        active_only => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Whether to only report "active" allocations (notably excluding "archived" allocations)',
        },
    ],
    doc => 'report on the disk usage of software results',
};

sub help_detail {
    return <<EOHELP
This command produces a report aggregating the disk usage of SoftwareResults by "sponsor".
The output consists of three columns:
 - the sponsor
 - the KB used by allocations owned by that sponsor
 - the TB used by allocations owned by that sponsor (equal to KB / (1024^3))

This command will take several minutes to run.
EOHELP
}

sub execute {
    my $self = shift;

    $self->allocation_chart;

    return 1;
}

sub _dbh {
    my $self = shift;

    unless (defined $self->{_dbh}) {
        $self->{_dbh} = Genome::Config::AnalysisProject->__meta__->data_source->get_default_handle;
    }

    return $self->{_dbh};
}

sub _table_row_function {
    my $self = shift;

    my $format = $self->format;
    if ($format eq 'fixed') {
        return sub {say join("\t", sprintf('%-100s', shift()), map { sprintf('%-9s', $_) } @_) };
    } elsif ($format eq 'tsv') {
        return sub {say join("\t", @_) };
    } else {
        #one would hope there's a case here for each valid value of the property!
        $self->fatal_message('We got an unexpected format: %s', $format);
    }
}

sub allocation_summary {
    my $self = shift;
    my $dbh = $self->_dbh();

    my @disk_group_names = (Genome::Config::get('disk_group_alignments'), Genome::Config::get('disk_group_models'));

    my $active_clause = ($self->active_only? q(AND status = 'active') : '');

    my $sth = $dbh->prepare(q(SELECT sru.user_id,SUM(da.kilobytes_requested / counts.num_users) usage
FROM result."user" sru
INNER JOIN result.software_result sr ON sr.id = sru.software_result_id
INNER JOIN (
    SELECT sr.id software_result_id, COUNT(sru.user_id) num_users FROM result.software_result sr
    INNER JOIN result."user" sru ON sr.id = sru.software_result_id
    WHERE sru.label = 'sponsor'
    GROUP BY sr.id
) counts ON sr.id = counts.software_result_id
INNER JOIN disk.allocation da ON da.owner_id = sr.id
WHERE sru.label = 'sponsor' AND da.disk_group_name IN (?, ?) ) . $active_clause . q(
GROUP BY sru.user_id
ORDER BY usage DESC;
))
        or die $dbh->errstr;

    $sth->execute(@disk_group_names) or die $sth->errstr;
    return $sth->fetchall_arrayref;
}

sub display_allocation_info {
    my $self = shift;

    my $allocation_data = shift;

    my $row_printer = $self->_table_row_function();

    $row_printer->('Owner', 'KB', 'TB');
    $row_printer->(('--') x 3);
    for my $data (@$allocation_data) {
        $row_printer->(
            $data->[0],
            sprintf('%d', $data->[1]),
            sprintf('%.3f', &_kb_to_tb($data->[1]))
        );
    }
}

sub _kb_to_tb {
    return shift() / 1024 / 1024 / 1024;
}

sub allocation_chart {
    my $self = shift;
    my $allocation_data = $self->allocation_summary;
    for my $data (@$allocation_data) {
        my $anp = Genome::Config::AnalysisProject->get($data->[0]);
        $data->[0] = $anp->name if $anp;
    }

    $self->display_allocation_info($allocation_data);
}

1;
