package Genome::SoftwareResult::Command::FillLookupHashes;

use strict;
use warnings;

use Genome;
use Genome::Utility::List "in";


class Genome::SoftwareResult::Command::FillLookupHashes {
    is => 'Command::V2',

    has_optional => [
        commit_size => {
            is => 'Number',
            default_value => '100',
            doc => 'Commit to the UR::Context with this many ' .
                'calculated hashes at a time.',
        },
        class_names => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Fill or replace the lookup_hash values for ' .
                'all SoftwareResults of this class.',
        },
        blacklist_class_names => {
            is => 'Text',
            is_many => 1,
            doc => 'No objects of these classes will have hashes filled.',
        },
        num_cpus => {
            is => 'Number',
            default_value => 1,
            doc => 'Use this many processors to fill.',
        },
        replace_existing => {
            is => 'Boolean',
            doc => 'If true, replace existing lookup hashesh with new calculation.A,'
        },
    ],
};

sub execute {
    my $self = shift;

    my @pids;
    for (my $i = 0; $i < $self->num_cpus; $i++) {
        my $pid = UR::Context::Process->fork();
        unless (defined($pid)) {
            die $self->error_message("Failed to fork on subprocess number $i");
        }

        if ($pid == 0) {
            # child
            $self->fill_hashes($self->num_cpus, $i);
            exit;
        } else {
            # parent
            push(@pids, $pid);
        }
    }

    for my $pid (@pids) {
        waitpid($pid, 0);
    }

    return 1;
}

sub fill_hashes {
    my ($self, $base, $remainder) = @_;

    my @blacklist = $self->blacklist_class_names();
    my $i = 0;
    my $iterator = $self->_get_iterator();

    while (my $sr = $iterator->next()) {
        next unless $sr->id % $base == $remainder;
        next if (in($sr->class, @blacklist));

        eval {
            $sr->lookup_hash($sr->calculate_lookup_hash());
            $i++;
        };
        if ($@ =~ m/Can't locate object method/) {
            $sr->lookup_hash("bad class " . $sr->id);
            $i++;
        }

        if ($i >= $self->commit_size) {
            $self->status_message("Committing chunk of $i rows");
            $i = 0;

            # Destroy iterator before commit, so we don't get a warning when we reset it below.
            undef $iterator;
            UR::Context->commit();

            $iterator = $self->_get_iterator();
        }
    }

    if ($i) {
        $self->status_message("Committing chunk of $i rows");

        # Destroy iterator before commit, so we don't get a warning when we reset it below.
        undef $iterator;
        UR::Context->commit();

        $iterator = $self->_get_iterator();
    }
}

sub _get_iterator {
    my $self = shift;

    my %params;
    unless ($self->replace_existing) {
        $params{lookup_hash} = undef;
    }

    if (defined($self->class_names) and @{$self->class_names}) {
        if (@{$self->class_names} == 1) {
            $params{class_name} = @{$self->class_names}[0];
        } else {
            $params{'class_name in'} = $self->class_names;
        }
    }

    return Genome::SoftwareResult->create_iterator(%params);
}


1;
