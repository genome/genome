package Genome::ProcessingProfile::Command::Diff;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::Command::Diff {
    is => 'Command::V2',
    has => [
        from => {
            is => 'Genome::ProcessingProfile',
            shell_args_position => 1,
            doc => 'the first profile',
        },
        to => {
            is => 'Genome::ProcessingProfile',
            shell_args_position => 2,
            doc => 'the second profile',
        },
    ],
};

sub help_brief {
    return 'compare two processing profiles';
}

sub help_synopsis {
    return <<EOS
genome processing-profile diff 1234 5678
EOS
}

sub help_detail {
    return <<EOS
Compare two processing profiles side-by-size.
EOS
}

sub execute {
    my $self = shift;

    my $from = $self->from;
    my $to = $self->to;

    my %from = map { $_->name => $_->value_id } $from->params();
    my %to = map { $_->name => $_->value_id } $to->params();
    my %both = (%from, %to);

    for my $name (sort keys %both) {
        my $before = $from{$name};
        my $after = $to{$name};
        printf("%50s: %s\t%s", $name, $before, $after);
        if ($before ne $after) {
            print "\t******";
        }
        print "\n";
    }

    return 1;
}

1;

