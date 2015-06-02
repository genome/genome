package Genome::Report::Command;
#:adukes check

use strict;
use warnings;

use Genome;
      
class Genome::Report::Command {
    is => 'Command',
        doc => 'work with reports',
        has => [ 
        report_directory => { 
            is => 'Text', 
            doc => 'Report directory.',
        },
    ],
};

#< Command >#
sub help_brief {
    my $class = ref($_[0]) || $_[0];
    return $class->__meta__->doc if not $class or $class eq __PACKAGE__;
    my ($func) = $class =~ /::(\w+)$/;
    return sprintf('%s a report', ucfirst($func));
}

sub help_detail {
    return help_brief(@_);
}

sub command_name {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name unless $class eq __PACKAGE__;
    return 'genome report';
}

sub command_name_brief {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name_brief unless $class eq __PACKAGE__;
    return 'report';
}
#<>#

#< Report >#
sub report {
    my $self = shift;

    unless ( $self->{_report} ) { 
        $self->{_report} = Genome::Report->create_report_from_directory($self->report_directory);
    }

    return $self->{_report};
}
#<>#

1;

#$HeadURL$
#$Id$
