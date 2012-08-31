package BAP::Config;

use strict;
use warnings;

# There can be only one!
use Memoize;
memoize('new');

sub new {

    my ($class) = @_;


    my $self = { };

    bless $self, $class;


    $self->{_major} = 2;
    $self->{_minor} = 10;
    $self->{_activity_db_file} = '/gscmnt/temp212/info/annotation/BAP_db/mgap_activity.db';
    
    return $self;

}

sub major {

    my ($self) = @_;

    
    return $self->{_major};
    
}

sub minor {

    my ($self) = @_;


    return $self->{_minor};

}

sub version {

    my ($self) = @_;


    return join('.', $self->major(), $self->minor());
    
}

sub activity_db_file {

    my ($self) = @_;


    return join('.', $self->{_activity_db_file});

}

1;
