package Finishing::Assembly::AGP::Reader;

use strict;
use warnings;

use base qw(Finfo::Reader);

use Finishing::Assembly::AGP::Utils;

sub _next
{
    my $self = shift;

    my $line;
    do
    {
        $line = $self->io->getline;
        return unless $line;
        chomp $line;
    } until defined $line and $line =~ /\w/;

    return Finishing::Assembly::AGP::Utils->instance->string_to_agp_hashref($line);
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Finishing/Assembly/AGP/Reader.pm $
#$Id: Reader.pm 30518 2007-11-30 22:45:36Z ebelter $
