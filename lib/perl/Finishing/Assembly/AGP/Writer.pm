package Finishing::Assembly::AGP::Writer;

use strict;
use warnings;

use base qw(Finfo::Writer);

use Finishing::Assembly::AGP::Utils;

sub _write_one
{
    my ($self, $agp) = @_;

    my $string = Finishing::Assembly::AGP::Utils->instance->agp_to_string($agp);

    return unless $string;
    
    return $self->io->print("$string\n");
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Finishing/Assembly/AGP/Writer.pm $
#$Id: Writer.pm 30518 2007-11-30 22:45:36Z ebelter $
