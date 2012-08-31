package Genome::Model::Tools::WuBlast::Xdformat::Verify;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::WuBlast::Xdformat::Verify {
    is => 'Genome::Model::Tools::WuBlast::Xdformat',
};

#< Standard command methods >#
sub help_brief {
    return "Verifies an xdformat database";
}

sub help_detail {
    return help_brief();
    return <<EOS
EOS
}

#< Additional Params #>
sub _additional_params_as_string {
    my $self = shift;

    return join(' ', $self->fasta_files);
}

#< Operation >#
sub _operation_name {
    return 'verify';
}

sub _operation_character {
    return 'V';
}

1;

#$HeadURL$
#$Id$
