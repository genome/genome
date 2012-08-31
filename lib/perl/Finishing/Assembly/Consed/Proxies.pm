################
#- READ PROXY -#
################

package Finishing::Assembly::Consed::AssembledReadProxy;

use base 'Finishing::Assembly::Proxy';

use Data::Dumper;

# TODO sync shared ace and phd attributes like base_string and time?
sub methods : CUMULATIVE
{
    return (qw/ ace phd /);
}

sub methods_for_ace : RESTRICTED
{
    return 
    (
        name => undef,
        base_string => undef,
        position => undef,
        complemented => undef,
        tags => 'read_tag',
        qual_clip_start => undef,
        qual_clip_stop => undef,
        align_clip_start => undef,
        align_clip_stop => undef,
        'time' => undef,
        chromat_file => undef,
        phd_file => undef,
        chem => undef,
        dye => undef,
        info_count => undef,
        'length' => undef,
    );
}

sub ace : RESTRICTED
{
    my ($self, %p) = @_;

    my $method = $p{method};

    return $self->source->ace_source->$method;
}

sub methods_for_phd : RESTRICTED
{
    return 
    (
        qualities => undef,
    );
}

sub phd : RESTRICTED
{
    my ($self, %p) = @_;

    my $method = $p{method};

    return $self->source->phd_source->{$method};
}

1;

#$HeadURL$
#$Id$
