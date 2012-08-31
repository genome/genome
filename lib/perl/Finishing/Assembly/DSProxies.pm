package Finishing::Assembly::DSProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Proxy';

use Finfo::ClassUtils 'class';
use Finishing::Assembly::DBIx::Proxies;

my %dbic_proxy :name(_dbic_proxy:p);

sub START
{
    my $self = shift;

    my ($source_class) = class( $self->source ) =~ /::(\w+)$/;
    my $dbic_proxy_class = sprintf('Finishing::Assembly::DBIx::%sProxy', $source_class);
    $self->_dbic_proxy
    (
        $dbic_proxy_class->new
        (
            source => $self->source, 
            object_constructor => $self->object_constructor,
        )
    );
    
    return 1;
}

sub get_method
{
    my ($self, $method, @args) = @_;

    # Priority: 
    # 1 defined method in ds proxy
    # 2 dbic_proxy get_method
    # 3 check methods in ds proxy (ace methods)
    
    if ( $self->can($method) )
    {
        return sub{ $self->$method(@args) };
    }

    if ( my $sub = $self->_dbic_proxy->get_method($method, @args) )
    {
        return $sub;
    }

    return $self->_check_methods($method, @args);
}

#############################################################################

package Finishing::Assembly::OrganismDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::AssemblyDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ScaffoldDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ContigDSProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DSProxy';

my %ace_contig :name(_ace_contig:p);

sub START
{
    my $self = shift;

    my $acefile;
    if ( $acefile = $self->acefile and -e $acefile)
    {
        my $ace_factory = Finishing::Assembly::Factory->connect('ace', $acefile);
        my $assembly = $ace_factory->get_assembly;
        my $ace_contig = $assembly->get_contig( $self->source->id );
        $self->_ace_contig($ace_contig);
    }
    return 1;
}

sub methods : CUMULATIVE
{
    return (qw/ ace /);
}

sub methods_for_ace
{
    return 
    (
        base_segments => undef,
    );  
}

sub ace
{
    my ($self, %p) = @_;

    my $method = $p{method};
    my $val = $self->_ace_contig->$method( @{ $p{args} });

    return unless defined $val;

    return $val unless $p{params};

    return $self->_construct_objects($p{params}, $val);
}

sub acefile
{
    my $self = shift;
    
    my $path = $self->source->assembly->file_path;

    return unless $path;
    
    return sprintf('%s/%d.ace', $path, $self->source->id);
}

sub tags
{
    my ($self, $new_tags) = @_;

    # dbic
    my $method = $self->_dbic_proxy->get_method('consensus_tags', $new_tags);
    my @tags = $method->() if $method;

    #ace
    if ( $self->_ace_contig and ( not @tags or $new_tags) )
    {
        my $ace_tags = $self->_ace_contig->tags($new_tags); 
        @tags = @$ace_tags unless @tags;
    }

    return \@tags;
}

sub base_string
{
    my ($self, $new_bases) = @_;

    # dbic
    my $method = $self->_dbic_proxy->get_method('base_string', $new_bases);
    my $bases = $method->() if $method;

    #ace
    if ( $self->_ace_contig and ( not $bases or $new_bases) )
    {
        $self->info_msg("Getting base_tring from ace");
        my $ace_bases = $self->_ace_contig->base_string($new_bases); 
        $bases = $ace_bases unless $bases;
    }

    return $bases;
}

sub qualities
{
    my ($self, $new_quals) = @_;

    # dbic
    my $method = $self->_dbic_proxy->get_method('qualities', $new_quals);
    my $quals = $method->() if $method;

    #ace
    if ( $self->_ace_contig and ( not $quals or $new_quals) )
    {
        $self->info_msg("Getting qualities from ace");
        my $ace_quals = $self->_ace_contig->qualities($new_quals); 
        $quals = $ace_quals unless $quals;
    }

    return $quals
}

sub length
{
    my $self = shift;

    # dbic
    my $method = $self->_dbic_proxy->get_method('length');
    my $length = $method->() if $method;

    # ace
    if ( not $length and $self->_ace_contig )
    {
        $self->info_msg("Getting length from ace");
        $length = $self->_ace_contig->length; 
    }
    
    return $length
}

#############################################################################

package Finishing::Assembly::AssembledReadDSProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DSProxy';

use Data::Dumper;

my %ace_read :name(_ace_read:p);

sub START
{
    my $self = shift;

    my $acefile;
    if ( $acefile = $self->acefile and -e $acefile)
    {
        my $ace_factory = Finishing::Assembly::Factory->connect('ace', $acefile);
        my $assembly = $ace_factory->get_assembly;
        my $ace_read = $assembly->get_assembled_read( $self->source->name );
        $self->_ace_read($ace_read);
    }

    return 1;
}

sub methods : CUMULATIVE
{
    return (qw/ ace /);
}

sub methods_for_ace
{
    return 
    (
        qual_clip_start => undef,
        qual_clip_stop => undef,
        align_clip_start => undef,
        align_clip_stop => undef,
        'time' => undef,
        chromat_file => undef,
        phd_file => undef,
        info_count => undef,
        tags => 'read_tag',
    );  
}

sub ace
{
    my ($self, %p) = @_;

    my $method = $p{method};
    my $val = $self->_ace_read->$method( @{ $p{args} });

    return unless defined $val;

    return $val unless $p{params};

    return $self->_construct_objects($p{params}, $val);
}


sub acefile
{
    my $self = shift;
    
    my $path = $self->source->assembly->file_path;
    return unless $path;
    
    return sprintf('%s/%d.ace', $path, $self->source->contig->id);
}

sub tags
{
    my ($self, $new_tags) = @_;

    # dbic
    my $method = $self->_dbic_proxy->get_method('read_tags', $new_tags);
    my @tags = $method->() if $method;

    #ace
    if ( $self->_ace_read and ( not @tags or $new_tags) )
    {
        my $ace_tags = $self->_ace_read->tags($new_tags); 
        @tags = @$ace_tags unless @tags;
    }

    return \@tags;
}

sub base_string
{
    my ($self, $new_bases) = @_;

    # dbic
    my $method = $self->_dbic_proxy->get_method('base_string', $new_bases);
    my $bases = $method->() if $method;

    #ace
    if ( $self->_ace_read and ( not $bases or $new_bases) )
    {
        my $ace_bases = $self->_ace_read->base_string($new_bases); 
        $bases = $ace_bases unless $bases;
    }

    return $bases;
}

sub qualities
{
    my ($self, $new_quals) = @_;

    # dbic
    my $method = $self->_dbic_proxy->get_method('qualities', $new_quals);
    my $quals = $method->() if $method;

    #ace
    if ( $self->_ace_read and ( not $quals or $new_quals) )
    {
        my $ace_quals = $self->_ace_read->qualities($new_quals); 
        $quals = $ace_quals unless $quals;
    }

    return $quals
}

sub length
{
    my $self = shift;

    # dbic
    my $method = $self->_dbic_proxy->get_method('length');
    my $length = $method->() if $method;

    # ace
    if ( not $length and $self->_ace_read )
    {
        $self->info_msg("Getting length from ace");
        $length = $self->_ace_read->length; 
    }
    
    return $length
}

#############################################################################

package Finishing::Assembly::ChromosomeFirstScaffoldDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ChromosomeDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ImprovementCorrelationDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ProjectDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ConsensusTagDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::ReadTagDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::AssemblyTagDSProxy;

use base 'Finishing::Assembly::DSProxy';

#############################################################################

package Finishing::Assembly::GapDSProxy;

use base 'Finishing::Assembly::DSProxy';
1;

#$HeadURL$
#$Id$
