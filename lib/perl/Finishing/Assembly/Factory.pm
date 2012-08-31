package Finishing::Assembly::Factory; 

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

#schemas
use Finishing::Assembly::Ace::Schema;
use Finishing::Assembly::Consed::Schema;
use Finishing::Assembly::GSC::Schema;
use Finishing::Assembly::DBIx::Schema;
use Finishing::Assembly::Phd::Ball;
use Finishing::Assembly::Phd::Directory;
use Finishing::Assembly::Phd::FastaAndQualDB;
use Finishing::Assembly::Source::Schema;

#objects
use Finishing::Assembly::Assembly;
use Finishing::Assembly::AssembledRead;
use Finishing::Assembly::Chromosome;
use Finishing::Assembly::ChromosomeFirstScaffold;
use Finishing::Assembly::Contig;
use Finishing::Assembly::Gap;
use Finishing::Assembly::ImprovementCorrelation;
use Finishing::Assembly::Organism;
#use Finishing::Assembly::Project;  
use Finishing::Assembly::Iterator;  
use Finishing::Assembly::Scaffold;
use Finishing::Assembly::Tags;

#proxies
use Finishing::Assembly::Ace::Proxies;
use Finishing::Assembly::Consed::Proxies;
use Finishing::Assembly::DBIx::Proxies;
use Finishing::Assembly::DSProxies;
use Finishing::Assembly::GSC::Proxies;
use Finishing::Assembly::Source::Proxies;

#others
use Data::Dumper;
use Finfo::ClassUtils 'class';
#use Finishing::Assembly::Project::Utils;
use XML::Simple;

#- ATTRIBUTES -#
my %name :name(name:r)
    :isa(string);
my %schema :name(schema:r)
    :isa([ 'object', ]);

#- CONFIGURATION -#
my $conf;
{
    my $conf_file = '/gsc/scripts/share/Assembly/factory-conf.xml';
    my $fh = IO::File->new('<' . $conf_file);
    __PACKAGE__->fatal_msg("Can't open conf file ($conf_file): $!")
        and return unless $fh;

    my $xs = XML::Simple->new
    (
        rootname => 'configuration',
        KeyAttr =>
        {
            schema => 'name',
            default => 'name',
            dbi => 'name',
            db => 'name',
        },
        ForceArray => [qw/ schema default dbi db on_connect_do /],
        ContentKey => '-content',
        VarAttr => 'name',
    );

    $conf = $xs->XMLin($fh);
    #print Dumper($conf);
}

sub available_schemas
{
    return sort { $a cmp $b } keys %{ $conf->{schema} };
}

sub available_dbs
{
    return sort { $a cmp $b } keys %{ $conf->{db} };
}

sub default_db
{
    return $conf->{default}->{db};
}

sub default_file_system_path
{
    return $conf->{default}->{file_system_path};
}

#- FACTORY IS A INSTANCE MANAGER AND PROXY FOR SCHEMA -#
#sub AUTOMETHOD
use Class::AutoloadCAN;
sub CAN {
    my ($starting_class, $requested_method, $self, @args) = @_;

    #my ($self, $id, @args) = @_;
    #my $requested_method = $_;

    return unless $self->schema->can($requested_method);

    # put this in conf??
    my $schema_methods = 
    {
        'Finishing::Assembly::Ace::Schema' =>
        {
            get_assembly => 'assembly',
        },
        'Finishing::Assembly::Consed::Schema' =>
        {
            get_assembly => 'assembly',
        },
        'Finishing::Assembly::DBIx::Schema' =>
        {
            create_project => 'project',
            get_project => 'project',
            get_or_create_project => 'project',
            get_organism =>  'organism',
            get_or_create_organism => 'organism',
            create_organism => 'organism',
            organisms => 'organism_iterator',
            get_assembly => 'assembly',
            #get_project => 'project',
        },
        'Finishing::Assembly::GSC::Schema' =>
        {
            get_project => 'project',
            create_project => 'project',
        },
        'Finishing::Assembly::Phd::Ball' =>
        {
            phd_names => 'null',
            phd_exists => 'null',
            get_phd => 'assembled_read',
            get_latest_phd => 'assembled_read',
            latest_phd_name => 'null',
            latest_phd_file => 'null',
            next_phd_name => 'null',
            next_phd_file => 'null',
        },
        'Finishing::Assembly::Phd::Directory' =>
        {
            phd_names => 'null',
            phd_exists => 'null',
            get_phd => 'assembled_read',
            get_latest_phd => 'assembled_read',
            latest_phd_name => 'null',
            latest_phd_file => 'null',
            next_phd_name => 'null',
            next_phd_file => 'null',
        },
        'Finishing::Assembly::Phd::FastaAndQualDB' =>
        {
            phd_names => 'null',
            phd_exists => 'null',
            get_phd => 'assembled_read',
            get_latest_phd => 'assembled_read',
            latest_phd_name => 'null',
            latest_phd_file => 'null',
            next_phd_name => 'null',
            next_phd_file => 'null',
        },
        'Finishing::Assembly::Source::Schema' =>
        {
            create_project => 'project',
            get_project => 'project',
            get_or_create_project => 'project',
            create_contig => 'contig',
            copy_contig => 'contig',
            create_consensus_tag => 'consensus_tag',
            copy_consensus_tag => 'consensus_tag',
            create_read_tag => 'read_tag',
            copy_read_tag => 'read_tag',
            create_assembly_tag => 'assembly_tag',
            copy_assembly_tag => 'assembly_tag',
			create_assembled_read => 'assembled_read',
			copy_assembled_read => 'assembled_read',
        },
    };

    my $schema_class = class($self->schema);
    $self->fatal_msg
    (
        "Valid schema method ($requested_method), but no information for wrapping source"
    ) unless exists $schema_methods->{$schema_class}->{$requested_method};

    my $type = $schema_methods->{$schema_class}->{$requested_method};
    my $wrap_method = sprintf
    (
        '_wrap_%s',
        $type,
    ); 
    $self->fatal_msg("Invalid source type ($type)") unless $self->can($wrap_method);

    return sub
    {
        my $source = $self->schema->$requested_method(@args); 
        return unless $source;
        return $self->$wrap_method($source);
    }
}

no warnings 'redefine';
sub new
{
    __PACKAGE__->fatal_msg("Use connect to create");
}
use warnings;
no warnings 'reserved';

my %instances;
sub connect
{
    my ($self, $db, $srcs, $connect_attr, $connect_mod) = @_;

    $db = $self->default_db unless defined $db;
    
    Finfo::Validate->validate
    (
        attr => 'db',
        value => $db,
        isa => [ 'in_list', $self->available_dbs ],
        msg => 'fatal',
    );
   
    my $schema_class = 'Finishing::Assembly::' . $conf->{db}->{$db}->{schema};
    my $dbi = $conf->{db}->{$db}->{dbi};
    my $string = '';
    if ( $srcs )
    {
        if ( ref($srcs) )
        {
            $dbi = $srcs;
            $string .= $srcs->{acefile};
        }
        else
        {
            $dbi .= $srcs;
            $string = $srcs;
        }
    }

    my $name = sprintf('%s=%s-%s', $schema_class, $db, $string);

    # Return schema, if already connected 
    return $instances{$name} if exists $instances{$name};

    my $user = $conf->{db}->{$db}->{user} || '';
    my $pw = $conf->{db}->{$db}->{pw} || '';
    if ( exists $conf->{db}->{$db}->{connect_attr} )
    {
        foreach my $attr ( keys %{ $conf->{db}->{$db}->{connect_attr} } )
        {
            next if exists $connect_attr->{$attr};
            $connect_attr->{$attr} = $conf->{db}->{$db}->{connect_attr}->{$attr};
        }
    }
    if ( exists $conf->{db}->{$db}->{connect_mod} )
    {
        foreach my $mod ( keys %{ $conf->{db}->{$db}->{connect_mod} } )
        {
            next if exists $connect_mod->{$mod};
            $connect_mod->{$mod} = $conf->{db}->{$db}->{connect_mod}->{$mod};
        }
    }

    #$self->info_msg(Dumper([ $dbi, $user, $pw, $connect_attr, $connect_mod ]));
    
    my $schema = $schema_class->connect
    (
        $dbi,
        $user,
        $pw,
        $connect_attr,
        $connect_mod,
    );
    $self->fatal_msg("Cannot connect to schema ($schema_class)") unless $schema;

    my $new_instance = Finfo::Std::new
    (
        __PACKAGE__, 
        name => $name,
        schema => $schema,
    );

    $instances{$name} = $new_instance;

    return $new_instance;
}

sub disconnect
{
    my $self = shift;

    delete $instances{ $self->name };


    $self->schema->disconnect if $self->schema->can('disconnect');
    $self = undef;

    return 1;
}

sub rollback
{
    my $self = shift;

    return 1 unless $self->schema->can('rollback');

    return $self->schema->rollback;
}

sub commit
{
    my $self = shift;

    return 1 unless $self->schema->can('commit');

    return $self->schema->commit;
}

#- CONSTRUCTOR -#
sub _object_constructor
{
    my ($self, %p) = @_;
    
    my $source = delete $p{source}
        or $self->fatal_msg("No object source to wrap", { caller_level => 1 });
    my $type = delete $p{type}
        or $self->fatal_msg("No object type to wrap", { caller_level => 1 });

    my $wrap_method = sprintf
    (
        '_wrap_%s',
        $type,
    ); 
    $self->fatal_msg("Invalid source type ($type)") unless $self->can($wrap_method);

    return $self->$wrap_method($source);
}

#- WRAP -#
sub _wrap_null
{
    my ($self, $source) = @_;

    return $source;
}
    
sub _wrap_organism
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('organism', $source)
        or return $source;

    return Finishing::Assembly::Organism->new(proxy => $proxy);
}

sub _wrap_organism_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor => sub 
        {
            my $source = shift; 
            return $self->_wrap_organism($source);
        },
    );
}

sub _wrap_assembly
{
    my ($self, $source) = @_; 
    
    my $proxy = $self->_get_proxy('assembly', $source)
        or return $source;

    return Finishing::Assembly::Assembly->new(proxy => $proxy);
}

sub _wrap_assembly_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor => sub
        { 
            my $source = shift; 
            return $self->_wrap_assembly($source);
        },
    );
}

sub _wrap_assembly_tag
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('assembly_tag', $source)
        or return $source;

    return Finishing::Assembly::AssemblyTag->new(proxy => $proxy);
}

sub _wrap_assembly_tag_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor =>  sub
        {
            my $source = shift; 
            return $self->_wrap_assembly_tag($source);
        },
    );
}

sub _wrap_consensus_tag
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('consensus_tag', $source)
        or return $source;

    return Finishing::Assembly::ConsensusTag->new(proxy => $proxy);
}

sub _wrap_consensus_tag_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor =>  sub
        {
            my $source = shift; 
            return $self->_wrap_consensus_tag($source);
        },
    );
}

sub _wrap_read_tag
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('read_tag', $source)
        or return $source;

    return Finishing::Assembly::ReadTag->new(proxy => $proxy);
}

sub _wrap_read_tag_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor =>  sub
        {
            my $source = shift; 
            return $self->_wrap_read_tag($source);
        },
    );
}

sub _wrap_scaffold
{
    my ($self, $source) = @_; 

    my $proxy = $self->_get_proxy('scaffold', $source)
        or return $source;
    
    return Finishing::Assembly::Scaffold->new(proxy => $proxy);
}

sub _wrap_scaffold_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor =>  sub
        {
            my $source = shift; 
            return $self->_wrap_scaffold($source);
        },
    );
}

sub _wrap_contig
{
    my ($self, $source) = @_;
    
    my $proxy = $self->_get_proxy('contig', $source)
        or return $source;

    return Finishing::Assembly::Contig->new(proxy => $proxy);
}

sub _wrap_contig_iterator
{
    my ($self, $rs) = @_;
    
    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor => sub 
        {
            my $source = shift; 
            return $self->_wrap_contig($source);
        },
    );
}

sub _wrap_gap
{
    my ($self, $source) = @_;
    
    my $proxy = $self->_get_proxy('gap', $source)
        or return $source;

    return Finishing::Assembly::Gap->new(proxy => $proxy);
}

sub _wrap_gap_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor => sub 
        {
            my $source = shift; 
            return $self->_wrap_gap($source);
        },
    );
}

sub _wrap_assembled_read
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('assembled_read', $source)
        or return $source;

    return Finishing::Assembly::AssembledRead->new(proxy => $proxy);
}

sub _wrap_assembled_read_iterator
{
    my ($self, $rs) = @_;
    
    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor => sub
        {
            my $source = shift; 
            return $self->_wrap_assembled_read($source);
        },
    );
}

sub _wrap_chromosome
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('chromosome', $source)
        or return $source;

    return Finishing::Assembly::Chromosome->new(proxy => $proxy);
}

sub _wrap_chromosome_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator => $rs,
        object_constructor => sub
        {
            my $source = shift;
            return $self->_wrap_chromosome($source);
        },
    );
}

sub _wrap_improvement_correlation
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('improvement_correlation', $source)
        or return $source;

    return Finishing::Assembly::ImprovementCorrelation->new(proxy => $proxy);
}

sub _wrap_improvement_correlation_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator =>$rs, 
        object_constructor => sub 
        {
            my $source = shift;
            return $self->_wrap_improvement_correlation($source);
        },
    );
}

sub _wrap_chromosome_first_scaffold
{
    my ($self, $source) = @_; 

    my $proxy = $self->_get_proxy('chromosome_first_scaffold', $source)
        or return $source;

    return Finishing::Assembly::ChromosomeFirstScaffold->new(proxy => $proxy);
}

sub _wrap_project
{
    my ($self, $source) = @_;

    my $proxy = $self->_get_proxy('project', $source)
        or return $source;
    
    return Finishing::Assembly::Project->new(proxy => $proxy);
}

sub _wrap_project_iterator
{
    my ($self, $rs) = @_;

    return Finishing::Assembly::Iterator->new
    (
        iterator =>$rs, 
        object_constructor => sub 
        {
            my $source = shift;
            return $self->_wrap_project($source);
        },
    );
}

#- PROXY -#
sub _get_proxy
{
    my ($self, $object_type, $source) = @_;

    my $object_class = join('', map { ucfirst } split(/_/, $object_type));
    my $source_type = ref($source);

    if ( $source_type =~ /^Finishing::Assembly::\w+$/ )
    {
        # this is the already wrapped behavior object
        return;
    }
    elsif ( $source_type =~ /^Finishing::Assembly::DBIx::Schema::/i )
    {
        my $proxy_class = sprintf('Finishing::Assembly::%sDSProxy', $object_class);
        return $proxy_class->new
        (
            source => $source,
            object_constructor => sub{ $self->_object_constructor(@_) }, 
        );
    }
    elsif ( $source_type =~ /^GSC::/ )
    {
        my $proxy_class = sprintf('Finishing::Assembly::GSC::%sProxy', $object_class);
        return $proxy_class->new
        (
            source => $source,
            object_constructor => sub{ $self->_object_constructor(@_) }, 
        );
    }
    else 
    {
        $source_type =~ /^Finishing::Assembly::(\w+)::/;
        my $proxy_class = sprintf('Finishing::Assembly::%s::%sProxy', $1, $object_class);
        return $proxy_class->new
        (
            source => $source,
            object_constructor => sub{ $self->_object_constructor(@_); }, 
        );
    }
}

1;

#$HeadURL$
#$Id$
