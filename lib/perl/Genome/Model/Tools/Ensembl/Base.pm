package Genome::Model::Tools::Ensembl::Base;

use strict;
use warnings;

use Genome;

my $DEFAULT_VERSION = '58_37c';

class Genome::Model::Tools::Ensembl::Base {
    is  => 'Command::V2',
    is_abstract => 1,
    has => [
        use_version => {
            default_value => Genome::Model::Tools::Ensembl::Base->default_version,
        },
        species => {
            doc => 'The Ensembl species name',
            valid_values => ['homo_sapiens', 'mouse'],
            default_value => 'homo_sapiens',
        },
        database => {
            is_calculated => 1,
            calculate_from => ['use_version','species'],
            calculate => sub {
                my $version = shift;
                my $species = shift;
                return $species .'_core_'. $version;
            },
        },
        password => {
            default_value => Genome::Model::Tools::Ensembl::Base->default_password,
            doc => 'The password for the Ensembl MySQL server.',
            is_optional => 1,
        },
        username => {
            default_value => Genome::Model::Tools::Ensembl::Base->default_username,
            doc => 'The username for the Ensembl MySQL server.',
            is_optional => 1,
        },
        hostname => {
            default_value => Genome::Model::Tools::Ensembl::Base->default_hostname,
            doc => 'The location of the Ensembl MySQL server.',
            is_optional => 1,
        }
    ],
    has_optional => [
        _api_version => {},
        _registry_loaded => {},
        _slice_adaptor => {},
    ],
};

sub default_version {
    return $DEFAULT_VERSION;
}

sub default_hostname {
    return Genome::Config::get('db_ensembl_host');
}

sub default_username {
    return Genome::Config::get('db_ensembl_user');
}

sub default_password {
    return Genome::Config::get('db_ensembl_pass');
}

sub api_version {
    my $self = shift;
    unless (defined($self->_api_version)) {
        my $version = $self->use_version;
        my $api_version;
        if ($version =~ /^(\d+)\_.*/){
            $api_version = $1;
        } else {
            die('Could not determine Ensembl API version from Ensembl database version:'. $version);
        }
        $self->_api_version($api_version);
    }
    return $self->_api_version;
}

sub create_dbh {
    my $self = shift;
    my $dbi_string = 'dbi:mysql:database='.$self->database .';host='.$self->hostname;
    my @connect_params = ( $dbi_string, $self->username, $self->password, { PrintError => 1 } );
    my $dbh = DBI->connect( @connect_params );
    my $connection_attempts = 1;

    #Check for failures to reconnect... if the connection failed, try again a few times
    my $connection_attempt_limit = 3;
    my $sleep_time = 15;
    while(!(defined($dbh)) && ($connection_attempts <= $connection_attempt_limit)){
        $connection_attempts++;
        $self->status_message('DBI connect failed ... sleeping for '. $sleep_time .' and try attempt: '. $connection_attempts);
        sleep $sleep_time;
        $dbh = DBI->connect( @connect_params );
    }
    unless ($dbh) {
        $self->error_message('Failed to connect using DBI with params: '. Data::Dumper::Dumper(@connect_params));
        die($self->error_message);
    }
    return $dbh;
}

sub load_registry {
    my $self = shift;
    unless ($self->_registry_loaded) {
        Bio::EnsEMBL::Registry->load_registry_from_db(
            -host => $self->hostname,
            -user => $self->username,
            -pass => $self->password,
            -db_version => $self->api_version
        );
        $self->_registry_loaded(1);
    }
    return 1;
}

sub slice_adaptor {
    my $self = shift;
    unless ($self->_registry_loaded) {
        $self->load_registry;
    }
    unless ($self->_slice_adaptor) {
        my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
            $self->species,
            'Core',
            'Slice'
        );
        $self->_slice_adaptor($slice_adaptor);
    }
    return $self->_slice_adaptor;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    my $bioperl_path = '/gsc/lib/perl5/bioperl/1.6.0/lib/perl5/';
    unshift(@INC, $bioperl_path);
    require Bio::EnsEMBL::DBSQL::DBAdaptor;

    return $self;
}

1;

