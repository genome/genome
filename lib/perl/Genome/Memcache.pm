package Genome::Memcache;

use strict;
use warnings;

use Cache::Memcached;

use UR;
use Genome;

class Genome::Memcache {
    is => 'UR::Singleton',
    doc => 'methods for accessing memcache',
    has => {
        memcache_server_location => {
            is => 'Text',
            default_value => $ENV{GENOME_SYS_SERVICES_MEMCACHE},
        },
        _memcache_server => {
            is => 'Cache::Memcached',
            is_transient => 1,
        },
        memcache_server => {
            calculate_from => ['_memcache_server', 'memcache_server_location'],
            calculate => q{
                 return $_memcache_server if $_memcache_server;
                 
                 $self->_memcache_server(
                    new Cache::Memcached {'servers' => [$memcache_server_location], 'debug' => 0, 'compress_threshold' => 10_000,}
                 );
                 
                 return $self->_memcache_server
            }
        },
    }
};

sub server {
    my ($class) = @_;
    my $server = $class->_singleton_object->memcache_server();
    return $server;
}

1;

