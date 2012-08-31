package Genome::Sys::Service::Set::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::Sys::Service::Set::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Sys::Service',
                    aspects => [
                        'id',
                        'name',
                        'host',
                        'restart_command',
                        'stop_command',
                        'log_path',
                        'status',
                        'pid_status',
                        'pid_name',
                        'url',
                    ],
                },
            ]
        }
    ],
};

1;
