#!/usr/bin/env genome-perl

use strict;
use lib '/gscuser/eclark/lib';
use above 'PAP';

use POE;
use POE::Component::IKC::Client;

my $instance_id = 88;

our $session = POE::Component::IKC::Client->spawn( 
    ip=>'blade2-4-2', 
    port=>13425,
    name=>'Controller',
    on_connect=>\&__build
);

POE::Kernel->run();

sub __build {
    our $controller = POE::Session->create(
        inline_states => {
            _start => sub {
                my ($kernel, $heap) = @_[KERNEL, HEAP];
                $kernel->alias_set("controller");
                $kernel->post('IKC','publish','controller',
                    [qw(got_plan_id got_instance_id complete error)]
                );

                my $kernel_name = $kernel->ID;

                $_[KERNEL]->post(
                    'IKC','call',
                    'poe://UR/workflow/resume',
                    [ 
                        $instance_id,
                        "poe://$kernel_name/controller/complete",
                        "poe://$kernel_name/controller/error"
                    ],
                    'poe:got_instance_id'
                );  
            },
            got_instance_id => sub {
                my ($kernel, $id) = @_[KERNEL, ARG0];
                print "Instance: $id\n";
#                $kernel->post('IKC'=>'shutdown');
            },
            complete => sub {
                my ($kernel, $arg) = @_[KERNEL, ARG0];
                my ($id, $instance, $execution) = @$arg;

                print "Complete: $id\n";
                $kernel->post('IKC'=>'shutdown');
            },
            error => sub {
                my ($kernel, $arg) = @_[KERNEL, ARG0];
                my ($id, $instance, $execution) = @$arg;

                print "Error: $id\n";
                $kernel->post('IKC'=>'shutdown');
            }
        }
    );
}
