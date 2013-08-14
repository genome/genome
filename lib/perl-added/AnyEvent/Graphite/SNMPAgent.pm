package AnyEvent::Graphite::SNMPAgent;

use warnings;
use strict;

use AnyEvent::Graphite;
use AnyEvent::SNMP;
use AnyEvent::DNS;
use Net::SNMP;

sub new {
    my($class, %args) = @_;

    # the rest of the args will get passed to the object
    # so remove host and port as we'll possibly set them manually here
    my $graphite = AnyEvent::Graphite->new(host => $args{host}, port => $args{port});
    delete $args{host};
    delete $args{port};
   
    # the hosts we'll be storing metrics in
    my $hosts = {};
    my $interval = delete $args{interval} || 60;

    my $timeout = delete $args{timeout} || 3;

    my $snmp_version = delete $args{snmp_version} || "1";

    # need to know about $self before defining the callback just below this
    my $self = bless {
        graphite => $graphite,
        hosts => $hosts,
        snmp_version => $snmp_version,
        %args,
    }, $class;

       
    # start the timer running
    $self->{timer} = AE::timer 0, $interval, sub { $self->gather_metrics; }; 

    return $self;
}

sub add_snmp {
    my ($self, %arg) = @_;
    # toss the arguments into the (per-host) queue
    # look up the host with AnyEvent::DNS so we can be async even here
    AnyEvent::DNS::a $arg{'host'}, sub {
        my (@addrs) = @_;
        push(@{$self->{hosts}{$addrs[0]}}, \%arg); 
    };
}

# this is called by $self->{timer} every $interval seconds
sub gather_metrics {
    my ($self) = @_;
    for my $host (keys %{$self->{hosts}}) {

        # skip any hosts that did not resolve
        next unless $host;

        # steal a community string from the first item in the list. They should all be the same
        my $community = $self->{hosts}{$host}[0]{community} || "public";
        my $session = Net::SNMP->session(
            -hostname => $host,
            -community => $community,
            -nonblocking => 1,
            -version => $self->{snmp_version}
        );

        # in this kind of context it's not clear what would be better to do with errors, here.
        next unless $session;

        # if you don't set a timeout, you can fill up your queues of outstanding processes
        # protects against 'lame' servers
        $session->timeout($self->{timeout});

        for my $metric (@{$self->{hosts}{$host}}) {
            $session->get_request( 
                -varbindlist => [$metric->{oid}], 
                -callback => sub {
                    my ($session) = @_;
                    my $result = $session->var_bind_list();
                    my $value = $result->{$metric->{oid}};
                    if(!defined($value)) {
                        warn "undefined reply from $host" . ":" . "$metric->{oid}: " . $session->error();
                    } else {
                        if($metric->{filter}) {
                            $value = $metric->{filter}->($value);
                        }
                        $self->{graphite}->send($metric->{graphite_key}, $value);
                    }
                });
        }
    }

}

=head1 NAME

AnyEvent::Graphite::SNMPAgent - An SNMP agent which does non-blocking streaming of data from an SNMP server

=head1 SYNOPSIS
    
    my $agent = AnyEvent::Graphite::SNMPAgent->new(
        host => '127.0.0.1',
        port => '2003',
        interval => 60,
        timeout => 5,
    );

    'host' and 'port' are for the graphite server
    'interval' is how many seconds should elapse between each time we try to fire off the queries.
    If you need multiple intervals create one AE::Graphite::SNMPAgent instance per set of metrics


    $agent->add_snmp(host => $host, oid => $oid, community => $community, graphite_key => $key, filter => sub { ... });

    print "Running forever. CTRL-C to interrupt\n";
    AnyEvent->condvar->recv;


=head1 AUTHOR

Joshua Barratt, C<< <josh at mediatemple.net> >>

=head1 COPYRIGHT & LICENSE

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of AnyEvent::Graphite
