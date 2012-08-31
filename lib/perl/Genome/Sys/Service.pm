package Genome::Sys::Service;
use strict;
use warnings;
use Genome;
use LWP::UserAgent;

class Genome::Sys::Service {
    is_abstract => 1,
    table_name => "NA", # you shouldn't have to have a value here for custom data sources: bug
    subclassify_by => '_subclass_name',
    has => [
        name => { is => 'Text' },
        _subclass_name => { calculate => q| 'Genome::Sys::Service::' . $name |, calculate_from => ['name'], },
    ],
    doc         => 'a service used by the genome system',

};

sub __display_name__ {
    my $self = shift;
    return $self->name . ' (' . $self->host . ')';
}

sub __load__ {
    my ($class, $bx) = @_;

    if ($class eq __PACKAGE__) {
        my $path = __FILE__;
        $path =~ s/.pm$//;
        my @files = glob("$path/*.pm");

        my @rows;
        for my $file (@files) {
            my $ext = File::Basename::basename($file);
            my $name = $ext;
            $name =~ s/.pm$//;
            my $class_name = __PACKAGE__ . '::' . $name;
            my $meta = UR::Object::Type->get($class_name);
            next unless ($class_name->isa(__PACKAGE__));
            unless ($class_name->isa("UR::Singleton")) {
                die "class $class_name should inherit from UR::Singleton!";
            }
            my @o = $class_name->get();
            push @rows,[$class_name, $name];
        }

        return ['id','name'],\@rows;
    }
    else {
        my $name = $class;
        $name =~ s/Genome::Sys::Service:://;
        return ['id','name'], [[$class, $name]];
    }
}

sub status {
    my $self = shift;

    return 'unknown' if !$self->url();

    my $ua = LWP::UserAgent->new();
    $ua->cookie_jar({ file => "/tmp/nomnomcookies.txt" });
    $ua->timeout(10);

    if ($ua->get($self->url())->is_success) {
        return 'running';
    } else {
        return 'stopped';
    }
}

sub pid_status {
    my $self = shift;
    my $pid_name = $self->pid_name();

    my $cmd = 'ssh ' . $self->host() . " 'ps aux | grep \"$pid_name\" | grep -v grep'";
    my $o = `$cmd`;

    if ($o) {
        return 'running';
    } else {
        return 'stopped';
    }
}

# dynamically create an accessor method for the service API, and emit a warning whenever it is not overridden in subclasses
# TODO: it would be better to check this when initializing the sub-class than when using the accessor ...but this
# may be easier to work with while we are still getting the API stabilized.
my @methods = (qw/host restart_command stop_command log_path pid_name url/);
for my $method (@methods) {
    my $sub = sub {
        my $self = shift;
        warn "service class " . $self->class . " does not implement static method $method!";
        return undef;
    };
    no strict 'refs';

    Sub::Install::install_sub(
        {
            code => $sub,
            into => __PACKAGE__,
            as => $method
        }
    );
}

1;

