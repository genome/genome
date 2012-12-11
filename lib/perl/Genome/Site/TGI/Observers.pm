package Genome::Site::TGI::Observers;

use strict;
use warnings;

my %observer_class_name = (
    'Command::V1' => 'Genome::Site::TGI::Observers::Command',
    'Genome::ModelGroup' => 'Genome::Site::TGI::Observers::ModelGroup',
    'Genome::Project' => 'Genome::Site::TGI::Observers::Project',
);

UR::Object::Type->add_observer(
    aspect => 'load',
    callback => sub {
        my $meta = shift;
        my $class_name = $meta->class_name;
        if (exists $observer_class_name{$class_name}) {
            my $observer_class_name = $observer_class_name{$class_name};
            my $observer_class_path = file_name($observer_class_name);
            require $observer_class_path;
            $observer_class_name->import($class_name);
        }
        die $@ if $@;
    },
);

sub file_name {
    my $class_name = shift;
    $class_name =~ s/::/\//g;
    $class_name .= '.pm';
    return $class_name;
}

1;
