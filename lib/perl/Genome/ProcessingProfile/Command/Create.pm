package Genome::ProcessingProfile::Command::Create;

use strict;
use warnings;
use Genome;

class Genome::ProcessingProfile::Command::Create {
    is => 'Command::SubCommandFactory',
    doc => 'create a processing profile',
};

sub _sub_commands_from { 'Genome::Model' }
sub _target_base_class { 'Genome::ProcessingProfile' }
sub _sub_commands_inherit_from { 'Genome::ProcessingProfile::Command::Create::Base' }

sub _command_subclass {
    my ($self, $class_name) = @_;
    my ($subclass_name) = $class_name =~ /^Genome::ProcessingProfile::Command::Create::(.*)$/;
    return $subclass_name;
}

sub _processing_profile_class {
    my ($self, $class_name) = @_;
    my $subclass_name = $self->_command_subclass($class_name);
    my $processing_profile_class = join('::', 'Genome::ProcessingProfile', $subclass_name);
    return $processing_profile_class;
}

sub _profile_string {
    my ($self, $class_name) = @_;
    my $subclass_name = $self->_command_subclass($class_name);
    my $profile_string = Genome::Utility::Text::camel_case_to_string($subclass_name);
    $profile_string =~ s/:://;
    return $profile_string;
}

sub _build_sub_command {
    my ($self, $class_name, @inheritance) = @_;

    class {$class_name} {
        is => \@inheritance,
        has => [
            $self->_sub_commands_inherit_from->_properties_for_class($self->_processing_profile_class($class_name)),
        ],
        doc => "Create a new profile for " . $self->_profile_string($class_name),
    };

    return $class_name;
}
