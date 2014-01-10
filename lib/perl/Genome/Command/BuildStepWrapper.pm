package Genome::Command::BuildStepWrapper;
use strict;
use warnings;
use Genome;

class Genome::Command::BuildStepWrapper {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => [
        wrapper_build       => { is => 'Genome::Model::Build', is_input => 1, },
        wrapper_build_label => { is => 'Text', is_input => 1, },
    ],
    has_optional_transient => [
        _wrapped_command    => { is => 'Command::V2' },
    ],
    is_abstract => 1,
};

sub _init_subclass {
    my $subclass_name = shift;
    my $meta = $subclass_name->__meta__;
    return 1;
}

sub shortcut {
    my $self = shift;
    my $delegate = $self->_wrapped_command;
    print STDERR "calling shortcut on $self to $delegate: @_\n";
    my $result = $delegate->shortcut(@_);
    print STDERR "shortcut result: $result\n";
    return $result;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    # determine what to run
    my $command_class = $self->class;
    $command_class =~ s/::BuildStepWrapper$//;
    my %has = Genome->_wrapper_has($command_class,__PACKAGE__);
    my %params = map { 
            if ($has{$_}{is_many}) {
                $_ => [$self->$_]
            }
            else {
                my $value = $self->$_;
                $_ => $value
            }
        } keys %has;

    my $command = eval { $command_class->create(%params); };
    unless ($command) {
        my $err = $@ || '';
        $err = "Failed to create wrapped command for $class:" . $err;
        $self->delete;
        die $class->error_message($err);
    }
    $self->_wrapped_command($command);

    return $self;
}

sub execute {
    my $self = shift;
    my $build = $self->wrapper_build;
    my $label = $self->wrapper_build_label;
    unless ($build) {
        die "no build on $self!";
    }

    # determine what to run
    my $command_class = $self->class;
    $command_class =~ s/::BuildStepWrapper$//;
    my $command = $self->_wrapped_command;
    die unless $command;
    
    # determine where to link its results on the filesystem
    my $dir = $build->data_directory;
    unless (-d $dir) {
        die $self->error_message("failed to find data directory for build " . $build->__display_name__);
    }
    my $subdir = "$dir/results";
    unless (-e $subdir) {
        Genome::Sys->create_directory($subdir);
    }
    my $link = $label;
    $link =~ s/[\s\:]/_/g;
    $link = "$subdir/$link";
    $self->debug_message("Running $command_class for build.  Linking results to $link.");
    
    # run
    my $retval = $command->execute(@_);
    my $result = $command->result();
    unless ($result and $result->isa("Genome::SoftwareResult")) {
        no warnings;
        die $self->error_message("expected $command to return a Genome::SoftwareResult, got: $result");
    }

    # link to the build in the database
    my $user_link = $result->add_user(label => $label, user => $build);

    # link to the data directory
    $self->debug_message("linking to $dir/$subdir: result " . $result->__display_name__ . " from " . $result->output_dir);
    Genome::Sys->create_symlink($result->output_dir,"$link");

    return $result;
}

1;

