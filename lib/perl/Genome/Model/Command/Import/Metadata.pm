#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use IO::File;
use Genome;

package Genome::Model::Command::Import::Metadata;

class Genome::Model::Command::Import::Metadata {
    is => 'Command::V2',
    has_input => [
        input_path => {
            is => 'FilesystemPath',
            default_value => '-',
            shell_args_position => 1,
            doc => 'the source of serialized model data (use "-" for standard input)',
        },
        log_path => {
            is => 'FilesystemPath',
            default_value => '-',
            doc => 'the path at which to place the log file for this import',
        },
    ],
    has_param => [
        update => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'if objects already exist, update them as needed instead of failing',
        },
        ignore_differences => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'if objects already exist and are different, ignore those differences instead of failing',
        },
    ],
    doc => 'import serialized Genome Models, typically from another system instance'
};

sub help_synopsis {
    return <<EOS
# on one system
genome model export metadata "name like 'myproject%'" >myfile.dat

# on another
genome model import metadata myfile.dat
EOS
}

sub execute {
    my $self = shift;
    my $f = $self->input_path;

    my $log_fh;
    if ($self->log_path eq '-') {
        $log_fh = 'STDOUT';
    }
    else {
        $log_fh = Genome::Sys->open_file_for_writing($self->log_path);
    }

    die "No database dump file specified!" unless $f;
    die "Failed to find $f!" unless -e $f or $f eq '-';

    my @dumps = grep { $_ !~ /^#/ } IO::File->new($f)->getlines;

    my %loaded;
    for my $dump (@dumps) {
        my $hash = do {
            no strict;
            no warnings;
            eval $dump;
        };
        die "NO HASH?: $dump" unless $hash;

        my $class = ref($hash);
        my $id = $hash->{id};

        eval "use $class";

        next if $class->isa("UR::Value");

        my $dbc = delete $hash->{db_committed};
        
        if ($dbc) {
            for my $key (%$dbc, 'id') {
                no warnings;
                unless ($key eq 'id' or $hash->{$key} eq $dbc->{$key}) {
                    die "data discrepancy: $hash with ID $hash->{id} has key $key with value $hash->{$key} but db_committed is $dbc->{$key}\n";
                }
                
                if ($hash->{$key} =~ m|gscmnt|) {
                    if ($hash->{$key} =~ s|gscmnt|opt/gms/fs|) {
                        $log_fh->print("updated $class $id $key to $hash->{$key}\n");
                    }
                    else {
                        die "error updating $class $id $key gscmnt content!";
                    }
                }
            }
        }

        $loaded{$class}{$id} = $hash;
    }

    for my $class (sort keys %loaded) {
        my $d = $loaded{$class};
        my @ids = sort keys %$d;

        # pre-cache queries below for speed
        my @objs = $class->get(id => \@ids);

        # get or create each
        for my $id (@ids) {
            # $id = "888.99" if $class->isa("Genome::Db");
            my $hash = $loaded{$class}{$id};
            my $prev = $class->get($id);
            if ($prev) { 
                $log_fh->print("## FOUND $class $id: " . $prev->__display_name__ . "\n");
            }
            else {
                if ($class->isa("Genome::Db")) {
                    # these exist because of filesystem data being in place
                    # the installation is slow, so just warn the user to install them
                    my ($source) = ($class =~ /Genome::Db::([^\:]+)/);
                    
                    # if we can speed this up, run the intaller below
                    # for now just tell the user what they must run
                    $source = lc($source);
                    $self->warning_message("EXTERNAL DATABASE NOT INSTALLED.  DO: genome db $source install $id");
                    
                    #my $installer_class = "Genome::Db::${source}::Command::Install";
                    #unless (UR::Object::Type->get($installer_class)) {
                    #    $self->warning_message("No installer $installer_class for $class!  Install manually!");
                    #    next;
                    #}
                    #eval {
                    #    $installer_class->execute(version => $id);
                    #};
                    #if ($@) {
                    #    $self->warning_message("errors installing $source $id: $@");
                    #}
                    #$prev = $class->get($id);
                    #unless ($prev) {
                    #    $self->warning_message("Failed to find $class $id!  Install manually.");
                    #}
                }
                else {
                    $log_fh->print("## IMPORTING $class $id: " . UR::Util::d($hash) . "\n");
                    my $entity = UR::Context->_construct_object($class,%$hash, id => $id);
                    die "failed create for $class $id\n" unless $entity;
                    $entity->__signal_change__('create');
                    $entity->{'__get_serial'} = $UR::Context::GET_COUNTER++;
                    $UR::Context::all_objects_cache_size++;
                    #print ">>> $hash $entity\n";
                }
            }
        }
    }

    my $err = sub { };
    for my $c (keys %loaded) {
        do {
            no strict;
            no warnings;
            my $m = $c . '::__errors__';
            *$m = $err;
        };

        my $h = $loaded{$c};
        for my $i (keys %$h) {
            my $o = $c->get($i);
            unless ($o){
                print "\nNo $c $i!";
                next;
            }
            #die "No $c $i!" unless $o;
            #print "$c $i $o\n";
            
            # TODO: standardize on a method in the class to initialize imported data
            # This shoudl match a companion method to export data.
            if ($o->isa("Genome::Disk::Volume")) {
                # this only changes when we do something special to mount a remote system r/w
                # and would require interaction with the remote GMS
                $o->can_allocate(0);
            }
        }
    }
    
    $self->status_message("import complete");
    return 1;
};

1;

