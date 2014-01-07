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
        verbose => {
            is => 'Boolean',
            default_value => 0,
            is_optional => 1,
            doc => 'verbose output',
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
                        if ($self->verbose){
                            $log_fh->print("updated $class $id $key to $hash->{$key}\n");
                        }
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
                if ($self->verbose){
                    $log_fh->print("## FOUND $class $id: " . $prev->__display_name__ . "\n");
                }
            }
            else {
                if ($class->isa("Genome::Db")) {
                    my ($source) = ($class =~ /Genome::Db::([^\:]+)/);
                    my $source_name = $source;
                    $source = lc($source);
                    my $installer_class = "Genome::Db::${source_name}::Command::Install";
                    unless (UR::Object::Type->get($installer_class)) {
                        die $self->error_message("No installer $installer_class for $class!");
                    }
                    $self->warning_message("EXTERNAL DATABASE NOT INSTALLED. Attempting to install source: $source, id: $id, of class: $class");
                    if ($id =~ /cosmic\/(\d+)\.(\d+)/){
                        #Set up import for 'cosmic' sources with IDs like: ('cosmic/61.1')
                        #create commands that look like:
                        #genome db cosmic install --branch=61_v1
                        my $branch = $1 . "_v" . $2;

                        eval {
                            $installer_class->execute(branch => $branch);
                        };
                        if ($@) {
                            die $self->error_message("errors installing $source $id: $@");
                        }
                    }elsif ($id =~ /^tgi\/(.*)\/(.*)\/(.*)\-(\d+)\.(\d+)/){
                        #Set up import for 'tgi' sources with IDs like: ('tgi/cancer-annotation/human/build37-20130401.1' or 'tgi/misc-annotation/human/build37-20130113.1')
                        #create commands that look like:
                        #genome db tgi install --subsource=cancer-annotation --species=human --branch=human-build37-20130401
                        #genome db tgi install --subsource=misc-annotation --species=human --branch=human-build37-20130113
                        my $subsource = $1;
                        my $species = $2;
                        my $build = $3;
                        my $date = $4;
                        my $version = $5;
                        my $branch = $species . "-" . $build . "-" . $date;
                        eval {
                            $installer_class->execute(subsource => $subsource, species => $species, branch => $branch);
                        };
                        if ($@) {
                            die $self->error_message("errors installing $source $id: $@");
                        }
                    }else{
                        die "could not parse Genome::Db source: $source";
                    }
                    $prev = Genome::Db->get($id);
                    unless ($prev) {
                        die $self->error_message("Failed to find $class $id!");
                    }
                    
                }
                else {
                    if ($self->verbose){
                        $log_fh->print("## IMPORTING $class $id: " . UR::Util::d($hash) . "\n");
                    }
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
                print "\nNo $c $i!" if $self->verbose;
                next;
            }
            
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

