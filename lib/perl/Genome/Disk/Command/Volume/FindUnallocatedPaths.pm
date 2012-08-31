package Genome::Disk::Command::Volume::FindUnallocatedPaths;

use warnings;
use strict;

use Genome;
use IO::Dir;

class Genome::Disk::Command::Volume::FindUnallocatedPaths{
    is => 'Command::V2',
    has_input => [
        volume => {
            is => 'Genome::Disk::Volume',
            shell_args_position => 1,
            doc => 'Identifier for disk volume on which to find unallocated paths.',
        },
    ],
    has_optional => [
        _allocated_paths => {
            is => 'HashRef',
            is_transient => 1,
        },
        _unallocated_paths => {
            is => 'Text',
            is_many => 1,
            is_output => 1,
            is_transient => 1,
        },
    ],
    doc => 'finds unallocated paths on the provided volume',
};

sub help_detail {
    return 'Scans the provided volume for paths that are not contained in an allocation';
}

sub execute{

    my $self = shift;
    my $mount_path = $self->volume->mount_path;
    my @allocations = Genome::Disk::Allocation->get(mount_path=>$mount_path);
    unless(@allocations) {
        $self->warning_message("No allocations on $mount_path.");
        return;
    }
    my %allocated_paths;
    for my $allocation (@allocations) {
        my @parts = split(/\/+/, $allocation->absolute_path);
        @parts = grep($_, @parts);
        my $dir = \%allocated_paths;
        for my $part(@parts) {
            unless(exists $dir->{$part}) {
                $dir->{$part} = {};
            }
            $dir = $dir->{$part};
        }
    }
    $self->_allocated_paths(\%allocated_paths);
    my ($allocated_subpaths, @unallocated_paths) = $self->find_unallocated_paths($mount_path);
    print join("\n", @unallocated_paths), "\n";
    $self->_unallocated_paths(\@unallocated_paths);
    return 1;
}

sub find_unallocated_paths{

    my ($self, $path) = @_;
    my @parts = split(/\/+/, $path);
    @parts = grep($_, @parts);
    my @unallocated_children;
    my $has_allocated_children = 0;
    my $allocated_paths = $self->_allocated_paths;
    my $dir = $allocated_paths;
    for my $part(@parts) {
        unless(exists $dir->{$part}) {
            return 0, $path;
        }
        $dir = $dir->{$part};
    }
    unless(keys %{$dir}) {
        return 1;
    }

    if (-l $path){
        return 0, $path;
    }
    if (-d $path){
        my $dh = IO::Dir->new($path);
        unless($dh){
            die $self->error_message("Could not open directory handle for $path.");
        }
        while(my $subpath = $dh->read()){
            next if $subpath =~ /^\.\.?$/;
            my ($allocated_subpaths, @unallocated_subpaths) = $self->find_unallocated_paths("$path/$subpath");
            if($allocated_subpaths) {
                $has_allocated_children = 1;
            }
            @unallocated_children = (@unallocated_children, @unallocated_subpaths);
        }
        if($has_allocated_children) {
            return 1, @unallocated_children;
        }else{
            return 0, "$path/";
        }
    }else{
        return 0, $path;
    }
}

1;
