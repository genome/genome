package Genome::Model::Command::Admin::FixOrhpanedAllocationsForDv2Results;

use warnings;
use strict;

use Genome;

use File::Find::Rule qw();
use List::MoreUtils qw(uniq);

class Genome::Model::Command::Admin::FixOrhpanedAllocationsForDv2Results {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'Build to fix.'
        },
        dry_run => {
            is => 'Boolean',
            default_value => 0,
        },
    ],
};

sub help_brief {
    return 'delete orphaned allocations and links for DV2 results'
}

sub help_detail {
    return 'Looks for orphaned allocations in the build directory (following symlinks) and then deletes them if it is believed to be safe.';
}

sub find_dir_symlinks {
    my $in = shift;

    my @l = grep { -d readlink $_ } File::Find::Rule->symlink()->in($in);
    return uniq @l, map { find_dir_symlinks(readlink $_) } @l;
}

sub execute {
    my $self = shift;

    my @symlinks = find_dir_symlinks($self->build->data_directory);

    for my $s (@symlinks) {
        my $t = readlink $s;
        my @d = File::Spec->splitdir($t);
        my $p = File::Spec->catdir(@d[4..$#d]);
        my $a = Genome::Disk::Allocation->get(allocation_path => $p);

        next unless $a && $a->owner_class_name->isa('Genome::Model::Tools::DetectVariants2::Result::Base');

        my $lookup_hash = (split /-/, $s)[-1];
        my $owner_lock = $a->owner_class_name->_lock($lookup_hash, undef);

        # the owner is generating if it's locked?
        next unless $owner_lock;

        if ($a->owner) {
            $a->owner_class_name->_unlock_resource($owner_lock);
            next;
        }

        my $owner_class_name = $a->owner_class_name;
        print join(' ', $a->id, $owner_class_name, $a->owner_id, $a->owner ? 'T' : 'F', $lookup_hash), "\n";
        print "Unlinking $s\n";
        print 'Deleting ' . $a->absolute_path . "\n";
        unless ($self->dry_run) {
            unlink $s;
            $a->delete;
        }
        $owner_class_name->_unlock_resource($owner_lock);
    }

    return 1;
}

1;
