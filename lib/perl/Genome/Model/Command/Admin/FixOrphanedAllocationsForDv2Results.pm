package Genome::Model::Command::Admin::FixOrphanedAllocationsForDv2Results;

use warnings;
use strict;

use Genome;

use File::Find::Rule qw();
use List::MoreUtils qw(uniq);

class Genome::Model::Command::Admin::FixOrphanedAllocationsForDv2Results {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Builds to fix.'
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

    for my $build ($self->builds) {
        my @symlinks = find_dir_symlinks($build->data_directory);

        for my $s (@symlinks) {
            my $t = readlink $s;
            my $a = Genome::Disk::Allocation->get_allocation_for_path($t);

            next unless $a && $a->owner_class_name->isa('Genome::Model::Tools::DetectVariants2::Result::Base');

            my $lookup_hash = eval { Genome::SoftwareResult::_validate_lookup_hash((split /-/, $s)[-1]) };
            next unless $lookup_hash;

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
            print 'Purging ' . $a->absolute_path . "\n";
            unless ($self->dry_run) {
                unlink $s;
                $a->purge(reason => $self->reason($build->id));
            }
            $owner_class_name->_unlock_resource($owner_lock);
        }
        print $build->id. " done\n";
    }

    return 1;
}

sub reason {
    my ($self, $build_id) = @_;
    return sprintf("Removing orphaned allocations for build %s", $build_id);
}

1;
