#!/usr/bin/env genome-perl

package Genome::InstrumentData::Command::MoveAlignmentDirectory;
use warnings;
use strict;
use File::Path 'rmtree';
use Genome;

class Genome::InstrumentData::Command::MoveAlignmentDirectory {
    is => 'Command',
    has_input => [
        from => {
            shell_args_position => 1,
            doc => 'all or part of the alignment allocation_path'
        },
        to => {
            shell_args_position => 2,
            doc => 'the disk to which the data should move'
        },
    ],
    doc => 'move alignment data'
};

sub help_synopsis {
    return <<EOS
 genome instrument-data move-alignment-directory maq0_7_1/NCBI-human-build36/abc12/1_12345 /gscmnt/sata835
EOS
}

sub execute {
    my $self = shift;
    my $new_location = $self->to;
  
    $self->error_message("permission denied");
    return;

    my $gda;
    my $from = $self->from;
    if ($from !~ /\D/) {
        $gda = Genome::Disk::Allocation->get($from);
    }
    unless ($gda) {
        my @gda = Genome::Disk::Allocation->get(allocation_path => $from);
        if (@gda != 1) {
            @gda = Genome::Disk::Allocation->get(
                'allocation_path like' => '%' . $from . '%'
            )
        }
        if (@gda > 1) {
            my @paths = map { $_->absolute_path } @gda;
            my $msg = "Mulitple paths matching $from!:\n" . join("\n",@paths) . "\n";
            $self->error_message($msg);
            return;
        }
        elsif (@gda == 0) {
            $self->error_message("No allocations found matching $from!");
            return;
        }
    }

    my $path = $gda->allocation_path;

    my $new_alloc_cmd = Genome::Disk::Allocation::Command::Allocate->create(allocation_path=>$gda->allocation_path,
                            disk_group_name=>$ENV{GENOME_DISK_GROUP_DEV},
                            kilobytes_requested=>$gda->kilobytes_requested * 2,
                            owner_class_name=>$gda->owner_class_name,
                            owner_id=>$gda->owner_id,
                            mount_path=>$new_location);

    unless ($new_alloc_cmd->execute()) {
            die "couldn't allocate new space";
    }

    my $new_gda = Genome::Disk::Allocation->get(allocation_path=>$gda->allocation_path, mount_path=>$new_location);

    my $old_path = $gda->absolute_path;
    my $new_path = $new_gda->absolute_path;

    my $cp_cmd = "rsync --verbose --links --size-only --perms --group --times -r $old_path $new_path";
    #my $cp_cmd = sprintf("cp -r %s/* %s", $old_path, $new_path);

    print "$cp_cmd\n";
    system($cp_cmd);

    my $original_sum_cmd = sprintf('md5sum %s/* | awk \'{print $1}\'',$old_path);
    my $new_sum_cmd = sprintf('md5sum %s/* | awk \'{print $1}\'',$new_path);

    my $orig_sum = `$original_sum_cmd`;
    my $new_sum = `$new_sum_cmd`;

    if ($orig_sum ne $new_sum) {
        die "MD5Sums don't match!! Aborting before I delete the old data and cause lots of trouble.\n$orig_sum\n\n$new_sum\n\n"
    } else {
        print "Sums look good... proceeding with deleting the original allocation\n$orig_sum\n\n$new_sum\n";
    }

    #$gda->deallocate() || die 'deallocate failed?';
    #
    #print "Deleting $old_path\n";
    #rmtree($old_path);
}

1;

