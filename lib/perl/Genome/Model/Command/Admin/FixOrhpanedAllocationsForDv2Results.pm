package Genome::Model::Command::Admin::FixOrhpanedAllocationsForDv2Results;

use warnings;
use strict;

use Genome;

require File::Find;
use File::Grep qw( fgrep fmap );

class Genome::Model::Command::Admin::FixOrhpanedAllocationsForDv2Results {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'Build to fix.'
        },
    ],
};

sub help_brief {
    return 'Delete orphaned allocations and links for dv results'
}

sub help_detail {
    return 'This command will parse the error files of a build to find and delete orphaned allocations for detect variants results. It then goes through the variant sub directory and recursively finds and removes broken symlinks. If no allocations are found, broken symlinks are removed.';
}

sub execute {
    my $self = shift;


    my $delete_orphaned_allocations = $self->_parse_error_logs_to_find_and_delete_orphaned_allocations;
    return if not $delete_orphaned_allocations;

    my $remove_broken_links = $self->_remove_broken_links_in_variants_directory;
    return if not $remove_broken_links;

    return 1;
}

sub _parse_error_logs_to_find_and_delete_orphaned_allocations {
    my $self = shift;

    $self->status_message('Get error logs...');
    $self->status_message('Log dir: '.$self->build->data_directory);
    my $log_directory = $self->build->log_directory;
    my @err_logs;
    File::Find::find(
        {
            wanted => sub { 
                push @err_logs, $_ if /\.err$/
            },
            no_chdir => 1,
        },
        $log_directory,
    );
    $self->status_message('Found '.@err_logs.' error logs');
    return 1 if not @err_logs; # ok

    $self->status_message('Find orphaned allocations in error logs...');
    my @owner_ids;
    for my $err_log ( @err_logs ) {
        push @owner_ids, grep { defined } fmap { /no software result for it's owner ID \((\d+)\)\./; $1; } $err_log;
    }
    $self->status_message('Found '.@owner_ids.' owner ids for orphaned allocations');
    return 1 if not @owner_ids;
    my @allocations = Genome::Disk::Allocation->get(owner_id => \@owner_ids);
    $self->status_message('Found '.@allocations.' possible orphaned allocations');

    $self->status_message('Delete existing orphaned allocations...');
    my $deleted = 0;
    for my $allocation ( @allocations ) {
        my $owner = $allocation->owner;
        next if $owner; # do not delete allcoations with owners
        $self->status_message('Owner: '.join(' ', map { $allocation->$_ } (qw/ id class /)) );
        $allocation->delete;
        $deleted++;
    }
    $self->status_message("Deleted $deleted orphaned allocations");

    return 1;
}

sub _remove_broken_links_in_variants_directory {
    my $self = shift;

    $self->status_message('Find broken links...');
    my $variants_directory = $self->build->data_directory.'/variants';
    $self->status_message('Variants dir: '.$variants_directory);
    my @broken_links;
    File::Find::find(
        {
            wanted => sub{ 
                return if not -l $_;
                return if -e $_;
                push @broken_links, $_;
            },
            no_chdir => 1,
            follow_fast => 1,
        },
        $variants_directory,
    );
    $self->status_message('Found '.@broken_links.' broken links');

    $self->status_message('Remove broken links...');
    for my $broken_link ( @broken_links ) {
        unlink $broken_link;
    }
    $self->status_message('Remove broken links...OK');

    return 1;
}

1;

