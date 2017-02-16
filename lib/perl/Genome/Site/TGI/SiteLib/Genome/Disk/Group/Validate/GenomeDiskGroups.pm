package Genome::Disk::Group::Validate::GenomeDiskGroups;

use strict;
use warnings;

use Carp qw(croak);
use List::MoreUtils qw(any uniq);

use Try::Tiny qw(try);

sub validate {
    my $class = shift;
    my $disk_group = shift;

    unless (is_genome_disk_group($disk_group)) {
        croak is_genome_disk_group_error($disk_group);
    }
    return;
};

sub is_genome_disk_group_error {
    my $self = shift;
    return sprintf("Disk group name (%s) not allowed.",
        $self->disk_group_name);
}

sub is_genome_disk_group {
    my $self = shift;
    my $name = $self->disk_group_name;
    return 1 if $ENV{UR_DBI_NO_COMMIT};
    return any { $_ eq $name } genome_disk_group_names();
}

sub is_cle_disk_group_name {
    my $name = shift;

    return 1 if $ENV{UR_DBI_NO_COMMIT};
    return any { $_ eq $name } cle_disk_group_names();
}

sub cle_disk_group_names {
    return (
        'cle_alignments',
        'cle_genome_models',
    );
}


sub genome_disk_group_names {
    return uniq cle_disk_group_names(), (
        # hard-coded for now because config is in a mess
        'cle_research_alignments',
        'cle_research_genome_models',
        'cle_research_analysis',
        'info_genome_models',
        'info_alignments',
        'apipe_ci',
        grep { defined $_ }
            map { try { Genome::Config::get($_) } }
            (
                qw(disk_group_archive
                disk_group_dev
                disk_group_references
                disk_group_alignments
                disk_group_models
                disk_group_research)
            )
    );
}

1;
