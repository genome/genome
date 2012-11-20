package Genome::Model::Tools::Annovar::Db;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

class Genome::Model::Tools::Annovar::Db {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        table_name => {
            is => 'String',
        },
        buildver => {
            is => 'String',
            valid_values => ["hg18", "hg19"],
            default_value => "hg19",
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);;
    return unless ($self);
    
    $self->_prepare_staging_directory;

    my $cmd = "Annovar.d/annotate_variation.pl --downdb --buildver ".$self->buildver." ".$self->table_name." ".$self->temp_staging_directory;
print "Running command $cmd\n";
    $self->status_message("Running command $cmd");
    my $rv = Genome::Sys->shellcmd(cmd => $cmd);
    
    unless ($rv) {
        $self->error_message("Couldn't download annovar db ".$self->table_name." for buildver ".$self->buildver);
        return $rv;
    }

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;
    return $self;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("annovardb-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    my $directory = join("/", "build_merged_alignments", $self->id, $base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

1;

