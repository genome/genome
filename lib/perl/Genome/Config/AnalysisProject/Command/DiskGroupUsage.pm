package Genome::Config::AnalysisProject::Command::DiskGroupUsage;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::DiskGroupUsage {
    is => 'Command::V2',
    has => [
        analysis_projects => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to report disk usage',
            is_many             => 1,
            shell_args_position => 1,
        },
        print_allocation_ids_for_disk_group => {
            is => 'Genome::Disk::Group',
            is_optional => 1,
        },
    ],
    doc => 'report on the disk usage of analysis projects per disk group',
};

sub help_detail {
    return <<EOHELP
This command will take several minutes to run.

Running builds will not be accounted for properly as the allocations are not complete.

The output is a tab-delimited disk usage summary for each disk group used per analysis project.

The columns of the output file are:
analysis project ID, disk group name, total KB

OPTIONAL: Each allocation ID is printed out, tab indented after the disk group summary line

EOHELP
}

sub execute {
    my $self = shift;

    my @headers = qw/
        analysis_project
        disk_group_name
        total_kb
    /;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
    );
    for my $anp ($self->analysis_projects) {
        my %disk_group_usage;
        for my $model ($anp->models) {
            for my $build ($model->builds) {
                my @allocations = $build->disk_usage_allocations;
                for my $allocation (@allocations) {
                    $disk_group_usage{$allocation->disk_group_name}{$allocation->id} = $allocation->kilobytes_used;
                }
            }
        }
        for my $disk_group_name (sort keys %disk_group_usage) {
            my $total_kb = 0;
            map { $total_kb += $_ } values %{$disk_group_usage{$disk_group_name}};
            my $data = {
                'analysis_project' => $anp->id,
                'disk_group_name' => $disk_group_name,
                'total_kb' => $total_kb,
            };
            $writer->write_one($data);
            if ($self->print_allocation_ids_for_disk_group && ($self->print_allocation_ids_for_disk_group->name eq $disk_group_name) ) {
                for my $allocation_id (sort keys %{$disk_group_usage{$disk_group_name}} ) {
                    print "\t". $allocation_id ."\n";
                }
            }
        }
    }
    $writer->output->close();

    return 1;
}


1;
