package Genome::Config::AnalysisProject::Command::DiskUsage;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::DiskUsage {
    is => 'Command::V2',
    has => [
        analysis_projects => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to report disk usage',
            is_many             => 1,
            shell_args_position => 1,
        },
    ],
    doc => 'report on the disk usage of analysis projects',
};

sub help_detail {
    return <<EOHELP
This command will take several minutes to run.

Running builds will not be accounted for properly as the allocations are not complete.

The output is a tab-delimited disk usage summary for each model type and configuration item per analysis project.

The sum of the model type disk usage is not necessarily the total disk usage for an Analysis Project.  Results can be shared between multiple model types.  
For instance, alignment files, BAMs, can be shared between SingleSampleGenotype, SomaticValidation and/or ReferenceAlignment depending on the processign profiles used.
The same goes for SomaticVariation and SomaticValidation results. VCF and variant files can be shared between the two model types.

The colums of the output file are:
analysis project ID, configuration item ID, subclass name, model count, build count, total KB, total base pair, bytes per base

EOHELP
}

sub execute {
    my $self = shift;

    my @lines;
    for my $anp ($self->analysis_projects) {
        for my $config_item ($anp->config_items) {
            my %model_type_disk_usage;
            for my $model ($config_item->models) {
                my $subclass_name = $model->subclass_name;
                $model_type_disk_usage{$subclass_name}{'model_count'}++;
                for my $build ($model->builds) {
                    $model_type_disk_usage{$subclass_name}{'build_count'}++;
                    my @allocations = $build->disk_usage_allocations;
                    for my $allocation (@allocations) {
                        $model_type_disk_usage{$subclass_name}{'allocations'}{$allocation->id} = $allocation->kilobytes_used;
                    }
                    my @instrument_data = $build->instrument_data;
                    for my $instrument_data (@instrument_data) {
                        $model_type_disk_usage{$subclass_name}{'instrument_data'}{$instrument_data->id} = $instrument_data->total_bases_read;
                    }
                }
            }
            for my $model_type (sort keys %model_type_disk_usage) {
                my $total_kb = 0;
                map { $total_kb += $_ } values(%{$model_type_disk_usage{$model_type}{'allocations'}});

                my $total_bp = 0;
                map { $total_bp += $_ } values (%{$model_type_disk_usage{$model_type}{'instrument_data'}});

                my $bytes_per_base = 'na';
                if ($total_bp) {
                    $bytes_per_base = ($total_kb * 1024) / $total_bp;
                }
                my $data = {
                    'analysis_project' => $anp->id,
                    'config_item' => $config_item->id,
                    'subclass_name' => $model_type,
                    'model_count' => $model_type_disk_usage{$model_type}{'model_count'} || 0,
                    'build_count' => $model_type_disk_usage{$model_type}{'build_count'} || 0,
                    'total_kb' => $total_kb,
                    'total_bp' => $total_bp,
                    'bytes_per_bp' => $bytes_per_base, 
                };  
                push @lines, $data; 
            }
        }
    }
    my @headers = qw/
        analysis_project
        config_item
        subclass_name
        model_count
        build_count
        total_kb
        total_bp
        bytes_per_bp
    /;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,    
    );
    for my $data (@lines) {
        $writer->write_one($data);
    }        
    return 1;
}


1;
