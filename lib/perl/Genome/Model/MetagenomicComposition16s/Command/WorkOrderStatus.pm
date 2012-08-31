package Genome::Model::MetagenomicComposition16s::Command::WorkOrderStatus; 

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require Carp;

class Genome::Model::MetagenomicComposition16s::Command::WorkOrderStatus {
    is => [qw/ Genome::Command::Base Genome::Report::GeneratorCommand /],
    has => [
        work_order => {
            is => 'Genome::WorkOrder',
            is_many => 0,
            shell_args_position => 1,
            doc => 'Work order id or name.',
        },
    ],
};

sub help_brief { 
    return 'Get the status for models in a work order';
}

sub help_detail {
    return help_brief();
}

sub execute {
    my $self = shift;
    
    my $work_order = $self->work_order;
    my @items = $work_order->items;
    if ( not @items ) {
        $self->error_message('No items for work order: '.$work_order->name);
        return;
    }

    my %stats;
    my @headers = (qw/ sample-name model-id sequencing_platform build-id status oriented-fastas /);
    my @data;
    for my $item ( @items ) {
        my @datum;
        push @data, \@datum;
        $stats{items}++;
        my $sample = $item->sample;
        if ( not $sample ) {
            print STDOUT 'NO SAMPLE FOR '.$item->id."\n";
            next;
        }
        push @datum, $sample->name;
        $stats{sample}++;

        my ($model) = grep { $_->class =~ /MetagenomicComposition16s|AmpliconAssembly/ } $sample->models;
        if ( not $model ) {
            push @datum, "NO_MODEL";
            next;
        }
        push @datum, $model->id;
        push @datum, $model->sequencing_platform;
        $stats{model}++;

        my @builds = $model->builds;
        if ( not @builds ) {
            push @datum, 'NO_BUILD';
            next;
        }
        push @datum, $builds[$#builds]->id, $builds[$#builds]->status;
        $stats{build}++;
        $stats{ lc($builds[$#builds]->status) }++;

        if ( $builds[$#builds]->status eq 'Succeeded' ) {
            push @datum, $builds[$#builds]->oriented_fasta_file;
        }
    }

    # Generate report
    my $report = $self->_generate_report_and_execute_functions(
        name => 'Work Order Status',
        description => 'Work Order Status for '.$work_order->name,
        row_name => 'sample',
        headers => \@headers,
        rows => \@data,
    ) or return;

    print STDERR Dumper(\%stats);

    return $report;
}

1;

