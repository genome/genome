package Genome::Model::Tools::BioSamtools::ParallelRefCov;

use strict;
use warnings;

use Genome;
use Workflow;

class Genome::Model::Tools::BioSamtools::ParallelRefCov {
    is => ['Workflow::Operation::Command','Genome::Model::Tools::BioSamtools'],
    workflow => sub {
        my $workflow = Workflow::Operation->create(
                                                   name => 'bio-samtools parallel ref-cov',
                                                   operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::BioSamtools::RefCov'),
                                               );
        $workflow->parallel_by('bed_file');
        return $workflow;
    },
    has => [
        regions => { is => 'Number', default_value => 10000, doc => 'The number of regions to include in each parallel process' },
        _stats_file => { is_optional => 1, },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless (defined($self->min_depth_filter)) {
        #This should look up a class variable like Genome::Model::Tools::OldRefCov::$DEFAULT_MIN_DEPTH_FILTER
        $self->min_depth_filter(1);
    }
    return $self;
}

sub pre_execute {
    my $self = shift;

    my $bed_file = $self->bed_file;

    my ($bam_basename,$bam_dirname,$bam_suffix) = File::Basename::fileparse($self->bam_file,qw/.bam/);
    my ($bed_basename,$bed_dirname,$bed_suffix) = File::Basename::fileparse($bed_file,qw/.bed/);

    $self->_stats_file($self->output_directory .'/'. $bam_basename .'_'. $bed_basename .'_STATS.tsv');

    my $regions = $self->regions;
    require Cwd;
    my $cwd = Cwd::cwd();
    chdir $self->output_directory;
    my $sub_bed_basename = $bed_basename .'_SUB_REGIONS';

    Genome::Sys->shellcmd(
                                          cmd => "split -a 4 -d -l $regions $bed_file $sub_bed_basename",
                                          input_files => [$bed_file],
                                      );
    chdir $cwd;
    my @files = glob($self->output_directory .'/'. $sub_bed_basename .'*');
    $self->bed_file(\@files);
    return 1;
}

sub post_execute {
    my $self = shift;

    my @failures = grep { $_ ne 1 } @{$self->result};
    if (@failures) {
        $self->error_message('One or more of the parallel commands failed');
        die($self->error_message);
    }

    Genome::Sys->cat(
        input_files => $self->stats_file,
        output_file => $self->_stats_file,
    );

    for my $file (@{$self->bed_file},@{$self->stats_file}) {
        unless (unlink $file) {
            $self->error_message('Failed to remove file '. $file .":  $!");
            die($self->error_message);
        }
    }
    $self->stats_file($self->_stats_file);
    return 1;
}


1;
