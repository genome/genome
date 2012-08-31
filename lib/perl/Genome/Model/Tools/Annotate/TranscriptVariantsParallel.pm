package Genome::Model::Tools::Annotate::TranscriptVariantsParallel;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Temp qw/ tempdir /;
use Genome::Utility::IO::SeparatedValueReader;
use Workflow;
use File::Basename;
use Cwd 'abs_path';
use Sys::Hostname;
use DateTime;

class Genome::Model::Tools::Annotate::TranscriptVariantsParallel{
    is => ['Workflow::Operation::Command','Genome::Model::Tools::Annotate::TranscriptVariants'],
    workflow => sub {
        my $workflow = Workflow::Operation->create(
            name => 'parallel transcript variants',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Annotate::TranscriptVariants')
        );
        $workflow->parallel_by('variant_file');
        return $workflow;
    },
    has => [
        split_by_chromosome => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'enable this flag to split the variants file by chromosome (default behavior)',
        },
        split_by_number => {
            is => 'Number',
            is_optional => 1,
            doc => 'enable this flag to split the variant file by number',
        },
        log_directory => {
            is => 'Path',
            is_optional => 1,
            doc => 'If provided, logs for each transcript variants subprocess are stored in the directory'
        },
        _output_file => {
            is => 'Text',
            is_optional => 1,
        },
        _variant_files => {
            is => 'File::Temp',
            is_many => 1,
            is_optional => 1,
        },
        _temp_dir => {
            is => 'Text',
            is_optional => 1,
        },
    ], 
};

sub pre_execute {
    my $self = shift;

    if ($self->output_file eq 'STDOUT') {
        $self->error_message("Must specify an output file") and die;
    }

    # Simple checks on command line args
    if (-s $self->output_file) {
        $self->error_message($self->output_file . " exists and has a size, exiting") and die;
    }
    $self->_output_file($self->output_file);

    unless (-s $self->variant_file) {
        $self->error_message($self->variant_file . " does not exist or has no size, exiting") and die;
    }

    # Useful information for debugging...
    my $dt = DateTime->now;
    $dt->set_time_zone('America/Chicago');
    my $date = $dt->ymd;
    my $time = $dt->hms;
    my $host = hostname;
    $self->status_message("Executing on host $host on $date at $time");

    # Make log directory and have workflow put child process output there
    if (defined $self->log_directory) {
        unless (-d $self->log_directory) {
            my $log_dir = Genome::Sys->create_directory($self->log_directory);
            unless (defined $log_dir) {
                $self->error_message("Error creating log directory at " . $self->log_directory);
                return;
            }
        }
        $self->_operation->log_dir($self->log_directory);
    }


    # Make temp dir
    my $temp_dir =  tempdir(
        'annotation_temp_XXXXX', 
        DIR => abs_path(dirname($self->output_file)), 
        CLEANUP => 1
    );
    unless (-d $temp_dir) {
        $self->error_message("Could not create temporary output directory at " . 
            abs_path(dirname($self->output_file)));
        die;
    }
    chmod(0775, $temp_dir);
    $self->_temp_dir($temp_dir);
    
    my @splitFiles;
    # Split by line number
    if ($self->split_by_number and not $self->split_by_chromosome) {
        $self->status_message("Splitting variant file " . $self->variant_file . " into " .
            $self->split_by_number . " line chunks");

        my $inputFileHandler = IO::File->new($self->variant_file, "r");
        
        my $done = 0;
        until ($done == 1) {
            my $temp = File::Temp->new(DIR => $self->_temp_dir);
            push @splitFiles, $temp;
            for (1..$self->split_by_number) {
                my $line = $inputFileHandler->getline;
                unless ($line) {
                    $done = 1;
                    last;
                }
                $temp->print($line);
            }
            $temp->close;
            $self->variant_file([map { $_->filename } @splitFiles]);
        }
        $inputFileHandler->close;
    }

    # Split file by chromosome
    else {
        $self->status_message("Splitting variant file " . $self->variant_file . " by chromosome");

        my @variant_columns = $self->variant_attributes;
        push @variant_columns, split(/,/, $self->extra_columns) if $self->extra_columns;
        my $reader = Genome::Utility::IO::SeparatedValueReader->new (
            input => $self->variant_file,
            headers => \@variant_columns,
            separator => '\t',
            is_regex => 1,
            ignore_extra_columns => 1,
        );
        unless ($reader) {
            $self->error_message("Could not get reader for " . $self->variant_file);
            die;
        }

        my $currChrom = '';
        my $fh;
        while (my $line = $reader->next) {
            my $chrom = $line->{chromosome_name};
            if ($chrom ne $currChrom and ($chrom !~ /[M|N]T/ or $currChrom !~ /[M|N]T/)) {
                $self->status_message("Making new file for chromosome $chrom");
                $currChrom = $chrom;
                $fh->close if $fh;
                $fh = File::Temp->new (DIR => $self->_temp_dir);
                chmod(0664, $fh->filename);
                push @splitFiles, $fh;
            }

            my @newline;
            foreach (@variant_columns) {
                push @newline, $line->{$_};
            }

            $splitFiles[-1]->print(join("\t", @newline)."\n");
        }
        $fh->close if $fh;
        $self->variant_file([map { $_->filename } @splitFiles]);
    }
    
    $self->status_message("Variant file split into " . scalar @splitFiles . " files");

    $self->_is_parallel(1);
    $self->no_headers(1);
    $self->_variant_files(\@splitFiles);
    return 1;
}

sub post_execute {
    my $self = shift;

    foreach my $error (@Workflow::Simple::ERROR) {
        print $error->error;
    }

    my @output_files;
    my @zero_size_files;
    for (@{$self->output_file}){
        print "Output file undefined\n" and next unless $_;
        print "Output file $_ doesn't exist\n" and next unless -e $_;

        unless (-s $_) {
            print "Output file $_ has no size\n";
            push @zero_size_files, $_;
        }
        else {
            push @output_files, $_;
        }
    }

    unless(@output_files){
        $self->error_message("No non-zero sized output files, exiting") and die;
    }
    Genome::Sys->cat(
        input_files => \@output_files, 
        output_file => $self->_output_file,
    );

    for my $file (@{$self->variant_file},@output_files, @zero_size_files) {
        unless (unlink $file) {
            $self->error_message('Failed to remove file '. $file .":  $!");
        }
    }

    unless (rmdir ($self->_temp_dir)) {
        $self->warning_message("Could not remove temporary annotation directory at " . $self->_temp_dir . "\n");
    }

    $self->output_file($self->_output_file);
    return 1;
}
1;

