package Genome::Model::Tools::Fastqc::GenerateReports;

use strict;
use warnings;
use version;

use Genome;

class Genome::Model::Tools::Fastqc::GenerateReports {
    is => 'Genome::Model::Tools::Fastqc',
    has_input => [
        input_files => {
            doc => 'The input files (fastq,sam,bam,etc.) to generate quality reports for. Comma delimited list.',
        },
        fastq_files => {
            doc => 'The input fastq files to generate quality reports for. Comma delimited list.',
            is_deprecated => 1,
            is_optional => 1,
        },
        report_directory => {
            is => 'Text',
            doc => 'Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it.  If this option is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.',
        },
        format => {
            doc => 'Bypasses the normal sequence file format detection and forces the program to use the specified format. (version 0.10.0 or greater)',
            is => 'Text',
            valid_values => ['bam','sam','bam_mapped','sam_mapped','fastq'],
            is_optional => 1,
        },
        casava => {
            doc => 'Files come from raw casava output. Files in the same sample group (differing only by the group number) will be analysed as a set rather than individually. Sequences with the filter flag set in the header will be excluded from the analysis. Files must have the same names given to them by casava (including being gzipped and ending with .gz) otherwise they won\'t be grouped together correctly. (version 0.10.0 or greater)',
            is => 'Boolean',
            default_value => 0,
            is_optional => 1,
        },
        extract => {
            doc => 'If set then the zipped output file will be uncomressed in the same directory after it has been created.  By default this option will be set if fastqc is run in non-interactive mode. (version 0.10.0 or greater)',
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
        },
        group => {
            doc => 'Disable grouping of bases for reads >50bp. All reports will show data for every base in the read.  WARNING: Using this option will cause fastqc to crash and burn if you use it on really long reads, and your plots may end up a ridiculous size. You have been warned! (version 0.10.0 or greater)',
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
        },
        threads => {
            doc => 'Specifies the number of files which can be processed simultaneously.  Each thread will be allocated 250MB of memory so you shouldn\'t run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine. (version 0.10.0 or greater)',
            is => 'Integer',
            default_value => 1,
            is_optional => 1,
        },
        contaminants => {
            doc => 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.  The file must contain sets of named contaminants in the form name[tab]sequence.  Lines prefixed with a hash will be ignored. (version 0.10.0 or greater)',
            is => 'Text',
            is_optional => 1,
        },
        quiet => {
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
            doc => 'Supress all progress messages on stdout and only report errors. (version 0.10.0 or greater)',
        },
    ],
};

sub help_brief {
    "FastQC - A high throughput sequence QC analysis tool";
}

sub help_detail {
    "FastQC reads a set of sequence files and produces from each one a quality control report consisting of a number of different modules, each one of which will help to identify a different potential type of problem in your data. Any bugs in fastqc should be reported either to simon.andrews\@babraham.ac.uk or in www.bioinformatics.bbsrc.ac.uk/bugzilla/";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    if ($self->fastq_files) {
        $self->warning_message('The option fastq_files is deprecated.  Please use input_files since FastQC now takes SAM/BAM format files.');
        $self->input_files($self->fastq_files);
    }
    return $self;
}

sub execute {
    my $self = shift;

    my @input_files;
    if (ref($self->input_files) eq 'ARRAY') {
        @input_files = @{$self->input_files};
    } else {
        @input_files = split(',',$self->input_files);
    }
    if ( version->parse($self->use_version) < version->parse('0.10.0') ) {
        my $cmd = $self->fastqc_path .' -Djava.awt.headless=true -Dfastqc.output_dir='. $self->report_directory .' uk.ac.bbsrc.babraham.FastQC.FastQCApplication '. join(' ',@input_files) ;
        $self->run_java_vm(
            cmd => $cmd,
            input_files => \@input_files,
        );
    } else {
        my $cmd = $self->fastqc_path .'/fastqc  --outdir '. $self->report_directory .' --threads '. $self->threads;
        my $files_string = join(' ',@input_files);
        if ($self->quiet) {
            $cmd .= ' --quiet';
        }
        if ($self->contaminants) {
            push @input_files, $self->contaminants;
            $cmd .= ' --contaminants '. $self->contaminants;
        }
        unless ($self->group) {
            $cmd .= ' --nogroup';
        }
        if ($self->extract) {
            $cmd .= ' --extract';
        } else {
            $cmd .= ' --noextract';
        }
        if ($self->casava) {
            $cmd .= ' --casava';
        }
        if ($self->format) {
            $cmd .= ' --format '. $self->format;
        }
        $cmd .= ' '. $files_string;
        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => \@input_files,
            set_pipefail => 1,
        );
    }
    return 1;
}
