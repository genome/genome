package Genome::Model::Tools::Sam::Debroadify;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::Sam::Debroadify {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        input_file => {
            is  => 'Text',
            doc => 'Full path to input SAM/BAM.',
        },
        output_file => {
            is => 'Text',
            doc => 'Full path to output SAM/BAM.',
        },
        reference_file => {
            is => 'Text',
            doc => 'Full path to the reference FASTA to build the output BAM/SAM against.',
            is_optional => 1,
        },
    ],
};

sub validate_inputs {
    my $self = shift;
    if ($self->output_file =~ /\.bam$/i && ! defined $self->reference_file) {
        die 'In order to output a BAM file, a reference_file must be provided.';
    }
    if (defined $self->reference_file && ! -e $self->reference_file) {
        die 'Could not find reference file at: ' . $self->reference_file;
    }
    unless (-e $self->input_file) {
        die 'Could not locate input file located at: ' . $self->input_file;
    }
    if ($self->output_file =~ /\.sam$/i && $self->reference_file) {
        my $message = join(' ',
            'ReorderSam seem to not work well with SAM output.',
            'Which means the resulting SAM will either have a seqdict that is out of order',
            'or if we bypass ReorderSam or will have information mistakenly stripped.',
            'Aborting.'
        );
        die $message;
    }
    return 1;
}

sub samtools_output_options {
    my $self = shift;
    my @samtools_output_options;
    if ($self->output_file =~ /\.bam$/i) {
        push @samtools_output_options, '-b';
    }
    if ($self->reference_file) {
        push @samtools_output_options, '-T ' . $self->reference_file;
    }
    return join(' ', @samtools_output_options);
}

sub samtools_input_options {
    my $self = shift;
    my @samtools_input_options;
    if ($self->input_file =~ /\.sam$/i) {
        push @samtools_input_options, '-S';
    }
    return join(' ', @samtools_input_options);
}

sub execute {
    my $self = shift;

    $self->validate_inputs();

    my $samtools_path = $self->path_for_samtools_version();
    my $samtools_output_options = $self->samtools_output_options();
    # If we have a reference_file we will be reordering it so save to tmp.
    my $chr_rename_output_file = ($self->reference_file ? Genome::Sys->create_temp_file_path() : $self->output_file);
    my $chr_rename_output = IO::File->new("| $samtools_path view -S $samtools_output_options -o $chr_rename_output_file -");
    unless ($chr_rename_output) {
        die 'Could not open the (temp) output pipe.';
    }

    my $samtools_input_options = $self->samtools_input_options();
    my $input = IO::File->new("$samtools_path view -h $samtools_input_options " . $self->input_file . ' |');
    unless ($input) {
        die 'Could not open the input pipe.';
    }
    $self->debug_message('Converting Broad chromosome references to TGI style.');
    $self->debug_message('output will be at: ' . $chr_rename_output_file);
    while (my $line = $input->getline) {
        if ($line =~ /^\@SQ/ && $line =~ /chr/) {
            $line =~ s/\@SQ\tSN:chrM/\@SQ\tSN:MT/;
            $line =~ s/\@SQ\tSN:chr(X|Y)/\@SQ\tSN:$1/;
            $line =~ s/\@SQ\tSN:chr(.*)\w/\@SQ\tSN:$1/;
        }
        if ($line !~ /^\@/ && $line =~ /chr/) {
            my @sam_fields = split("\t", $line);
            $sam_fields[2] =~ s/^chrM$/MT/;
            $sam_fields[2] =~ s/^chr(.*)(_|\w)/$1$2/;
            $sam_fields[6] =~ s/^chrM$/MT/;
            $sam_fields[6] =~ s/^chr(.*)(_|\w)/$1$2/;
            $line = join("\t", @sam_fields);
        }
        print $chr_rename_output $line;
    }
    close $input;

    close $chr_rename_output;
    unless (-s $chr_rename_output_file) {
        die 'output has no size: ' . $chr_rename_output_file;
    }

    # This needs to be done because Broad's reference is sorted 1-22,X,Y but ours is 1-9,X,Y,10-22.
    if ($self->reference_file) {
        $self->debug_message('Reordering BAM to match reference.');
        $self->debug_message('output will be at: ' . $self->output_file);
        my $reorder_cmd = Genome::Model::Tools::Picard::ReorderSam->create(
            input_file => $chr_rename_output_file,
            output_file => $self->output_file,
            reference_file => $self->reference_file,
        );
        unless ($reorder_cmd->execute) {
            die 'failed to execute reorder_cmd';
        }
    }
    unless (-s $self->output_file) {
        die 'output has no size: ' . $self->output_file;
    }

    return 1;
}
