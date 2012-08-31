package Genome::Model::Tools::Bed::ChromSort;

use strict;
use warnings;

use Sort::Naturally;
use Genome;
use Data::Dumper;
use Carp qw/confess/;

use constant BUFSZ => 65536;

class Genome::Model::Tools::Bed::ChromSort {
    is => ['Command'],
    doc => 'Sort a BED file by chromosome write it to a new file',
    has_input => [
        input => {
            is => 'File',
            shell_args_position => 1,
            doc => 'The input BED file to be sorted',
        },
        output => {
            is => 'File',
            shell_args_position => 2,
            doc => 'Where to write the output BED file',
        },
    ],
    has_transient_optional => [
        _input_fh => {
            is => 'IO::File',
            doc => 'Filehandle for the source BED file',
        },
        _output_fh => {
            is => 'IO::File',
            doc => 'Filehandle for the output BED file',
        },
    ]
};

sub help_brief {
    "Sort BED files by chromosome (assumes start and end positions within each are sorted).";
}

sub help_synopsis {
    "gmt bed sort a.bed b.bed";
}

sub help_detail {
    "Sort BED files by chromosome (assumes start and end positions within each are sorted).";
}

sub execute {
    my $self = shift;

    $self->{tmpfiles} = {};
    $self->{tmpfile_count} = 0;
    $self->{open_tmpfile_array} = [];

    eval { $self->initialize_filehandles; };
    if ($@) {
        $self->error_message($@);
        $self->close_filehandles;
        return;
    }

    my $retval = $self->sort;
    $self->close_filehandles;
    return $retval;
}

sub initialize_filehandles {
    my $self = shift;

    if (!defined $self->_input_fh) {
        $self->_input_fh(Genome::Sys->open_file_for_reading($self->input));
    }

    if (!defined $self->_output_fh) {
        $self->_output_fh(Genome::Sys->open_file_for_writing($self->output));
    }

    return 1;
}

sub close_filehandles {
    my $self = shift;

    while (my ($chrom, $tmpfile) = each %{$self->{tmpfiles}}) {
        $tmpfile->{fh}->close;
    }

    close($self->_input_fh) if $self->_input_fh;
    close($self->_output_fh) if $self->_output_fh;

    return 1;
}

sub sort  {
    my $self = shift;
    while (my $line = $self->_input_fh->getline) {
        #svs.hq could contain # lines
        if ($line =~ /^#/) {
            $self->_output_fh->print($line);
            next;
        }
        my ($chrom, $junk) = split("\t", $line, 2);
        my $tmpfile = $self->_get_filehandle_for_chrom($chrom);
        $tmpfile->{fh}->print($line);
    }

    my @chrom_order = nsort keys %{$self->{tmpfiles}};
    for my $chrom (@chrom_order) {
        my $tmpfile = $self->_get_filehandle_for_chrom($chrom);
        my $fh = $tmpfile->{fh};
        $fh->seek(0, 0);

        my $buf;
        while ($fh->read($buf, 65536) > 0) {
            $self->_output_fh->print($buf)
        }

    }
    for my $chrom (keys %{$self->{tmpfiles}}){
        my $fh = $self->{tmpfiles}->{$chrom}->{fh};
        $fh->close if $fh;
        delete $self->{tmpfiles}->{$chrom};
    }

    return 1;
}

sub _get_filehandle_for_chrom {
    my ($self, $chrom) = @_;
    #keep a max of 500 filehandles open here, certain metagenomic references have many thousand contigs and will exceed the max number of available file descriptors
    if (exists $self->{tmpfiles}->{$chrom}) {
        #we have opened a fh for this chrom before
        unless ($self->{tmpfiles}->{$chrom}->{fh}){
            #do nothing unless we have previously closed this filehandle 
            if ($self->{tmpfile_count} > 500){
                #exceeded max open fh count
                #close oldest filehandle
                my $old_chrom = shift @{$self->{open_tmpfile_array}};
                $self->{tmpfiles}->{$old_chrom}->{fh}->close;
                $self->{tmpfile_count}--;
                delete $self->{tmpfiles}->{$old_chrom}->{fh};
            }
            my $path = $self->{tmpfiles}->{$chrom}->{path};
            my $fh = new IO::File("+>>$path") or confess "failed to open temp file $path for chromosome $chrom";
            $self->{tmpfiles}->{$chrom}->{fh} = $fh;
            $self->{tmpfile_count}++;
            push @{$self->{open_tmpfile_array}}, $chrom;
        }
    }else{
        #first time opening an fh for this chrom
        if ($self->{tmpfile_count} > 500){
            #exceeded max open fh count
            #close oldest filehandle
            my $old_chrom = shift @{$self->{open_tmpfile_array}};
            $self->{tmpfiles}->{$old_chrom}->{fh}->close;
            $self->{tmpfile_count}--;
            delete $self->{tmpfiles}->{$old_chrom}->{fh};
        }
        #open new fh
        my $path = Genome::Sys->create_temp_file_path;
        my $fh = new IO::File("+>$path") or confess "failed to open temp file $path for chromosome $chrom";
        $self->{tmpfiles}->{$chrom} = {
            path => $path,
            fh => $fh,
        };
        $self->{tmpfile_count}++;
        push @{$self->{open_tmpfile_array}}, $chrom;
    }
    return $self->{tmpfiles}->{$chrom};
}

1;
