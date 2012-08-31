package Genome::Model::Tools::Snp::Match;

use strict;
use warnings;

use Genome;
use Command;
use FileHandle;

class Genome::Model::Tools::Snp::Match {
    is => 'Command',
    has => [
    f => {
        type => 'String',
        is_optional => 0,
        doc => "File of (minimally) chromosome and position",
    },
    input_files => {
        shell_args_position => 1,
        is_many => 1,
        doc => 'The files to query',
    },

    ],
    has_optional => [
    v => {
        type => 'Boolean',
        doc => 'when used, report only those positions not in the file-of-queries',
        default => 0,
    },
    fdelimiter => {
        type => 'Character',
        doc => 'The field delimiter separating fields in the file-of-queries',
        default => "\t",
    },
    idelimiter => {
        type => 'Character',
        doc => 'The field delimiter separating files in the input files to be searched',
        default => "\t",
    },
    list_file_of_match => {
        type => 'Boolean',
        doc => 'Whether or not to report the name of the file a match is found in',
        default => '0',
    },
    start_and_stop => {
        type => 'Boolean',
        doc => 'assume that there is a third column which is the stop position. Match on both start and stop',
        default => 0,
    },
    clean_whitespace => {
        type => 'Boolean',
        doc => 'remove any whitespace in the lines',
        default => 0,
    }

    ],
};

sub help_brief {
    "Searches files by chromosome and position assuming that all are of a Samtools pileup or annotation input like format. Uses a hash. Be aware of memory requirements.";
}

sub help_detail {
}

sub execute {
    my $self=shift;
    my $fseparator = $self->fdelimiter;
    my $iseparator = $self->idelimiter;
    my $list_file_of_match = $self->list_file_of_match;
    my $clean_whitespace = $self->clean_whitespace;
    my $not_in_f = $self->v;
    my $has_start_and_stop = $self->start_and_stop;

    my %query_positions;

    #read in the query positions
    my $handle = new FileHandle;
    $handle->open($self->f,"r");
    unless($handle) {
        $self->error_message("Unable to open query file: " . $self->f);
        return;
    }

    while(my $line = $handle->getline) {
        chomp $line;
        my ($chromosome, $start, $stop,) = split $fseparator, $line;
        $stop = "" unless defined $stop;
        map {$_ =~ s/\s+//g } ($chromosome, $start, $stop,) if $clean_whitespace;
        $query_positions{$chromosome}{$start}{$stop} = 1;
    }
    $handle->close;

    foreach my $input ($self->input_files) {
        $handle->open($input,"r");
        unless($handle) {
            $self->error_message("Unable to open input file: $input");
        }

        if($list_file_of_match) {
            print $input, ":\n";
        }

        while(my $line = $handle->getline) {
            chomp $line;
            #assuming a snp file
            my ($chromosome, $start, $stop,) = split $iseparator, $line;
            map {$_ =~ s/\s+//g } ($chromosome, $start, $stop) if $clean_whitespace;

            if($not_in_f) {
                if($has_start_and_stop) {
                    print $line, "\n" if(!exists($query_positions{$chromosome}{$start}{$stop}));
                }
                else {
                    print $line, "\n" if(!exists($query_positions{$chromosome}{$start}));
                }
            }
            else {
                if($has_start_and_stop) {
                    print $line, "\n" if(exists($query_positions{$chromosome}{$start}{$stop}));
                }
                else {
                    print $line, "\n" if(exists($query_positions{$chromosome}{$start}));
                }
            }
        }
        $handle->close;
    }

    return 1;
}


1;




