package Genome::Model::Tools::Nimblegen::CalculateDesignTargetSpace;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use POSIX qw(ceil floor); 

class Genome::Model::Tools::Nimblegen::CalculateDesignTargetSpace {
    is => 'Command',
    has => [
    reference_index => {
        type => 'String',
        is_optional => 0,
        doc => "Samtools index of the reference sequence (To check chromosomal bounds of regions). This is the .fai file associated with a reference sequence fa(sta) file",
    },
    array_design_file => {
        type => 'String',
        is_optional => 1,
        doc => "The array design file in the format we send to Nimblegen. Essentially a 1-based bed-like format. Reads from stdin if unspecified.",
    },
    output_file => {
        type => 'String',
        is_optional => 1,
        doc => "Output file containing results. Output goes to stdout if unspecified.",
    },
    minimum_size => {
        type => 'Integer',
        default => 120,
        is_optional => 1,
        doc => "Minimum size of a target. Less than this will be padded out.",
    },

    ]
};


sub execute {
    my $self=shift;

    my $reference_index = $self->reference_index;

    my $fh = IO::File->new($reference_index,"r");
    unless($fh) {
        $self->error_message("Unable to open the reference sequence index: $reference_index");
        return;
    }

    #read in the index to get the chromosome lengths
    my %chromosome_lengths;
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $length) = split /\t/, $line;
        unless($chr =~ /^chr/) {
            $chr = "chr$chr";
        }
        $chromosome_lengths{$chr} = $length;
    }
    $fh->close;

    my $output_fh;
    if(defined $self->output_file) {
        $output_fh = IO::File->new($self->output_file,"w");
        unless($output_fh) {
            $self->error_message("Unable to open file " . $self->output_file . " for writing.");
            return;
        }
    }
    else {
        $output_fh = IO::File->new_from_fd(fileno(STDOUT),"w");
        unless($output_fh) {
            $self->error_message("Unable to open STDOUT for writing.");
            return;
        }
    }

    my $input_fh;
    if(defined $self->array_design_file) {
        $input_fh = IO::File->new($self->array_design_file,"r");
        unless($input_fh) {
            $self->error_message("Unable to open file ". $self->array_design_file . " for reading.");
            return;
        }
    }
    else {
        $input_fh = IO::File->new_from_fd(fileno(STDIN),"r");
        unless($input_fh) {
            $self->error_message("Unable to open STDIN for reading.");
            return;
        }
    }

    #open a temporary file to sanitize input to bed tools
    my ($temp_file, $temp_file_path) = Genome::Sys->create_temp_file;
    unless($temp_file) {
        $self->error_message("Unable to open a temp file for input file sanitation.");
        return;
    }

    $DB::single = 1;
    my $errors = 0;
    my $non_unique_coverage = 0;
    while(my $line = $input_fh->getline) {
        chomp $line;
        #this should be 1-based.
        my ($chr, $start, $stop) = split /\t/, $line;

        my $length = $stop - $start + 1;
        if( $length < $self->minimum_size) {
            #pad out around feature
            #larger side will be to the left. Arbitratry dlarson decision. Don't know what Nimblegen or any other Vendor actually does in the case of uncentered variants.
            #this may still break if the minimum size is ever smaller than the chromosome.

            my $length_for_pad = $self->minimum_size - $length;
            my $left_pad = ceil($length_for_pad / 2);
            my $right_pad = floor($length_for_pad / 2);
            if($start <= $left_pad) {
                my $hangoff = $left_pad - $start + 1;
                $left_pad -= $hangoff;
                $right_pad += $hangoff;
            }
            if(($chromosome_lengths{$chr} - $stop) < $right_pad) {
                my $hangoff = $right_pad - ($chromosome_lengths{$chr} - $stop);
                $right_pad -= $hangoff;
                $left_pad += $hangoff;
            }
            $start -= $left_pad;
            $stop += $right_pad;
        }

        unless(defined $chr && defined $start && defined $stop) {
            $self->error_message("Missing chromosome, start or start field at line $. from $line\n");
            $errors++;
            next;
        }
        #check and make sure that every line starts with chr
        unless($chr =~ /^chr[0-9XY][0-9]{0,1}/) {
            $self->error_message("Line $. has an invalid chromosome: $line.\nDid you forget to prepend chr?\n");
            $errors++;
            next;
        }
        unless(exists $chromosome_lengths{$chr}) {
            $self->error_message("Unknown chromosome $chr in line $. from $line.\n");
            $errors++;
        }

        #check that every start is less than every stop
        unless($start < $stop) {
            $self->error_message("Line $. contains a start coordinate that is greater than it's stop: $line\n");
            $errors++;
        }

        #check that all coordinates are on chromosomes
        unless(!exists($chromosome_lengths{$chr}) || ($start > 0 && $start <= $chromosome_lengths{$chr})) {
            $self->error_message("$start is not within the limits of $chr in line $.: $line\n");
            $errors++;
        }
        unless(!exists($chromosome_lengths{$chr}) || ($stop > 0 && $stop <= $chromosome_lengths{$chr})) {
            $self->error_message("$stop is not within the limits of $chr in line $.: $line\n");
            $errors++;
        }
        print $temp_file "$chr\t$start\t$stop\n" if $temp_file;
        $non_unique_coverage += ($stop - $start + 1);
        #$non_unique_coverage += ($stop - $start); #for true BED file
    }
    $temp_file->close;
    my $unique_coverage = 0;
    my $temp_output_path = Genome::Sys->create_temp_file_path;
    if(Genome::Model::Tools::BedTools::Merge->execute(input_file => $temp_file_path, output_file => $temp_output_path)) {
        my $merged_output_fh = IO::File->new($temp_output_path);
        unless($merged_output_fh) {
            $self->error_message("Unable to open temporary bedtools mergeBed output file");
            return;
        }
        while(my $line = $merged_output_fh->getline) {
            chomp $line;
            my ($chr, $start, $stop) = split /\t/, $line;
            $unique_coverage += ($stop - $start + 1);
        }
        $merged_output_fh->close;

    }
    else {
        $self->error_message("Error running bedtools mergeBed to determine unique region target space");
        return;
    }

    print $output_fh "File had $errors error(s)\n";
    printf $output_fh "Non-unique region coverage of %0.02f Mbp\n",$non_unique_coverage / 1_000_000;
    printf $output_fh "Unique region coverage of %0.02f Mbp\n",$unique_coverage / 1_000_000;
    $input_fh->close;
    $output_fh->close;

    return 1;

}

sub help_brief {
    "Calculates the Mbp of target space in the Nimblegen Design file.";
}

sub help_detail {
    return <<EOS
This script checks for several common errors in the production of Nimblegen Design files, mostly regarding start and stop coordinates on the chromosomes and the addition of the chr onto the chromosome names.  More usefully, it calculates the target space of the array. The Non-unique region coverage is the total number Mbp without accounting for overlap between regions. The Unique region coverage is the number of Mbp of the genome requested for capture. The current theoretical limit for a Nimblegen custom capture chip is 35 Mbp. 

EOS
}

1;
