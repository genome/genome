package Genome::Model::Tools::LiftOverMultipleColumns;

use warnings;
use strict;
use IO::File;
use Genome;

class Genome::Model::Tools::LiftOverMultipleColumns {
    is => 'Command',
    has => [
    source_file => {
        is => 'Text',
        doc => 'The file to be translated',
    },
    destination_file => {
        is => 'Text',
        doc => 'Where to output the translated file',
    },
    columns_to_lift => {
        is => 'Text',
        doc => '1-based, underscore-separated, consecutive column numbers to lift, either in the form of 1_2 or 1_2_3 (multiple locations in the file may be submitted in a comma-delimited format, such as "1_2,7_8_9")',
    },
    ],
    has_optional => [
    chain_file => {
        is => 'Text',
        doc => 'The liftOver "chain" file that maps from the source reference to the destination reference. Needed if the lift-direction param is not specified',
    },
    lift_direction => {
        is => 'Text',
        doc => 'shorthand for common lift operations one of: [hg19ToHg18,hg18ToHg19]',
    },
    header_id_string => {
        is => 'Text',
        doc => 'A string to be used in a perl line match which will identify the header lines',
        default => '#'
    },
    ],
};

sub help_detail {
    return <<EOS
    This tool will perform liftOver on multiple sets of consecutive coordinate columns in your input "source" file (ANNOTATION FORMAT IS ASSUMED), and rewrite the input file into the output "destination" file with the lifted coordinates replacing the original coordinates. The coordinates can either be described by colums for chr and then pos (2 columns), or by chr and then two position columns (3 columns total). The locations for the lift to be performed should be indicated by 1-based column numbers separated by underscores, given in the parameter --columns-to-lift, like this: a typical annotation format file would have coordinate columns listed like --columns-to-lift "1_2_3". A varscan file would be lifted with this --columns-to-lift "1_2". A file with multiple lift locations necessary, such as perhaps a BreakDancer output file might look like: --columns-to-lift "3_4_5,7_8_9".
EOS
}

sub execute {

    # parse input arguments
    my $self = shift;
    my $source_file = $self->source_file;
    my $output_file = $self->destination_file;
    my $direction = $self->lift_direction;
    my $chain_file = $self->chain_file;
    my $sets = $self->columns_to_lift;
    my $header_identifier = $self->header_id_string;

    # checks on input params
    unless (defined $direction or defined $chain_file) {
        $self->error_message("Please define either a direction for lifting or a chain file to be used by liftOver.");
        return;
    }

    # declarations
    my @sets = split ",",$sets;
    my @zb_chr_cols;
    my %files;
    # $files{$set}{'old_wfh'} = file to be written to and then lifted
    # $files{$set}{'lifted_rfh'} = file to be read after liftover for printing
    my $temp_dir = Genome::Sys->create_temp_directory();

    # create filehandles to write old coordinates into before lifting over
    for my $set (@sets) {
        $files{$set}{'old_wfh'} = new IO::File $temp_dir."/".$set,"w";
        my ($chr_col) = split "_",$set;
        push @zb_chr_cols,$chr_col-1;
    }

    # read through input file first time
    my $source_fh = new IO::File $source_file,"r";
    while (my $line = $source_fh->getline) {

        # if header, next
        if ($line =~ /$header_identifier/i) { next; }

        # otherwise, write necessary columns to files for lifting over
        chomp $line;
        my @fields = split /\t/,$line;
        for my $set (@sets) {
            my ($chr_col) = split "_",$set;
            $chr_col--;
            my $count = ($set =~ s/_/_/g);
            if ($count == 1) {
                $files{$set}{'old_wfh'}->print(join("\t",$fields[$chr_col],$fields[$chr_col+1],$fields[$chr_col+1]+1,"C","T") . "\n");
            }
            elsif ($count == 2) {
                $files{$set}{'old_wfh'}->print(join("\t",$fields[$chr_col],$fields[$chr_col+1],$fields[$chr_col+2],"C","T") . "\n");
            }
        }
    }

    # close open filehandles
    $source_fh->close;
    for my $set (@sets) { $files{$set}{'old_wfh'}->close; }

    # run liftover on all temp files
    for my $set (@sets) {
        my $infile = $temp_dir."/".$set;
        my $outfile = $infile.".lifted";
        my $unmapped = $infile.".unmapped";
        my $lift_cmd = Genome::Model::Tools::LiftOver->create(
            source_file => $infile,
            destination_file => $outfile,
            lift_direction => $direction,
            input_is_annoformat => 1,
            unmapped_file => $unmapped,
            chain_file => $chain_file,
        );
        unless ($lift_cmd->execute) {
            $self->error_message("Lifted file compromised. Aborting...");
            return;
        }
        $lift_cmd->delete;

        # check to make sure that all sites were mapped. it would require a lot more code to make sure that upon finding a failure,
        # if there are multiple "sets" of columns in particular, the entire row will correctly be removed... an option might be to just
        # not include this row by putting the line in a hash and when seeing it the second time, doing "next;". This still would require more code...
        # the current, lazy check:
        my $unmapped_wc = `wc -l $unmapped`;
        ($unmapped_wc) = split /\s/,$unmapped_wc;
        if ($unmapped_wc > 0) {
            $self->error_message("Some sites unmapped from set $set. Aborting operation.");
            next;
        }
    }

    # open all lifted files
    for my $set (@sets) {
        $files{$set}{'lifted_rfh'} = new IO::File $temp_dir."/".$set.".lifted","r" or die "Couldn't read lifted file for set $set\n";
    }

    # open output filehandle
    my $output_fh = new IO::File $output_file,"w";
    unless ($output_fh) {
        $self->error_message("Could not open output file $output_file.");
        return;
    }

    # read through input file a second time, printing output to output file
    $source_fh = new IO::File $source_file,"r";
    while (my $line = $source_fh->getline) {

        # handle header
        if ($line =~ /$header_identifier/) { print $output_fh $line; next; }

        for my $set (@sets) {
            my $line = $files{$set}{'lifted_rfh'}->getline;
            my $count = ($set =~ s/_/_/g);
            if ($count == 1) { my ($chr,$start) = split /\t/,$line; $files{$set}{'current_line'} = join("\t",$chr,$start); }
            if ($count == 2) { my ($chr,$start,$stop) = split /\t/,$line; $files{$set}{'current_line'} = join("\t",$chr,$start,$stop); }
        }

        # write necessary columns to files for lifting over
        chomp $line;
        my @fields = split /\t/,$line;
        for (my $i=0; $i <= $#fields; $i++) {
            unless (grep { /^$i$/ } @zb_chr_cols) {
                print $output_fh $fields[$i];
            }
            else {
                my $possible_set_1 = ($i+1) . "_" . ($i+2);
                my $possible_set_2 = $possible_set_1 . "_" . ($i+3);
                if (exists $files{$possible_set_1}) { print $output_fh $files{$possible_set_1}{'current_line'}; $i+=1; }
                if (exists $files{$possible_set_2}) { print $output_fh $files{$possible_set_2}{'current_line'}; $i+=2; }
            }
            print $output_fh "\t" unless ($i >= $#fields);
        }
        print $output_fh "\n";
    }

    return 1;
}

1;
