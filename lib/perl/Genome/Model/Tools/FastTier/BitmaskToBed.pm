package Genome::Model::Tools::FastTier::BitmaskToBed;

use strict;
use warnings;

use Genome;
use IO::File;
use Bit::Vector;

class Genome::Model::Tools::FastTier::BitmaskToBed {
    is => 'Command',
    has => [
        output_file => {
            type => 'Text',
            is_output => 1,
            doc => 'Output file for the bed tier',
        },
        bitmask => {
            type => 'Text',
            is_input => 1,
            doc => 'The input bitmask',
        },
    ],
};


sub help_brief {
    "Used to turn bitmask files into bed files"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools fast-tier bitmask-to-bed ...    
EOS
}

sub execute {
    my $self = shift;

    my $output = $self->output_file;
    my $temp_output = Genome::Sys->create_temp_file_path;
    my $bitmask = $self->bitmask;
    my $out_fh = Genome::Sys->open_file_for_writing($temp_output);




    unless($bitmask) {
        die("nobitmask");
    }
    #do some stuff to read this from a file without making it suck
    my $in_fh = IO::File->new($bitmask,"<:raw");
    unless($in_fh) {
        die("cant.read.from.".$self->bitmask);
    }
    my $read_string;
    sysread $in_fh, $read_string, 4;
    my $header_length = unpack "N", $read_string;
    sysread $in_fh, $read_string, $header_length;
    my $header_string = unpack "a*",$read_string;
    my %genome = split /\t/, $header_string; #each key is the name, each value is the size in bits

    #now read in each one

    for my $chr (sort keys %genome) {
        my $index=0;
        $genome{$chr} = Bit::Vector->new($genome{$chr}); #this throws an exception if it fails. Probably should be trapped at some point in the future
        sysread $in_fh, $read_string, 4;
        my $chr_byte_length = unpack "N", $read_string;
        my $chr_read_string;
        sysread $in_fh, $chr_read_string, $chr_byte_length;
        my $chr_length = $genome{$chr}->Size();
        my $start=0;
        my $chr_length2=  $chr_byte_length*8;
        $genome{$chr}->Block_Store($chr_read_string);
        while(1) {
            my $bit = $genome{$chr}->bit_test($index);
            if($bit && $start) {
                $DB::single=1 if ($start==67167422); 
            }
            elsif($bit) {
                $start=$index;
            }
            elsif(!$bit && $start) {
                my $end=$index-1;
                $start = $start - 1;
                if($start<0){
                    $start = 0;
                }
                $out_fh->print("$chr\t$start\t$end\n");
                $start = 0;
            }
            if(($index+1 >= $chr_length) && $start) {
                $out_fh->print("$chr\t$start\t$chr_length\n");
                $start = 0;
                $index=0;
                last;
            }
            elsif($index+1 >= $chr_length) {
                $index=0;
                last;
            }
            $index++;
        }
    }
    $in_fh->close;
    $out_fh->close;

    if(-s $temp_output){
        my @input = ($temp_output);
        my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create(input_files => \@input, output_file => $output);
        unless($sort_cmd->execute){
            die $self->error_message("Could not complete sorting operation on bed output!");
        }
    }
    else {
        my $cmd = "touch $output";
        Genome::Sys->shellcmd(cmd => $cmd);
    }

    return 1;
}
