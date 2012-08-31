package Genome::Model::Tools::Sam::FilterSoftclipMismatch;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Sam::FilterSoftclipMismatch{
    is => 'Genome::Model::Tools::Sam',
    has => [
        input_file => {
            is => 'Path',
            doc => 'SAM or BAM input file to filter.  File extension controls format - a .bam is expected to be a BAM file.  Files with no extension are assumed to be SAM.'
        },
        output_file => {
            is => 'Text',
            doc => 'Merged SAM or BAM file to write.  If the fjjjjjjjjilename ends with .bam, it is written in BAM format.  SAM format is used otherwise.  Use - to write to STDOUT.'
        },
        threshold => {
            is => 'Int',
            doc => 'percent threshold to filter out reads with mismatches and softclipping.  If the mismatch + softclipped bases is greater than this percent of read length, it will not appear in the output file',
        },
    ],
};

sub execute{
    my $self = shift;

    unless(-e $self->input_file){
        $self->error_message("input file ".$self->input_file." doesn't exist!");
        die $self->error_message();
    }

    my $in_fh = $self->open_bamsam_in($self->input_file);
    unless ($in_fh){
        $self->error_message("couldn't create input filehandle for ".$self->input_file);
        die $self->error_message();
    }

    my $out_fh = $self->open_bamsam_in($self->output_file);
    unless ($out_fh){
        $self->error_message("couldn't create output filehandle for ".$self->output_file);
        die $self->error_message();
    }

    my $name_col = 0;
    my $flag_col = 1;
    my $cigar_col = 5;
    my $seq_col = 9;

    while (my $line = $in_fh->getline){
        if(substr($line, 0, 1) eq '@'){
            $out_fh->print($line);
            next;
        }

        my @s = split(/\t/, $line);

        my $nm;
        unless($s[$flag_col] & 4) {
            # mapped
            if ($self->softclip_mismatch){
                #if the softclip_mismatch arg is present, add the softclip(and hardclip) bases from the cigar string to the mismatch number
                my $cigar = $s[$cigar_col];
                my ($start_clip, $stop_clip) = $cigar =~ /^(\d+)[hs].*(\d+)[hs]$/;
                $start_clip ||=0;
                $stop_clip ||=0;
                $nm += $start_clip+$stop_clip;
            }
            my $nm_pos = -1;
            foreach my $tag (@s[11..$#s]) {
                $nm_pos = index($tag,'nm:i:');
                if ($nm_pos != -1) {
                    $nm += substr($tag,$nm_pos+5);
                    last;
                }
            }
            if ($nm_pos == -1) {
                $self->error_message('alignment index ' . $s[$name_col] . ' in "' . $self->input_file . "\" is mapped but lacks an nm tag\n$line");
                die $self->error_message();
            }
            my $len = length($s[$seq_col]);
            if ($nm/$len < $self->threshold){
                $out_fh->print 
            }
        }
    }
}
    
1;

