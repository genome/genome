package Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
};

sub help_brief {
    "Tools to convert samtools indel format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert indel samtools-to-bed --source indels_all_sequences --output indels_all_sequences.bed
EOS
}

sub help_detail {
    return <<EOS
    This is a small tool to take indel calls in samtools format and convert them to a common BED format (using the first four columns).
EOS
}

sub process_source {
    my $self     = shift;
    my $input_fh = $self->_input_fh;

    while (my $line = <$input_fh>) {
        next if $line =~ /^#/;  #mpileup vcf gets header
        my ($chromosome, $position, $consensus_quality, $read_depth, $reference_bases, $variant_bases);
        my %indel_calls;

        my @tokens = split /\s+/, $line;

        if ($tokens[4] =~ /^[A-Z]+/) {  #mpileup
            unless ($tokens[7] =~ /INDEL/) {
                $self->error_message("Indel line:\n$line does not have INDEL as token in 8th column");
                die;
            }

            ($chromosome, $position, $reference_bases, $variant_bases) = map{$tokens[$_]}qw(0 1 3 4);
            $consensus_quality = sprintf "%2.f", $tokens[5]; 

            my @indel_list = split /,/, $variant_bases;
            for my $indel_str (@indel_list) {
                my ($conv_indel_str, $pos) = $self->convert_indel_string($position, $reference_bases, $indel_str);
                unless ($conv_indel_str and $pos) {
                    $self->warning_message("Failed to convert indel string: $indel_str from $line");
                    next;
                }
                $indel_calls{$conv_indel_str} = $pos;
            }

            if ($tokens[7] =~ /DP=(\d+)/) {
                $read_depth = $1;
            } 
            else { 
                $self->warning_message("read depth not found on line $line");
                next;
            }
        }
        else {  #pileup
            ($chromosome, $position, $consensus_quality, $read_depth) = map{$tokens[$_]}qw(0 1 4 7);
            next unless $tokens[2] eq '*'; #pileup indel format includes reference lines as well
            %indel_calls = (
                $tokens[8] => $position, 
                $tokens[9] => $position,
            );
        }

        for my $indel (sort keys %indel_calls) {
            next if $indel eq '*'; # For pileup indel: Indicates only one indel call...and this isn't it!
            my ($reference, $variant, $start, $stop);

            #samtools pileup reports the position before the first deleted base or the inserted base ... so the start position is already correct for bed format
            $start = $indel_calls{$indel};
            if (substr($indel, 0, 1) eq '+') {
                $reference = '*';
                $variant   = substr($indel, 1);
                $stop      = $start; #Two positions are included-- but an insertion has no "length" so stop and start are the same
            } 
            elsif (substr($indel, 0, 1) eq '-') {
                $reference = substr($indel, 1);
                $variant   = '*';
                $stop      = $start + length($reference);
            } 
            else {
                $self->warning_message("Unexpected indel format encountered ($indel) on line:\n$line");
                #return; skip wrong indel format line instead of failing for now
                next;
            }
            $self->write_bed_line($chromosome, $start, $stop, $reference, $variant, $consensus_quality, $read_depth);
        }
    }
    return 1;
}

#vcf indel string need to be converted to bed friendly and handle rare cases like below.
#The rule to get the correct indels: chop the most consecutive matching bases against 
#ref bases from left side and chop the most consecutive matching bases from right sides,
#the leftover is the indels.

#17    71072555    .    AGGG    AAGGGG,AGGGG,AAAGGGG

#17    71072555    71072555    */AG
#17    71072558    71072558    */G
#17    71072555    71072555    */AAG

sub convert_indel_string {
    my ($self, $pos, $ref_str, $indel_str) = @_;
    my ($ref_length, $indel_length) = (length($ref_str), length($indel_str));

    my $left_match_count  = 0; 
    my $right_match_count = 0;

    if ($ref_length == $indel_length) {
        $self->warning_message("ref and indel str have the same length: $ref_str, $indel_str");
        return; #skip this for now
    }

    my ($min_length, $max_str, $type);

    if ($ref_length > $indel_length) {
        ($min_length, $max_str, $type) = ($indel_length, $ref_str, '-');
    }
    else {
        ($min_length, $max_str, $type) = ($ref_length, $indel_str, '+');
    }

    my @ref_bases   = split //, $ref_str;
    my @indel_bases = split //, $indel_str;

    for my $i (0..$min_length-1) {
        last unless $ref_bases[$i] eq $indel_bases[$i];
        $left_match_count++;
    }

    my $conv_pos = $pos + $left_match_count - 1; #the pos right before the indel

    if ($left_match_count == $min_length) {  #match all the bases of shorter str
        return ($type . substr($max_str, $left_match_count), $conv_pos);
    }
    else {  
        my @rev_ref_bases   = reverse @ref_bases;
        my @rev_indel_bases = reverse @indel_bases;

        #right match count <= str leftover after chop left match count
        my $right_match_index_limit = $min_length - $left_match_count - 1;

        for my $j (0..$right_match_index_limit) {
            last unless $rev_ref_bases[$j] eq $rev_indel_bases[$j];
            $right_match_count++;
        }

        my $max_length   = length $max_str;
        my $indel_length = $max_length - $left_match_count - $right_match_count;
        return ($type . substr($max_str, $left_match_count, $indel_length), $conv_pos);
    }
}


1;
