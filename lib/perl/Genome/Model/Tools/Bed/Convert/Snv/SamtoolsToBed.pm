package Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;


class Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
};

sub help_brief {
    "Tools to convert samtools SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv samtools-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {
    return <<EOS
    This is a small tool to take SNV calls in samtools format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self     = shift;
    my $input_fh = $self->_input_fh;

    while (my $line = <$input_fh>) {
        next if $line =~ /^#/;  #mpileup vcf file contains header while pileip output not
        my @tokens = split /\s+/, $line;
        my ($chromosome, $position, $reference, $consensus, $quality, $depth, $gt_info);

        if ($tokens[4] =~ /^[A-Z]+/) {   # 5th column, mpileup is alt variant bases, while pileup is consensus quality 
            ($chromosome, $position, $reference, $consensus, $gt_info) = map{$tokens[$_]}qw(0 1 3 4 9);
            my @vars = split /,/, $consensus; #mpileup use -A option to spit out alternative variant calls
            my $id = 0;
            my %id_vars = (0 => $reference);
            for my $var (@vars) {
                $id++;
                $var = uc $var;
                next if $var eq 'X';  # X should always be the last one in the list
                $id_vars{$id} = $var;
            }
            # 10th column in mpileup vcf output includes allele genotype info 0/1 and 1/1,
            # 0 -> ref, 1 -> 1st in ALT colum, 2 -> 2nd in ALT column
            my ($a1, $a2) = $gt_info =~ /^(\d)\/(\d):/;
            my $str = join '', sort ($id_vars{$a1}, $id_vars{$a2});
            $consensus = Genome::Info::IUB->string_to_iub($str);
            unless ($consensus) {
                $self->warning_message("Failed to get proper variant call from line: $line");
                next;
            }
                
            $quality = sprintf "%2.f", $tokens[5];

            if ($tokens[7] =~ /DP=(\d+)/) {
                $depth = $1;
            } 
            else { 
                $self->warning_message("read depth not found on line $line");
                next;
            }
        }
        else { #pileup format
            ($chromosome, $position, $reference, $consensus, $quality, $depth) = map{$tokens[$_]}qw(0 1 2 3 4 7);
        }
        #position => 1-based position of the SNV
        #BED uses 0-based position of and after the event
        #use consensus quality in bed
        $self->write_bed_line($chromosome, ($position - 1), $position, $reference, $consensus, $quality, $depth);
    }
    return 1;
}


1;
