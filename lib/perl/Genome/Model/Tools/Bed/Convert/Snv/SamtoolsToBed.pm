package Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use Genome::Utility::Vcf qw(parse_vcf_line get_vcf_header get_samples_from_header);


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
    my ($header, @header_array, @sample_names);

    while (my $line = <$input_fh>) {
        # if this is a vcf file, grab the header for later
        if ($line =~ /^#/) {
            unless ($header) {
                $header = get_vcf_header($self->source);
                @header_array = split "\n", $header;
                @sample_names = get_samples_from_header(\@header_array);
            }
            next;
        }
        my @tokens = split /\s+/, $line;
        my ($chromosome, $position, $reference, $consensus, $quality, $depth, $gt_info);

        # If this is a VCF file...
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
                die $self->error_message("Failed to get proper variant call from line: $line");
            }

            $quality = sprintf "%2.f", $tokens[5];

            my $parsed_line = parse_vcf_line($line, \@sample_names);

            # Check both info and sample fields for depth
            if ($parsed_line->{info} ne "." and $parsed_line->{info}{"DP"}) {
                $depth = $parsed_line->{info}{"DP"};
            } elsif ($parsed_line->{sample}{$sample_names[0]}{"DP"}) {
                $depth = $parsed_line->{sample}{$sample_names[0]}{"DP"};
            }
            else { 
                die $self->error_message("Read depth not found on line $line");
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
