package Genome::Model::Tools::Array::CompareTwoGenotypeFiles;

use warnings;
use strict;

use IO::File;
use Genome;

class Genome::Model::Tools::Array::CompareTwoGenotypeFiles {
    is => 'Command',
    has => [
    genotype_file1 => {
        type => 'String',
        is_optional => 0,
        doc => "Genotype file 1, format [chr pos genotypes], no header",
    },
    genotype_file2 => {
        type => 'String',
        is_optional => 0,
        doc => "Genotype file 2, format [chr pos genotypes], no header",
    },
    detailed_report => {
        type => 'String',
        is_optional => 1,
        doc => "A probe-by-probe comparison report between the two files",
    },
    summary_report => {
        type => 'String',
        is_optional => 1,
        doc => "A simple summary of the total number of probes, and how many matched, were LOH, etc.",
    },
    ]
};

sub help_brief {
    "Compare two genotype files of our standard format: [chromosome  position  alleles(usually 2, such as AA,CC,etc.)]. Genotype files typically have no header, so remove headers or adjust the script, or it will probably just call them a 'match'."
}

sub help_detail {
    "Compare two genotype files of our standard format: [chromosome  position  alleles(usually 2, such as AA,CC,etc.)]. Genotype files typically have no header, so remove headers or adjust the script, or it will probably just call them a 'match'."
}

sub execute {
    my $self = shift;

    #parse inputs
    my $geno1 = $self->genotype_file1;
    my $geno2 = $self->genotype_file2;
    my $detailed_report = $self->detailed_report;
    my $output_summary = $self->summary_report;

    #print detailed report?
    my $print_details = int( defined( $detailed_report ));

    #hash for recording results
    my %results;
    my $results = \%results;

    #parse geno1
    my %geno1;
    my $geno1fh = IO::File->new( $geno1 );
    while (my $line = $geno1fh->getline) {
        chomp $line;
        #remove extraneous control characters
        $line =~ s/\s+$//;
        my ($chr, $pos, $alleles) = split /\t/,$line;
        $alleles = sort_alleles($alleles);
        $geno1{$chr}{$pos} = $alleles;
    }
    $geno1fh->close;

    #if printing details, prepare to print
    my $det_rep_fh;
    if ($print_details) {
        $det_rep_fh = IO::File->new( $detailed_report, ">" );
        print $det_rep_fh "Chr\tPos\tGeno_file1\tGeno_file2\tStatus\n";
    }

    #parse geno2 and compare
    my $geno2fh = IO::File->new( $geno2 );
    while (my $line = $geno2fh->getline) {
        chomp $line;

        #remove extraneous control characters
        $line =~ s/\s+$//;

        my ($chr,$pos,$alleles) = split /\t/,$line;

        $alleles = sort_alleles($alleles);
        my $rev_comp_alleles = rev_comp($alleles);
        $rev_comp_alleles = sort_alleles($rev_comp_alleles);

        #proceed with comparison if this chr and pos exists in geno1
        if (exists($geno1{$chr}{$pos})) {

            if ($print_details) {
                my $file_line = join("\t",$chr,$pos,$geno1{$chr}{$pos},$alleles);
                print $det_rep_fh "$file_line\t";
            }

            if ($geno1{$chr}{$pos} eq $alleles) {
                print $det_rep_fh "match\n" if $print_details;
                log_results($results,"matched");
                next;
            }
            elsif ($geno1{$chr}{$pos} eq "--" && $alleles eq "--") {
                print $det_rep_fh "no data file1,file2\n" if $print_details;
                log_results($results,"had no data in either file");
                next;
            }
            elsif ($geno1{$chr}{$pos} eq "--") {
                print $det_rep_fh "no data file1\n" if $print_details;
                log_results($results,"had no data in file1");
                next;
            }
            elsif ($alleles eq "--") {
                print $det_rep_fh "no data file2\n" if $print_details;
                log_results($results,"had no data in file2");
                next;
            }
            elsif ($geno1{$chr}{$pos} eq $rev_comp_alleles) {
                print $det_rep_fh "match rev-complement\n" if $print_details;
                log_results($results,"matched rev-complemented");
                next;
            }
            elsif (is_LOH($geno1{$chr}{$pos},$alleles)) {
                print $det_rep_fh "LOH\n" if $print_details;
                log_results($results,"showed LOH");
                next;
            }
            else {
                print $det_rep_fh "alleles different\n" if $print_details;
                log_results($results,"did not match");
            }
        }
        else {
            if ($print_details) {
                my $file_line = join("\t",$chr,$pos,"Not_Found",$alleles,"site not found in $geno1");
                print $det_rep_fh "$file_line\n";
            }
        }
    }
    $geno2fh->close;
    if ($print_details) {
        $det_rep_fh->close;
    }

    #print results summary
    if (defined($output_summary)) {
        my $summaryfh = IO::File->new( $output_summary, ">" );
        print $summaryfh "There were $results{'total'} probes shared between the 2 files.\n";
        for my $type (sort keys %results) {
            next if $type eq "total";
            my $percent = $results{$type} * 100 / $results{'total'};
            my $rounded = sprintf("%.1f",$percent);
            print $summaryfh "$results{$type} probes ($rounded%) $type.\n";
        }
    }

    return 1;
}

sub log_results {
    my $results = shift;
    my $key = shift;
    $results->{'total'}++;
    $results->{$key}++;
}

sub is_LOH {
    my ($alleles1,$alleles2) = @_;
    my ($a11,$a12) = split //,$alleles1;
    my ($a21,$a22) = split //,$alleles2;
    if ($a11 eq $a12 && ($a11 eq $a21 || $a11 eq $a22)) {
        return 1;
    }
    if ($a21 eq $a22 && ($a21 eq $a11 || $a21 eq $a12)) {
        return 1;
    }
    else {
        return 0;
    }
}

sub rev_comp {
    my $alleles = shift;
    $alleles = reverse($alleles);
    $alleles =~ tr/actgACTG/TGACTGAC/;
    return $alleles;
}

sub sort_alleles {
    my $alleles = shift;
    if($alleles =~ m/NO|NoCall/) {
        return "--";
    }
    my @alleles = split(//,$alleles);
    @alleles = sort @alleles;
    my $sorted_alleles = join("",sort @alleles);
    return $sorted_alleles;
}

1;
