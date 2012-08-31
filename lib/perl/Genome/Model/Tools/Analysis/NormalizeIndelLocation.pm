package Genome::Model::Tools::Analysis::NormalizeIndelLocation;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Cwd qw( abs_path );
use POSIX;

class Genome::Model::Tools::Analysis::NormalizeIndelLocation {
    is => 'Command',
    has => [
    annotation_input_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => 'List of sites in annotation input format to normalize',
        default => '',
    },
    reference_fasta =>
    {
        type => 'String',
        is_optional => 0,
        doc => 'File from which to grab reference sequence',
    },
    output_file => {
        is => 'Text',
        is_output => 1,
        doc => "Output of left-shifted indels",
    },

    ]
};
#Code should operate as follows
#Look at the indel
#See if we can permute the indel to the left by one base
#Continue until we can no longer do it.
#Report

sub execute {
    my $self=shift;
    my $output_file = $self->output_file;

    # Check that we're on a 64-bit system and can run with the deployed samtools
    unless (POSIX::uname =~ /64/) {
        $self->error_message("Must run on a 64 bit machine");
        die;
    }
    my $indel_file = $self->annotation_input_file;

    my $fh = IO::File->new($indel_file, "r");
    unless($fh) {
        $self->error_message("Couldn't open $indel_file: $!"); 
        return;
    }

    my $chromosome_sequence = undef;   #array of bases in current chromosome
    my $current_chr = undef;
    my $fasta = $self->reference_fasta;

    open(OUTFILE, ">" . $output_file) or die "Can't open outfile: $!\n";
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $ref, $var, @rest) = split /\t/, $line;   #not checking type directly here
        unless($current_chr && $chr eq $current_chr) {
            #load the chromosomal sequence
            $chromosome_sequence = '';
            unless(open(FAIDX,"samtools faidx $fasta $chr |")) {
                $self->error_message("Couldn't pipe samtools faidx");
                return;
            }
            while(my $fasta_line = <FAIDX>) {
                next if $fasta_line =~ /^>/;
                chomp $fasta_line;
                $chromosome_sequence .= $fasta_line;
            }
            unless(close(FAIDX)) {
                $self->error_message("Error reading from samtools faidx pipe");
                return;
            }
            $current_chr = $chr;
        }

        #we should now have a big string containing the reference as well as the position on that string (1-based) from the annotation file
        if($ref =~ /[ACTGN]/ && $var !~ /[ACTGN]/) { 
            #if deletion on reference test with sequence
            #AGTGTGTGTGTA
            #AGTGTGT--GTA want A--GTGTGTGT 
            #original line would be 
            #test	8	9   GT  0

            my @allele = split //,$ref;
            my $preceding_ref_base = substr $chromosome_sequence,$start - 2,1;
            while($allele[-1] eq $preceding_ref_base) {
                #TODO do some error checking here on the coordinates or we will probably screw some poop up
                unshift @allele, pop @allele; #rotate the allele string
                $start--; $stop--;
                $preceding_ref_base = substr $chromosome_sequence,$start - 2,1;
            }
            $ref = join(q{},@allele);
        }
        elsif($var =~ /[ACTGN]/ && $ref !~ /[ACTGN]/) {
            
            #if insertion like say NPM1
            #
            #ATGCAT****GCATG
            #ATGCATGCATGCATG
            #
            #the original line would be
            #test	6	7	0	GCAT
            my @allele = split //,$var;
            my $preceding_ref_base = substr $chromosome_sequence,$start - 1,1;
            while($allele[-1] eq $preceding_ref_base) {
                #TODO do some error checking here on the coordinates or we will probably screw some poop up
                unshift @allele, pop @allele; #rotate the allele string
                $start--; $stop--;
                $preceding_ref_base = substr $chromosome_sequence,$start - 1,1;
            }
            $var = join(q{},@allele);
        }
        else {
            $self->error_message("This line was not an indel: $line");
        }
        print OUTFILE join("\t",$chr,$start,$stop,$ref,$var,@rest),"\n";
    }
    return 1;
}


1;

sub help_brief {
    "Scans a snp file and finds adjacent sites. Then identifies if these are DNPs"
}

sub help_detail {
    <<'HELP';
    Need to fill in help detail
HELP
}

