package Genome::Model::PhenotypeCorrelation::Command::FilterCorrelationResults;

use Digest::MD5;
use Genome;
use POSIX;
use Sort::Naturally qw/nsort/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::FilterCorrelationResults {
    is => "Genome::Command::Base",
    doc => "Remove lines from a variant matrix based on MAF",
    has => [
        input_file => {
            is => "File",
            doc => "The clinical correlation output to filter.",
        },
        output_file => {
            is => "File",
            doc => "Path for the output file.",
        },
        per_site_report_file => {
            is => "File",
            doc => "Path for the per site report from a run of joinx vcf-report. This will be used to get the MAF.",
        },
        delimiter => {
            is => 'String',
            doc => 'file delimiter to parse out fields since they are not consistent.',
            default => "\t",
        },
        minimum_maf => {
            is => "Number",
            doc => "The minimum MAF for a site to be retained. If there are more than two alleles the higher of the two minor allele frequencies will be tested."
        },
    ],
};

sub execute {
    my ($self) = @_;

    $DB::single=1;
    my $delim = $self->delimiter;
    my $input_file = $self->input_file;
    my $output_file = $self->output_file;

    die "Can't filter an empty correlation result!" unless -s $input_file;
    my %markers_to_retain = map { $_ => 1 } $self->markers_to_retain;

    my $ifh = Genome::Sys->open_file_for_reading($input_file);
    my $header_line = <$ifh>;
    chomp $header_line;
    my @header_fields = split($delim, $header_line);
    #open output file
    my $ofh = Genome::Sys->open_file_for_writing($output_file);
    print $ofh $header_line, "\n";

    while (<$ifh>) {
        chomp;
        my %fields;
        @fields{@header_fields} = split($delim);
        if( exists($markers_to_retain{$fields{'x'}}) ) {
            print $ofh $_,"\n";
        }
    }

    $ifh->close();
    return 1;
}


sub markers_to_retain {
    my $self = shift;
    my $per_site_report_file = $self->per_site_report_file;
    my $minimum_maf = $self->minimum_maf;
    
    my @sites_to_retain;

    my $fh = Genome::Sys->open_file_for_reading($per_site_report_file);
    my $header = $fh->getline;
    chomp $header;
    my @header_fields = split /\t/, $header;

    while(my $l = $fh->getline) {
        chomp $l;
        my %line;
        @line{@header_fields} = split "\t", $l;

        my ($highest, $second_highest) = (0,0);
        for my $af (split ",", $line{'ByAltAlleleFreq'}) {
            if($af > $highest) {
                ($highest, $second_highest) = ($af, $highest);
            }
            elsif($af > $second_highest) {
                $second_highest = $af;
            }
        }
        if($second_highest >= $minimum_maf) {
            #going to retain this marker
            my $id = $self->_generate_variant_id(@line{qw( Chrom Pos Ref Alt )});
            push @sites_to_retain, $id;
        }
    }
    return @sites_to_retain;
}


sub _generate_variant_id {
    my ($self, $chrom, $pos, $ref, $alts) = @_;
    if($chrom =~ /^(\d+)$/) {
        $chrom = "X$chrom";
    }
    my $id = join("_",$chrom, $pos,$ref,$alts);
    #remove out any commas
    $id =~ s/,/./g;
    return $id;
}

1;
