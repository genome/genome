package Genome::Model::PhenotypeCorrelation::Command::FilterCorrelationResults;

use Digest::MD5;
use Genome;
use POSIX;
use Sort::Naturally qw/nsort/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::FilterCorrelationResults {
    is => "Command::V2",
    doc => "Remove lines from a variant matrix based on MAF",
    has_input => [
        input_file => {
            is => "File",
            doc => "The clinical correlation output to filter.",
        },
        output_file => {
            is => "File",
            is_output => 1,
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
            doc => "The minimum MAF for a site to be retained. If there are more than two alleles the higher of the two minor allele frequencies will be tested.",
            is_optional => 1,
        },
        maximum_missing_rate => {
            is => "Number",
            doc => "The maximum rate of missing calls for a site to be retained. This does not differentiate between alleles.",
            is_optional => 1,
        },
        maximum_filtered_rate => {
            is => "Number",
            doc => "The maximum rate of filtered calls for a site to be retained. This does not differentiate between alleles.",
            is_optional => 1,
        },
        maximum_excluded_rate => {
            is => "Number",
            doc => "The maximum rate of filtered + missing calls for a site to be retained. This does not differentiate between alleles.",
            is_optional => 1,
        },
        variant_id_column => {
            is => "Text",
            doc => "The name of the variant id column in the input file",
            default_value => "x",
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
    my $var_id_col = $self->variant_id_column;
    unless (grep {$_ eq $var_id_col} @header_fields) {
        die "Variant id column '$var_id_col' not found in input file header ($input_file): $header_line.\n";
    }

    #open output file
    my $ofh = Genome::Sys->open_file_for_writing($output_file);
    print $ofh $header_line, "\n";


    while (<$ifh>) {
        chomp;
        my %fields;
        @fields{@header_fields} = split($delim);

        #restrict the variant id to just the chromosome and position
        ($fields{$var_id_col}) = $fields{$var_id_col} =~ m/(^[^_]+_\d+)_/;
        
        if( exists($markers_to_retain{$fields{$var_id_col}}) ) {
            print $ofh $_,"\n";
        }
    }

    $ifh->close();
    return 1;
}


sub markers_to_retain {
    my $self = shift;
    my $per_site_report_file = $self->per_site_report_file;
    
    my @sites_to_retain;

    my $fh = Genome::Sys->open_file_for_reading($per_site_report_file);
    my $header = $fh->getline;
    chomp $header;
    my @header_fields = split /\t/, $header;

    while(my $l = $fh->getline) {
        chomp $l;
        my %line;
        @line{@header_fields} = split "\t", $l;

        #apply filters
        next if defined $self->minimum_maf and $self->_fails_maf($self->minimum_maf, %line);
        next if defined $self->maximum_missing_rate and $self->_fails_missingness($self->maximum_missing_rate, %line);
        next if defined $self->maximum_filtered_rate and $self->_fails_filtered_rate($self->maximum_filtered_rate, %line);
        next if defined $self->maximum_excluded_rate and $self->_fails_excluded_rate($self->maximum_excluded_rate, %line);

        #going to retain this marker
        my $id = $self->_generate_variant_id(@line{qw( Chrom Pos Ref Alt )});
        push @sites_to_retain, $id;
    }
    return @sites_to_retain;
}


sub _generate_variant_id {
    my ($self, $chrom, $pos, $ref, $alts) = @_;
    if($chrom =~ /^(\d+)$/) {
        $chrom = "X$chrom";
    }
    my $id = join("_",$chrom, $pos);
    return $id;
}

sub _fails_maf {
    my ($self, $minimum_maf, %line) = @_;

    #verify we can do what we want!
    die "No ByAltAlleleFreq column in header\n" unless(exists($line{'ByAltAlleleFreq'}));
    
    my ($highest, $second_highest) = (0,0);
    for my $af (split ",", $line{'ByAltAlleleFreq'}) {
        if($af > $highest) {
            ($highest, $second_highest) = ($af, $highest);
        }
        elsif($af > $second_highest) {
            $second_highest = $af;
        }
    }
    return ($second_highest < $minimum_maf);
}

sub _fails_missingness {
    my ($self, $max_missingness, %line) = @_;
    #verify we can do what we want!
    die "No TotalSamples and/or NumberMissing column in header\n" unless(exists($line{'TotalSamples'}) and exists($line{'NumberMissing'}));

    my $missingness = $line{'NumberMissing'} / $line{'TotalSamples'};

    return ($missingness > $max_missingness);
}

sub _fails_filtered_rate {
    my ($self, $max_filtered_rate, %line) = @_;
    #verify we can do what we want!
    die "No TotalSamples and/or NumberFiltered column in header\n" unless(exists($line{'TotalSamples'}) and exists($line{'NumberFiltered'}));

    my $filtered_rate = $line{'NumberFiltered'} / $line{'TotalSamples'};

    return ($filtered_rate > $max_filtered_rate);
}


sub _fails_excluded_rate {
    my ($self, $max_excluded_rate, %line) = @_;
    #verify we can do what we want!
    die "No TotalSamples and/or NumberMissing NumberFiltered column in header\n" unless(exists($line{'TotalSamples'}) and exists($line{'NumberFiltered'}) and exists($line{'NumberMissing'}));

    my $excluded_rate = ($line{'NumberMissing'} + $line{'NumberFiltered'}) / $line{'TotalSamples'};

    return ($excluded_rate > $max_excluded_rate);
}

1;
