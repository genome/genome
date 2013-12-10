package Genome::Model::Tools::Annotate::AddEvsMaf;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::AddEvsMaf {
    is => 'Command',
    has_input => [
        anno_file => {
            doc => 'The path to the annotation file.',
        },
        vcf_file => {
            doc => 'The path to the VCF file downloaded from the exome variant server',
            example_values => ["/gscmnt/gc13001/info/build_merged_alignments/detect-variants--linus226.gsc.wustl.edu-cmiller-9722-0a144ca3948f434f8f991e524da9ec89/snvs.hq"],
        },
        output_file => {
            doc => 'The path to use for the output file.'
        },
    ],
};

sub help_synopsis {
    q(gmt annotate add-evs-maf --anno-file test.annotate.top --vcf-file=snvs.annotated.vcf.gz --output-file=test.annotate.evs.top);
}

sub help_detail {
    q(This tool appends allele frequency information from the Exome Variant Server to the end of an annotation file as three extra columns. Missing or undefined values are denoted by '-'.);
}

sub execute {
    my $self = shift;

    my $anno_file   = $self->anno_file;
    my $vcf_file    = $self->vcf_file;
    my $output_file = $self->output_file;
    if (-e $output_file) {
        die "output_file already exists: $output_file\n";
    }

    # we have to handle the case where there are multiple anno lines at the same position, so
    # we store every annotation position (chr/pos), store the dbsnp info associated with it,
    # then re-read the annotation file to append the dbsnp info.  Also has the benefit of
    # ensuring that the annotation file output is in the same order as the info

    my $annotation = store_annotation($anno_file);

    #read through the VCF, storing values for any annotation lines that match
    my %vcf_vals;
    my $vcf_fh;
    #handle gzipped or non-gzipped VCF
    if($vcf_file =~ /.gz$/){
        $vcf_fh = Genome::Sys->open_gzip_file_for_reading($vcf_file);
    } else {
        $vcf_fh = Genome::Sys->open_file_for_reading($vcf_file);
    }

    while (my $line = $vcf_fh->getline) {
        chomp $line;
        next if $line =~ /^\#/;

        my ($chr, $pos, $rsID, $ref, $var, $INFO) = (split(/\t/, $line))[0..4, 7];

        my $MAF = "-";
        my $AA_MAF = "-";
        my $EA_MAF = "-";
        if($INFO =~ /[^_]MAF=([0-9.]+)/){
            $MAF = "$1"
        }
        if($INFO =~ /EA_MAF=([0-9.]+)/){
            $EA_MAF = "$1"
        }
        if($INFO =~ /AA_MAF=([0-9.]+)/){
            $AA_MAF = "$1"
        }

        my $key = site_key($chr, $pos);

        if(defined( $annotation->{$key}) ){
            # there is some assumed symmetry between $var and $dbSNPinfo
            my @var_alleles = split(/, /, $var);

            foreach my $var_allele (@var_alleles){
                if(defined($annotation->{$key}->{$var_allele})){
                    if($annotation->{$key}->{$var_allele} ne "0"){
                        $vcf_vals{$key}{$var_allele}{"MAF"} = $MAF;
                        $vcf_vals{$key}{$var_allele}{"EA_MAF"} = $EA_MAF;
                        $vcf_vals{$key}{$var_allele}{"AA_MAF"} = $AA_MAF;
                    }
                }
            }
        }
    }
    $vcf_fh->close();

    print_annotation($anno_file, $output_file, \%vcf_vals);

    return 1;
}

sub store_annotation{
    my $anno_file = shift;

    my $annotation = {};
    $annotation->{"header"}->{"header"} = 0;

    my $anno_fh = Genome::Sys->open_file_for_reading($anno_file);
    while (my $line = $anno_fh->getline) {
        chomp $line;

        if($line =~ /^chromosome/){
            $annotation->{"header"}->{"header"} = $line;
            next;
        }

        my @list = split(/\t/, $line);
        my ($chr, $pos, $var_allele) = (split(/\t/, $line))[0, 1, 4];
        my $key = site_key($chr, $pos);

        $annotation->{$key}->{$var_allele} = $line;
    }
    $anno_fh->close;
    return $annotation;
}

sub print_annotation{
    my $anno_file = shift;
    my $output_file = shift;
    my $vcf_vals = shift;

    my $output_fh = Genome::Sys->open_file_for_writing($output_file);

    my $anno_fh = Genome::Sys->open_file_for_reading($anno_file);
    while (my $line = $anno_fh->getline) {
        chomp $line;

        if($line =~ /^chromosome/){
            print $output_fh $line . "\tEA_MAF\tAA_MAF\tAll_MAF\n";
            next;
        }

        my @list = split(/\t/, $line);
        my ($chr, $pos, $var) = (split(/\t/, $line))[0, 1, 4];
        my $key = site_key($chr, $pos);

        my $suffix = "-\t-\t-";
        if(defined($vcf_vals->{$key}->{$var})){
            $suffix = $vcf_vals->{$key}->{$var}->{"EA_MAF"};
            $suffix .= "\t" . $vcf_vals->{$key}->{$var}->{"AA_MAF"};
            $suffix .= "\t" . $vcf_vals->{$key}->{$var}->{"MAF"};
        }
        print $output_fh $line . "\t" . $suffix . "\n";

    }
    $anno_fh->close;
    $output_fh->close;
}


sub site_key {
    my ($chr, $pos) = @_;
    return join('_', ($chr, $pos));
}

1;
