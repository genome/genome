package Genome::Model::Tools::Annotate::AddRsid;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::AddRsid {
    is => 'Command',
    has_input => [
        anno_file => {
            doc => 'The path to the annotation file.',
        },
        vcf_file => {
            doc => 'The path to the VCF file containing the RS ID information.',
        },
        output_file => {
            doc => 'The path to use for the output file.'
        },
    ],
};

sub help_synopsis {
    q(gmt annotate add-rsid --anno-file test.annotate.top --vcf-file=snvs.annotated.vcf.gz --output-file=test.annotate.rsid.top);
}

sub help_detail {
    q(This tool appends RS ID and GMAF information to the end of the annotation file as 2 extra columns.  It grabs the dbSNP id and GMAF information from the annotated snv vcf file in the somatic varation build variants subdirectory.  Missing or undefined values are denoted by '-'.);
}

sub execute {
    my $self = shift;

    my $anno_file   = $self->anno_file;
    my $vcf_file    = $self->vcf_file;
    my $output_file = $self->output_file;
    if (-e $output_file) {
        die "output_file already exists: $output_file\n";
    }

    my $RSid = store_RSid($vcf_file);

    my $out_fh  = Genome::Sys->open_file_for_writing($output_file);
    my $anno_fh = Genome::Sys->open_file_for_reading($anno_file);
    while (my $line = $anno_fh->getline) {
        chomp $line;

        next if $line =~ /^chromosome/;

        my @list = split(/\t/, $line);
        my ($chr, $pos, $RSid_var_allele) = (split(/\t/, $line))[0, 1, 4];

        my $key = RSid_key($chr, $pos);
        my $RSid_info = $RSid->{$key}->{$RSid_var_allele};
        my $rsID = $RSid_info->{rsID} || '-';
        my $GMAF = $RSid_info->{GMAF} || '-';
        push @list, $rsID, $GMAF;

        my $str = join("\t", @list);
        $out_fh->print($str, "\n");
    }
    $anno_fh->close;
    $out_fh->close;

    return 1;
}

sub store_RSid {
    my $vcf = shift;

    my $RSid = {};
    my $vcf_fh = Genome::Sys->open_gzip_file_for_reading($vcf);
    while (my $line = $vcf_fh->getline) {
        chomp $line;

        next if $line =~ /^\#/;

        my ($chr, $pos, $rsID, $ref, $var, $INFO) = (split(/\t/, $line))[0..4, 7];

        next if $rsID eq '.'; # skip if RSid not defined

        my $GMAF = ($INFO =~ /(GMAF=[0-9.]+)/)[0] || '-';
        my $key = RSid_key($chr, $pos);

        # there is some assumed symmetry between $var and $dbSNPinfo
        my @var_alleles = split(/, /, $var);
        my @dbSNPids = split_dbSNPBuildID($INFO);

        for (my $i = 0; $i < @dbSNPids; $i++) {
            next unless $dbSNPids[$i] =~ /^\d+$/;

            my $RSid_var_allele = $var_alleles[$i];
            $RSid->{$key}->{$RSid_var_allele}->{rsID} = $rsID;
            $RSid->{$key}->{$RSid_var_allele}->{GMAF} = $GMAF
        }
    }
    $vcf_fh->close();
    return $RSid;
}

sub split_dbSNPBuildID {
    my $INFO = shift;
    die 'split_dbSNPBuildID without argument' unless defined $INFO;
    my $dbSNPinfo = ($INFO =~ /dbSNPBuildID=([0-9, .]+)/)[0];
    return unless defined $dbSNPinfo;
    my @dbSNPids = split(/,\s*/, $dbSNPinfo);
    return @dbSNPids;
}

sub RSid_key {
    my ($chr, $pos) = @_;
    return join('_', $chr, $pos);
}

1;
