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


    my $annotation = store_annotation($anno_file);

    #use a temp file, then sort it into the output file at the end
    my ($tmp_fh,$tmpfile) = Genome::Sys->create_temp_file;
    unless($tmp_fh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }

    #output the header line, if it exists;
    my $out_fh = Genome::Sys->open_file_for_writing($output_file);
    if($annotation->{"header"}->{"header"}){
        print $out_fh $annotation->{"header"}->{"header"} . "\trsid\tGMAF\n";
        $annotation->{"header"}->{"header"} = 0;
    }
    $out_fh->close;


    #read through the VCF, outputting any annotation lines that match
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

        if($rsID eq '.'){
            $rsID = "-";
        }

        my $GMAF = ($INFO =~ /(GMAF=[0-9.]+)/)[0] || '-';

        my $key = RSid_key($chr, $pos);

        if(defined( $annotation->{$key}) ){
            # there is some assumed symmetry between $var and $dbSNPinfo
            my @var_alleles = split(/, /, $var);
            my @dbSNPids = split_dbSNPBuildID($INFO);

            for (my $i = 0; $i < @dbSNPids; $i++) {
                next unless $dbSNPids[$i] =~ /^\d+$/;

                my $RSid_var_allele = $var_alleles[$i];


                if(defined($annotation->{$key}->{$RSid_var_allele})){
                    if($annotation->{$key}->{$RSid_var_allele} ne "0"){                    
                        print $tmp_fh $annotation->{$key}->{$RSid_var_allele} . "\t$rsID\t$GMAF\n";
                        $annotation->{$key}->{$RSid_var_allele} = 0;
                    }
                }
            }
        }
    }
    $vcf_fh->close();

    #take care of any positions/alleles not in the VCF
    for my $k (keys(%{$annotation})){
        for my $v (keys(%{$annotation->{$k}})){
            if($annotation->{$k}->{$v}){
                print $tmp_fh $annotation->{$k}->{$v} . "\t-\t-\n";
            }
        }
    }
    $tmp_fh->close;

    #should probably do a proper shell out wrapper here...
    `joinx sort -i $tmpfile >>$output_file;`;

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
        my ($chr, $pos, $RSid_var_allele) = (split(/\t/, $line))[0, 1, 4];
        my $key = RSid_key($chr, $pos, );

        $annotation->{$key}->{$RSid_var_allele} = $line;
    }
    $anno_fh->close;
    return $annotation;
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
    return join('_', ($chr, $pos));
}

1;
