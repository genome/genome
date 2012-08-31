package Genome::Model::Tools::Array::ConvertHapMapGenotypeToGoldSnp;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bio::DB::Fasta;

class Genome::Model::Tools::Array::ConvertHapMapGenotypeToGoldSnp {
    is => 'Command',
    has => [
    genotype_file => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "Input file of Hapmap Genotype data for a single individual",
    },
    output_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "Output file name for converted Genotype Data. Will be in a format that mimics a Gold SNP file",
    },        
    REFDIR =>
    {
        #This is the path to the reference sequence used for aligning the model
        type => 'String',
        is_optional => 0,
        default => "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c",
    },        
    refdb =>
    {
        type => 'Reference',
        is_optional => 1,
    },
    ],
};


sub execute {
    my $self=shift;

    #Check on the file names
    unless(-f $self->genotype_file) {
        $self->error_message("Hapmap genotype file is not a file: " . $self->genotype_file);
        return;
    }

    #Check and open filehandles
    my $genotype_fh=IO::File->new($self->genotype_file);
    unless($genotype_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->genotype_file );
        return;
    }

    my $output_fh=IO::File->new($self->output_file,"w");
    unless($output_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->output_file );
        return;
    }

    $self->refdb(Bio::DB::Fasta->new($self->REFDIR));

    while(my $line = $genotype_fh->getline) {
        chomp $line;
        next if $line =~ /chrom/;
        my @fields = $self->convert($line);
        next if(grep {$_ eq 'N'} @fields[(3,4)]);
        print $output_fh (join "\t", @fields,"\n");
    }

    return 1;
}

1;

sub help_detail {
    return "This module takes a Hapmap Genotype file for a single patient and converts it to Gold SNP style format for use with established tools";
}

sub help_brief {
    return "Convert a Hapmap Genotype file to Gold SNP style";
}


#Create hashes of gold SNPs

sub convert {
    my ($self,$line) = @_;
    my $refdb = $self->refdb;
    my ($snp_id,$dbsnp_alleles,$chr,$pos,$strand,$genotype) = split /\s+/, $line; 

    #adjust chromosome
    $chr =~ s/^chr(.*)/$1/i;

    my @alleles = split //, $genotype;

    my $ref_allele = $refdb->seq($chr,$pos => $pos);
    $ref_allele = uc($ref_allele);

    my @calls;
    my $allele_num;
    
    for($allele_num = 0; $allele_num < 2; $allele_num++) {
        $calls[$allele_num] = $alleles[$allele_num] eq $ref_allele ? 'ref' : 'SNP';
    }

    return ($chr, $pos, $pos, @alleles[(0,1)],@calls[(0,1,0,1)]); 
}

