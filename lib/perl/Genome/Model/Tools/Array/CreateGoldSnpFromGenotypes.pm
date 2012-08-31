package Genome::Model::Tools::Array::CreateGoldSnpFromGenotypes;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bio::DB::Fasta;

class Genome::Model::Tools::Array::CreateGoldSnpFromGenotypes {
    is => 'Command',
    has => [
    genotype_file1 => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "Input file of Genotypes from one platform",
    },
    genotype_file2 => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "Input file of Genotypes from additional platform",
    },
    output_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "a Gold SNP file",
    },        
    reference_fasta_file => {
        type => 'String',
        doc => 'The path to the reference sequence fasta file to use',
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
    unless(-f $self->genotype_file1) {
        $self->error_message("First Genotype file is not a file: " . $self->genotype_file1);
        return;
    }

    #Check and open filehandles
    my $genotype_fh1=IO::File->new($self->genotype_file1);
    unless($genotype_fh1) {
        $self->error_message("Failed to open filehandle for: " .  $self->genotype_file1 );
        return;
    }
    
    #Check on the file names
    unless(-f $self->genotype_file2) {
        $self->error_message("Second Genotype file is not a file: " . $self->genotype_file2);
        return;
    }

    #Check and open filehandles
    my $genotype_fh2=IO::File->new($self->genotype_file2);
    unless($genotype_fh2) {
        $self->error_message("Failed to open filehandle for: " .  $self->genotype_file2 );
        return;
    }

    my $output_fh=IO::File->new($self->output_file,"w");
    unless($output_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->output_file );
        return;
    }

    $self->refdb(Bio::DB::Fasta->new($self->reference_fasta_file));

    my ($chr2, $pos2, $genotype2) = (1,1,q{}); #expecting 

    while(my $line1 = $genotype_fh1->getline) {
        chomp $line1;

        my ($chr1, $pos1, $genotype1) = split /\s+/, $line1;

        my $line2;
        while((($pos1 > $pos2 && $chr1 eq $chr2) || ($pos1 < $pos2 && $chr1 ne $chr2))  && ($line2 = $genotype_fh2->getline)) {
            ($chr2, $pos2, $genotype2) = split /\s+/, $line2;
        }

        if($chr2 eq $chr1 && $pos1 == $pos2) {
            #intersecting position
            #check genotypes
            my @alleles = split //, uc($genotype1);
            my $ref = uc $self->refdb->seq($chr1, $pos1 => $pos1);

            if($genotype1 ne '--' && $genotype1 =~ /[ACTGN]/ && $genotype2 =~ /[ACTGN]/ && ($genotype1 eq $genotype2 || $genotype1 eq reverse($genotype2))) {

                #print genotypes with call
                my $type1 = ($alleles[0] eq $ref) ? 'ref' : 'SNP';
                my $type2 = ($alleles[1] eq $ref) ? 'ref' : 'SNP';
                print $output_fh "$chr1\t$pos1\t$pos1\t$alleles[0]\t$alleles[1]\t$type1\t$type2\t$type1\t$type2\n";
            }
        }

    }

    return 1;
}

1;

sub help_detail {
    return "This module takes two genotype files and converts the intersection of their genotypes and positions to Gold SNP file format";
}

sub help_brief {
    return "Take two genotype files and make a gold SNP file format";
}


