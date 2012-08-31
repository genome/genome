package Genome::Model::Tools::Somatic::FilterLoh;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::FilterLoh{
    is => 'Command',
    has => [
        tumor_snp_file => {
            is  => 'String',
            is_input  => 1,
            doc => 'The list of tumor SNPs in maq-like format (chromosome position reference variant)... often the somatic sniper output from the somatic pipeline',
        },
        normal_snp_file => {
            is  => 'String',
            is_input  => 1,
            doc => 'The list of normal SNPs in maq-like format...often the normal build filtered_snp_file',
        },
        output_file => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => "Variants that have successfully passed the LOH filter (are not LOH events)"
        },
        loh_output_file => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            doc => "Variants that have failed the LOH filter (are LOH events). (optional) If this is not provided, no output will be given for variants failing the filter."
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
    ],
};

sub help_brief {
    "Separate LOH calls from non-LOH calls",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic filter-loh --tumor tumor.snp --normal normal.snp -o non_loh.out
gmt somatic filter-loh --tumor tumor.snp --normal normal.snp --output non_loh.out --loh-output loh.out  
EOS
}

sub help_detail {                           
    return <<EOS 
This filters out SNVs that are likely to be the result of Loss of Heterozygosity (LOH) events. The Somatic Pipeline will pass these through on its own as they are tumor variants that differ from the normal. This script defines a variant as LOH if it is homozygous, there is a heterozygous SNP at the same position in the normal and the tumor allele is one of the two alleles in the normal.  
EOS
}

sub execute {
    my $self = shift;

    unless(-f $self->tumor_snp_file) {
        $self->error_message($self->tumor_snp_file . " is not a file");
        die;
    }

    unless(-f $self->normal_snp_file) {
        $self->error_message($self->normal_snp_file . " is not a file");
        die;
    }

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $out_fh = IO::File->new(">".$self->output_file);
    unless($out_fh) {
        $self->error_message("Unable to open " . $self->output_file . " for writing");
        die;
    }

    my $loh_fh;
    if (defined $self->loh_output_file) {
        $loh_fh = IO::File->new(">".$self->loh_output_file);
        unless($loh_fh) {
            $self->error_message("Unable to open " . $self->loh_output_file . " for writing");
            die;
        }
    }

    #MAKE A HASH OF NORMAL SNPS!!!!!!!!!!!!!
    #Assuming that we will generally be doing this on small enough files (I hope). I suck.

    my $normal_snp_fh = IO::File->new($self->normal_snp_file,"r");
    unless($normal_snp_fh) {
        $self->error_message("Unable to open " . $self->normal_snp_file);
        die;
    }

    my %normal_variants;

    while(my $line = $normal_snp_fh->getline) {
        chomp $line;
        my ($chr, $pos, $ref, $var_iub) = split /\t/, $line;
        #first find all heterozygous sites in normal
        next if($var_iub =~ /[ACTG]/);
        my @alleles = Genome::Info::IUB->iub_to_alleles($var_iub);
        $normal_variants{$chr}{$pos} = join '',@alleles;
    }

    $normal_snp_fh->close;
    
    my $tumor_snp_fh = IO::File->new($self->tumor_snp_file,"r");
    unless($tumor_snp_fh) {
        $self->error_message("Unable to open " . $self->tumor_snp_file);
        die;
    }

    while(my $line = $tumor_snp_fh->getline) {
        chomp $line;

        # Getting adapted intput here so 
        my ($chr, $pos, $ref, $var_iub) = split /\t/, $line;
        
        #now compare to homozygous sites in the tumor
        if($var_iub =~ /[ACTG]/ && exists($normal_variants{$chr}{$pos})) {
            if(index($normal_variants{$chr}{$pos},$var_iub) > -1) {
                #then they share this allele and it is LOH
                if (defined $loh_fh) {
                    $loh_fh->print("$line\n");
                }
            }
            else {
                $out_fh->print("$line\n");
            }
        }
        else {
            $out_fh->print("$line\n");
        }
    }

    $tumor_snp_fh->close;

    return 1;
}

1;
