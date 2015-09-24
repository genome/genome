package Genome::Model::Tools::EpitopePrediction::BindingFilter;
# Based on myParseReportsTopGenesTopAlleleStringent_4.pl by TWYLIE (Feb 2012)

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);

use Genome::Utility::IO::SeparatedValueReader qw();
use Genome::Utility::IO::SeparatedValueWriter qw();

class Genome::Model::Tools::EpitopePrediction::BindingFilter {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    has_input => [
        fof_file => {
            is  => 'FilePath',
            doc => 'FOF containing list of parsed epitope files for different allele-length combinations (same sample)',
        },
    ],
    has_optional_input => [
        variant_file => {
            is  => 'FilePath',
            doc => 'Original variant file in TGI annotation format. If provided, this will create a second output file that appends binding prediction information to the variant annotation lines and outputs those that pass filters',
        },
        minimum_fold_change => {
            is => 'Number',
            doc => 'Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT)',
            default => 0,
        },
        binding_threshold => {
            is => 'Number',
            doc => 'report only epitopes where the mutant allele has ic50 binding scores below this value',
            default => 500,
        },

    ],
    has_output => [
        output_file => {
            is  => 'FilePath',
            doc => 'Output .xls file containing list of filtered epitopes based on binding affinity for each allele-length combination per gene',
        }
    ],
};

sub help_brief {
    "Takes in a FOF with parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinity."
}

sub execute {
    my $self = shift;


    # if provided, grab the variants in TGI format and hash them by gene and mutation
    my %variants;
    if(defined($self->variant_file)){
        my $var_fh = Genome::Sys->open_file_for_reading($self->variant_file);
        while (my $line = $var_fh->getline) {
            chomp $line;
            #header?
            if($line=~/chromosome_name/){
                $variants{"header"} = $line;
                next;
            }
            my @F = split("\t",$line);            
            my $gene = $F[6];
            my $amino_acid_change = $F[15];
            $amino_acid_change =~ s/p\.//; #change p.E435G to E435G like the netMHC output
            $variants{join("\t",($gene,$amino_acid_change))} = $line;
        }
    }



    # read in the netMHC predictions and filter
    my $fof_fh = Genome::Sys->open_file_for_reading($self->fof_file);
    my $out_fh = Genome::Utility::IO::SeparatedValueWriter->create(
        output    => $self->output_file,
        separator => "\t",
        headers   => [
            'Mode',
            'Sample',
            'Length',
            'Gene Name',
            'Allele',
            'Point Mutation',
            'Sub Peptide Position',
            'MT Score',
            'WT Score',
            'MT Epitope Seq',
            'WT Epitope Seq',
            'Fold Change',
        ],
    ) or die "Unable to create SeparatedValueWriter\n";

    my %prediction;
    my $threshold = $self->binding_threshold;
    while (my $file = $fof_fh->getline) {
        chomp $file;
        my $basename = basename( $file);
        my @f      = split( /\./, $basename );
        my $sample = $f[0];
        $sample =~ s/_netmhc//g;
        my $allele = $f[1];
        my $length = $f[2];

        my $mode = 'filtered';
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            input     => $file,
            separator => "\t",
            headers   => [
                'Gene Name',
                'Point Mutation',
                'Sub Peptide Position',
                'MT Score',
                'WT Score',
                'MT Epitope Seq',
                'WT Epitope Seq',
                'Fold Change',
            ],
        );
        $reader->next; # skip headers
        while (my $gene_name = $reader->next) {
            $gene_name->{Allele} = $allele;
            push( @{ $prediction{$mode}->{$sample}->{$length}->{genes} }, $gene_name );
        }
    }
    close ($fof_fh);

    my %best;
    foreach my $mode (sort keys %prediction) {
        foreach my $sample (sort keys %{ $prediction{$mode} }) {
            foreach my $length (sort keys %{ $prediction{$mode}->{$sample} }) {
                foreach my $gene (sort @{ $prediction{$mode}->{$sample}->{$length}->{genes} }) {
# BEST
                    if ( $best{$sample}->{$gene->{'Gene Name'}}->{SCORE} ) {
                        if ($gene->{'MT Score'} > $best{$sample}->{$gene->{'Gene Name'}}->{SCORE}) {
                            next;
                        }
                        if ($gene->{'MT Score'} < $best{$sample}->{$gene->{'Gene Name'}}->{SCORE}) {
                            $best{$sample}->{$gene->{'Gene Name'}}->{GENES} = [];
                        }
                    }

                    $best{$sample}->{$gene->{'Gene Name'}}->{SCORE} = $gene->{'MT Score'};
                    $gene->{Sample} = $sample;
                    $gene->{Length} = $length;
                    $gene->{Mode}   = $mode;
                    push( @{ $best{$sample}->{$gene->{'Gene Name'}}->{GENES} }, $gene );
                }
            }
        }
    }

# REPORTING

    my $varfile;
    if(defined($self->variant_file)){
        my $filename = $self->output_file . ".varanno";    
        $varfile = Genome::Sys->open_file_for_writing($filename);

        #print header for variant output, if needed        
        if(defined($variants{"header"})){
            print $varfile join("\t",($variants{"header"},
                                     "GeneName","MutationEffect","HLAallele","PeptideLength",
                                     'SubPeptidePosition','MT Score',
                                     'WT Score','MT Epitope Seq',
                                     'WT Epitope Seq','Fold Change')) . "\n";
        }
    }
                                 
    
    foreach my $sample (sort keys %best) {
        foreach my $gene (sort keys %{ $best{$sample} }) {
            foreach my $entry (@{ $best{$sample}->{$gene}->{GENES} }) {
                if (($entry->{'MT Score'} < $threshold) && ($entry->{'Fold Change'} > $self->minimum_fold_change)){
                    #write out epiptopes
                    $out_fh->write_one($entry);

                    #also write out variants with epitopes appended
                    if(defined($self->variant_file)){
                        my $key = join("\t",$gene,$entry->{'Point Mutation'});
                        if(defined($variants{$key})){                            
                            print $varfile join("\t", ($variants{$key},
                                                      $gene,$entry->{'Point Mutation'},$entry->{Allele},
                                                      $entry->{Length},$entry->{'Sub Peptide Position'},$entry->{'MT Score'},
                                                      $entry->{'WT Score'},$entry->{'MT Epitope Seq'},
                                                      $entry->{'WT Epitope Seq'},$entry->{'Fold Change'})) . "\n";

                        } else {
                            
                            $self->warning_message("Couldn't find variant for " . $gene . " " . $entry->{'Point Mutation'} . "in variant_file");
                        }
                    }
                }
            }
        }
    }

    return 1;
}
