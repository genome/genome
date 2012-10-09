package Genome::Model::Tools::Vcf::RareDelTest;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationRate - Calculate the mutation rate (per megabase) given a list of mutations (e.g. tier1 SNVs) and a set of regions (e.g. coding space)
#
#    AUTHOR:        Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#    CREATED:    04/22/2011 by D.K.
#    MODIFIED:    04/22/2011 by D.K.
#
#    NOTES:
#
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it


## Pre-define a ranking system for VEP annotation, where higher = more severe ##


## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Vcf::RareDelTest {
    is => 'Command::V2',

    has_input => [                                # specify the command's single-value properties (parameters) <---
        gene_file    => { is => 'Text', doc => "Input rare-deleterious gene table" , is_optional => 0},
        output_file    => { is => 'Text', doc => "Output file for genes with FET p-value" , is_optional => 0, is_output => 1},
        output_significant    => { is => 'Text', doc => "Output file for genes with significant FET p-value" , is_optional => 0, is_output => 1},
        output_details    => { is => 'Text', doc => "Output file for details of significant genes" , is_optional => 1},
        p_value_threshold    => { is => 'Text', doc => "Default p-value threshold to report details for genes" , is_optional => 0, default => 0.05},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs a Fisher's Exact Test on a rare deleterious table of genes"
}

sub help_synopsis {
    return <<EOS
This command produces a rare/deleterious table by gene from VCF and VEP files
EXAMPLE:    gmt vcf rare-del-test --gene-file my.deltable.tsv --variant-file my.deltable.variants.tsv --output-file my.deltable.fisher.tsv
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    my $gene_file = $self->gene_file;
    my $output_file = $self->output_file;
    my $output_significant = $self->output_significant;
    my $p_threshold = $self->p_value_threshold;

    my $R_code = <<EOS;
# Load data from gmt vcf rare-del-table
x <- read.table("$gene_file", header=TRUE);

# function to compute Fisher's exact test p-value given a row of the rare-del-table
test.row <- function(x) {
    v <- as.numeric(x[c("controls_without_var", "controls_with_var", "cases_without_var", "cases_with_var")]);
    m <- matrix(v, nrow=2);
    fisher.test(m)\$p.value;
}

fet_p_value <- apply(x, 1, test.row);

# Append p-value column, find significant entries, and write result files
y <- cbind(x, fet_p_value);
signif <- y[fet_p_value < $p_threshold,];
write.table(y, "$output_file", sep="\t", quote=FALSE, row.names=FALSE);
write.table(signif, "$output_significant", sep="\t", quote=FALSE, row.names=FALSE);
EOS

    my ($fh, $path) = Genome::Sys->create_temp_file;
    $fh->print($R_code);
    $fh->close();
    Genome::Sys->shellcmd(
        cmd => "Rscript $path",
        input_files => [$gene_file],
        output_files => [$output_file, $output_significant],
    );

    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;
