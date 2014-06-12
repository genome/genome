package Genome::Model::Tools::Analysis::RemoveContaminatingVariants;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw(firstidx);
use File::Basename;

class Genome::Model::Tools::Analysis::RemoveContaminatingVariants {
    is => 'Command',
    has_input => [
        input_file => {
            doc => 'The path to the input file. This file is assumed to have a header',
        },
        output_file => {
            doc => 'The path to the output file containing only passing variants.',
        },
        reference_read_col => {
            doc => 'The position of the column containing reference supporting read counts. If unspecified, tries to infer this by looking for "Tumor_ref_count" in header',
            is_optional => 1,
        },
        variant_read_col => {
            doc => 'The position of the column containing variant supporting read counts. If unspecified, tries to infer this by looking for "Tumor_var_count" in header',
            is_optional => 1,
        },
        fdr_cutoff => {
            doc => 'Sites with corrected p-values above this range will be removed',
            is_optional => 1,
            default => 0.05,
        }

    ],
};

sub help_synopsis {
    q(gmt analysis remove-contaminating-variants --input-file snvs.indels.annotated --output-file snvs.indels.annotated.filtered);
}

sub help_detail {
    "This tool uses a fisher's exact test to determine whether a variant's read counts are significantly different than one of equal depth that is purely heterozygous (50% variant allele frequency or homozygous (100% variant allele frequency). This is useful for removing host contamination from xeno/allografts or tumors grown in culture";
}

sub execute {
    my $self = shift;
    my $input_file   = $self->input_file;

    my $refcol = 0;
    my $varcol = 0;
    if($self->reference_read_col){
        $refcol = $self->reference_read_col;
    }
    if($self->variant_read_col){
        $varcol = $self->variant_read_col;
    }


    #read the file's header
    open my $file, '<', $input_file; 
    my $line = <$file>; 
    chomp($line);
    my @header = split("\t",$line);
    close $file;
    #if ref and var cols not defined, infer them from the header
    unless($refcol){
        $refcol = firstidx{ $_ eq "Tumor_ref_count" } @header;
        $varcol = firstidx{ $_ eq "Tumor_var_count" } @header;
        #increment them, since R's arrays are one-based
        $refcol++;
        $varcol++;
    }
    
    if(!($refcol) || !($varcol)){
        die("Couldn't infer reference and variant count columns from the header. Please provide them as arguments");
    }


    #create temp directory for munging file
    my $tempdir = Genome::Sys->create_temp_directory();
    unless($tempdir) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }

    my $dirname = dirname(__FILE__);
    my $rfile = "$dirname/RemoveContaminatingVariants.R";

    #add pvals as trailing columns
    my $cmd = "Rscript $rfile $input_file $refcol $varcol $tempdir/pvals";
    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
        );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }


    #now, print all sites below fdr threshold
    my $inFh = IO::File->new( "$tempdir/pvals" ) || die "can't open file\n";
    my $outfile = Genome::Sys->open_file_for_writing($self->output_file);

    #print the header
    print $outfile join("\t",(@header,"type","pval","qval")) . "\n";

    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @F = split("\t",$line);
        #last column is FDR value
        unless ($F[$#F] > $self->fdr_cutoff){
            print $outfile $line . "\n";
        }
    }
    close($inFh);

    1;
}
