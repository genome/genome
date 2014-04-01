package Genome::Model::Tools::Music::Cosmic;

use warnings;
use strict;
use Genome;
use IO::File;
use FileHandle;
use Text::CSV_XS;
# default cutoff for GDSC
use constant PVALUE => 0.0137033;
use constant IC50   => 1.0;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Cosmic {
    is => 'Command::V2',
    has_input => [
        var_file    => { is => 'String', doc => "List of variants using TCGA MAF specification v2.3, VCF v4.0/v4.1 or WU (Washington University) annotation" },
        output_file => { is => 'String', doc => "Output MAF, VCF or WU annotation file with additional data" },
        ref_build   => { is => 'Text', example_values => ['Build36','Build37'], valid_values => ['Build36', 'Build37'], doc => "Specify the genome build used by genomic loci in the provided variant list" },
    ],
    has_optional_input => [
        gdsc_dir       => { is => 'Path', doc => "Folder contains of the association information of drugs and targets" },
        cosmic_dir     => { is => 'Path', doc => 'Cosmic mutation database folder' },
        omimaa_dir     => { is => 'Path', doc => 'omim amino acid mutation database folder' },
        dgene_dir      => { is => 'Path', doc => 'dGene database folder' },
        wu_annotation  => { is => 'Boolean', default => 0, doc => "Use this if input MAF contains WUSTL annotation format headers" },
        vcf_annotation => { is => 'Boolean', default => 0, doc => "Use this if input file is VCF v4.0/v4.1 format" },
        nuc_range      => { is => 'Integer', default => '5', doc => "Set how close a 'near' match is when searching for nucleotide position near matches" },
        aa_range       => { is => 'Integer', default => 2, doc => "Set how close a 'near' match is when searching for amino acid near matches" },
        strand_option  => { is => 'Boolean', default => 0, doc => "Set optional strand for Cosmic amino acid when searching for amino acid matches, default is strand sensitive" },
        p_value        => { is => 'Number', default => '0.0137033', doc => "Set a significance value cutoff from GDSC Genomic correlations with MANOVA, default is 20% FDR" },
        ic_50          => { is => 'Number', default => '1.0', doc => "Set a IC50 Effect cutoff from GDSC Genomic correlations with MANOVA, default is 1.0" },
    ],
    doc => "Match a list of variants to those in COSMIC, and highlight druggable targets",
};

sub help_detail {
    return <<HELP;
This module compares variants in the MAF, WU annotation or VCF format file user submitted 
to known variants in COSMIC V61, OMIM (2012.09.20) and GDSC V1.1 databases with 
six conditions following:

Position match:

1. Identify variants that exactly match a known variant in COSMIC for both nucleotide AND Amino Acid (when AA is available);
2. Identify variants that match the same genomic locus but report a different nucleotide change;
3. Identify variants that match the same Amino Acid but report a different Amino Acid change;

Proximity match:

4. Identify variants that are within N bps of distance to a known variant in COSMIC (default N = 5);
5. Identify variants that are within M AA of distance to a known variant in COSMIC and OMIM (default M = 2);

This module also check drugable genes using dGgene database and the drugs information of 
annotated genes from MAF file, using gene & drug mapping data from GDSC database. For each 
gene that a variant is annotated to, the drugs that can target it and the drugs based on 
given p-value and IC_50 value from Multivariate ANOVA for all compounds will be reported. 
Any site without a match in a particular databases is reported as "NA" with respect to 
that database.

The output of this module returns each row the original input MAF and WU annotation file 
with 4 columns appended to the end of each, the first two columns for COSMIC matches and 
one column for each of the OMIM and GDSC databases. For VCF format, the output information 
will be added into "INFO" part of VCF file. 

The detailed appended columns format for MAF and WU annotaion is described in the following:

------------------------------------------

COSMIC_Nuc/COSMIC_AA/OMIM   GDSC  dGene  

exact|position|proximity  Target_drugs|drugs_based_on_cutoff    dGene_description

** 
Exact, position and proximity mean different match types, each separated by a "|" symbol. 
For "COSMIC_Nuc" column, the number of samples and tissue types will be reported for each type of match.
For "COSMIC_AA" and "OMIM" columns, only amino acid change will be reported for each type of match.
For "GDSC" column, list drugs will be reported for each type of match, each separated by a "|" symbol.
For "dGene" column, the matched drugable gene description will be reported.

The detailed appended "INFO" format for VCF format will be shown in the outputfile header.
By preprocessing of lift-over, this module compares both build 36 and build 37 coordinates 
that are specified in COSMIC to the coordinates in your variants file.

In addition to the standard version 2.3 MAF headers, there needs to be 3 columns appended. 
Only nucleotide matches from COSMIC database will be reported, if no appended columns provided. 
These column headers in the MAF must have these names in the header in order for the tool to find them:

   transcript_name - the transcript name, such as NM_000028
 amino_acid_change - the amino acid change, such as p.R290H
           strand  - the strand of transcript, such as -1/+1
HELP
}

sub help_synopsis {
    return <<EOS;
 ... music cosmic \\
        --var-file input_dir/myMAF.tsv \\
        --output-file output_dir/myMAF_output.tsv \\

 ... music cosmic \\
        --var-file input_dir/myMAF.tsv \\
        --output-file output_dir/myMAF_output.tsv \\
        --omimaa-dir omim_dir/ \\
        --cosmic-dir cosmic_dir/ \\
EOS
}
sub _doc_credits {
    return <<EOS

This tool depends on copies of data from the following databases, packaged in a form useable for quick analysis:

 * COSMIC - http://www.sanger.ac.uk/genetics/CGP/cosmic/
 * OMIM - http://www.ncbi.nlm.nih.gov/omim
 * GDSC - http://www.cancerrxgene.org 

EOS
}
sub _doc_authors {
    return <<EOS
        Beifang Niu, Ph.D.
EOS
}
sub doc_copyright_years {
    my $year = (localtime(time))[5] + 1900;
    return ($year);
}

#########################################################################

sub execute {
    my $self = shift;
    ## Based on preprocessing of COSMIC db
    my $wu_anno      = $self->wu_annotation;
    my $maf_file     = $self->var_file;
    my $omim_dir     = $self->omimaa_dir;
    my $aa_range     = $self->aa_range;
    my $vcf_anno     = $self->vcf_annotation; 
    my $nuc_range    = $self->nuc_range;
    my $strand_op    = $self->strand_option;
    my $ref_build    = $self->ref_build;
    my $mano_ic_50   = $self->ic_50;
    my $cosmic_dir   = $self->cosmic_dir;
    my $output_file  = $self->output_file;
    my $mano_p_value = $self->p_value;
    # GDSC MANOVA results
    my $gdsc_dir     = $self->gdsc_dir;
    my $gdsc_file;
    my $omim_file;
    my $cosmic_mu;
    my $gdsc_manova_file;
    # dGene results
    my $dgene_dir    = $self->dgene_dir;
    my $dgene_file;
    # Default pvalue and ic50 cutoff
    my $default_cutoff = 1;
    my $default_ic50   = IC50;
    my $default_pvalue = PVALUE;
    $default_cutoff = 0 if ($mano_p_value != PVALUE or $mano_ic_50 != IC50);
    ## Preprocessing
    # Para check
    unless($ref_build =~ m/build36/i xor $ref_build =~ m/build37/i) {
        $self->error_message("You must either specify reference_build as either \"Build36\" or \"Build37\"");
        die $self->error_message;
    }
    unless($nuc_range > 0 && $nuc_range <= 100) {
        $self->error_message("You must specify nucletide range N bps (0< N <=100)");
        die $self->error_message;
    }
    # Input file check
    $omim_file = "$omim_dir/OMIM_aa_20120921.csv"              if (-d $omim_dir);
    $cosmic_mu = "$cosmic_dir/Cosmic.Parse.v60.var.v1"         if (-d $cosmic_dir);
    $gdsc_file = "$gdsc_dir/drugs_targets_associations.info"   if (-d $gdsc_dir);
    $gdsc_manova_file = "$gdsc_dir/gdsc_manova_output_w2.csv"  if (-d $gdsc_dir);
    $dgene_file = "$dgene_dir/dgene.tsv"                       if (-d $dgene_dir);
    ## Main proccessing
    # parse var file
    # maximal length of variant
    my $max_span = 0;
    my $build = "b36";
    $build = "b37" if ($ref_build =~ m/build37/i);
    my $fh = new FileHandle;
    die "Could not open COSMIC mutation file\n" unless($fh->open($cosmic_mu));
    $self->debug_message("loading and parsing COSMIC database ...\n");
    my $cosmics = ParseCosmicMutation($fh, $build, \$max_span, $strand_op);
    $fh->close;
    # parse GDSC target drugs file
    my $gfh = new FileHandle;
    die "Could not open GDSC drug targets file\n" unless($gfh->open($gdsc_file));
    $self->debug_message("loading and parsing GDSC drug response database ...\n");
    my $gene_drugs = ParseGdscFile($gfh);
    $gfh->close;
    # parse GDSC drug manova results file
    my $gfh_manova = new FileHandle;
    die "Could not open GDSC MANOVA results file\n" unless($gfh_manova->open($gdsc_manova_file));
    $self->debug_message("loading and parsing GDSC MANOVA results database ...\n");
    my $manova_pout = ParseGdscManovaFile($gfh_manova, $default_ic50, $default_pvalue);
    $gfh_manova->close;
    # parse dGene database file
    my $dgenefh = new FileHandle;
    die "Could not open dGene database file\n" unless($dgenefh->open($dgene_file));
    $self->debug_message("loading and parsing dGene database ...\n");
    my $dgene_pout = ParsedGeneFile($dgenefh);
    $dgenefh->close;
    # parse OMIM database
    my $omimaa = {};
    my $fh_omim = new FileHandle;
    die "Could not open OMIM file\n" unless($fh_omim->open($omim_file));
    $self->debug_message("loading and parsing GDSC MANOVA results database ...\n");
    ParseOMIMFile($fh_omim, $omimaa);
    $fh_omim->close;
    # parse maf file and find hits
    my @cosmic_ps = ($nuc_range, $max_span, $aa_range, $strand_op);
    my @manova_ps = ($mano_ic_50, $mano_p_value, $default_cutoff);
    my @hash_ps   = ($cosmics, $gene_drugs, $dgene_pout, $manova_pout, $omimaa, $wu_anno, $vcf_anno);
    my $mfh   = new FileHandle;
    my $outfh = new FileHandle;
    die "Could not open maf file\n"       unless($mfh->open($maf_file));
    die "Could not create output file\n"  unless($outfh->open(">$output_file"));
    $self->debug_message("Comaring to COSMIC database ...\n");
    # Find COSMIC/OMIM/GDSC matches
    FindMatches($mfh, $outfh, \@cosmic_ps, \@manova_ps, \@hash_ps);
    $mfh->close;
    $outfh->close;
    $self->debug_message("Processing Done. \n");
    return 1;
}

#############################################################

## Parse COSMIC mutation file
sub ParseCosmicMutation {
    my($cos_mu_file, $bd, $max_sp, $strandop) = @_;
    my $refc = {};
    $$max_sp = 0;
    while (my $ll = <$cos_mu_file>) {
        chomp($ll);
        my @cols = split(/\t/, $ll);
        my ($chr, $b36start, $b36stop, $b37start, $b37stop, $ref, $var, $sp_ts, $aach) = @cols;
        my $start = ($bd eq "b36" ? $b36start : $b37start);
        my $stop  = ($bd eq "b36" ? $b36stop  : $b37stop);
        my $span = $stop - $start;
        $$max_sp = $span if ($span > $$max_sp);
        # parse smaple id and tissue
        my @ta0 = map{[split/,/]} split /:/, $sp_ts;
        foreach my $tai0 (@ta0) {
            my($s_id, $tissue) = @$tai0;
            # Nucleotide exact or same position Cosmic match
            $refc->{nues}->{$chr}->{$start}->{$stop.$ref.$var}->{$s_id} = 1;
            $refc->{nuet}->{$chr}->{$start}->{$stop.$ref.$var}->{$tissue}++;
            # Nucleotide proximity
            $refc->{nups}->{$chr}->{$start}->{$stop}->{$s_id} = 1;
            $refc->{nupt}->{$chr}->{$start}->{$stop}->{$tissue}++;   
        }
        # parse transcript info
        my @ta1 = map{[split/,/]} split /:/, $aach;
        foreach my $tai1 (@ta1) {
            my($gene, $trans, $strand, $aa) = @$tai1;
            next if ($aa eq "NULL");
            my($res1, $rstart, $res2, $rstop, $new_res) = AA_Check($aa);
            if ($rstart) {
                my $strand_char = ($strand eq "+1" ? "+" : "-");
                unless($strandop){ $strand_char = "="; }
                my $t_aa = $res1.$rstart.$res2;
                # Amino Acid exact or same position Cosmic match
                $refc->{resp}->{$strand_char}->{$trans}->{$rstart}->{$res1.$rstop.$res2}->{$t_aa} = 1;
                # Amino Acid proximity
                $refc->{respp}->{$strand_char}->{$trans}->{$rstart}->{$t_aa} = 1;
            }
        }
    }
    return ($refc);
}
## Parse GDSC drug targets file
sub ParseGdscFile {
    my($g_file) = @_;
    # Drug targets parse results
    my $g_s = {};
    while (my $ll = <$g_file>) {
        chomp($ll);
        my @cols = split(/\t/, $ll);
        my ($drug, $gene) = @cols[0,1];
        $g_s->{$gene}->{$drug} = 1;
    }
    return ($g_s);
}

## Parse dGene database file
sub ParsedGeneFile {
    my($dg_file) = @_;
    # dGene parse results
    my $dg_s = {};
    while (my $ll = <$dg_file>) {
        next if ($ll =~ /^tax_id/);
        chomp($ll);
        my @cols = split(/\t/, $ll);
        my ($gene, $chr, $desc) = @cols[2,6,8];
        $dg_s->{$gene}->{$chr} = $desc;
    }
    return ($dg_s);
}

## Parse GDSC manovo results
sub ParseGdscManovaFile {
    my($g_file, $def_ic, $def_p) = @_;
    # MANOVA parse results
    my $m_pout = {};
    while (my $ll = <$g_file>) {
        next if ($ll =~ /^drug,gene,/);
        chomp($ll);
        my @t = split(/,/, $ll);
        my ($drug, $gene, $pvalue, $ic50) = @t[0,1,6,47];
        next if ($pvalue eq "NaN" or $ic50 eq "NaN");
        if ($pvalue <= $def_p and $ic50 <= $def_ic) { $m_pout->{def}->{$gene}->{$drug} = 1; }
        $m_pout->{dlist}->{$gene}->{$pvalue}->{$ic50}->{$drug} = 1;
    }
    my @t0 = keys %{$m_pout->{dlist}};
    foreach my $item (@t0) {
        my @t1 = keys %{$m_pout->{dlist}->{$item}};
           @t1 = sort { $a <=> $b } @t1;
        $m_pout->{psort}->{$item} = \@t1;
    }
    return($m_pout);
}
## Parse OMIM database
sub ParseOMIMFile {
    my($omim_fh, $omim_aa) = @_;
    while (my $ll = <$omim_fh>) {
        next if ($ll =~ /^gene/);
        chomp($ll);
        my @t = split(/\t/, $ll);
        # maybe we need report disorder id if there is
        my ($gene, $pos, $aao, $aamu) = @t[0,2,3,4];
        $omim_aa->{$gene}->{$pos}->{$aao.$pos.$aamu} = 1;
    }
}
## Match hits finding
sub FindMatches {
    my($smfh, $outfh, $cosmic_ps, $m_ps, $h_ps) = @_;
    # cosmic and omim paras
    my($dis, $max_sp, $aa_ran, $strand_op) = @$cosmic_ps;
    my($cosmics, $g_d, $dg_pout, $m_pout, $omim, $wu_anno, $vcf_anno) = @$h_ps;
    my($def_ic, $def_p, $def_or) = @$m_ps;
    my($total, $header, $rcols, $info_col) = FileHeadParse($smfh, $wu_anno, $vcf_anno);
    #print $header;
    print $outfh $header;
    # content parse
    my $hited = {};
    while (my $ll = <$smfh>) {
        # This part should be improve based one different file format 
        chomp($ll);
        # VCF format 
        if ($vcf_anno) {
           my ($nucinfo, $resinfo) = ParseCols($ll, $total, $rcols, $wu_anno, $vcf_anno);
           #nucleotide match
           my $vcf_nuc_cosmic_results = "";
           my $vcf_aa_cosmic_results  = "";
           my $vcf_gdsc_results = "";
           my $vcf_dgene_results = "";
           foreach my $nuc_info_line (@{$nucinfo}) {
               my ($mafchr, $mafstart, $mafstop, $mafref, $mafvar) = split(/\|/, $nuc_info_line);
               my @vcf_ps0 = ($mafchr, $mafstart, $mafstop, $mafref, $mafvar);
               my $vcf_nuc_results = COSMIC_Compare_nucleotide($cosmics, \@vcf_ps0, $dis, $max_sp);
               $vcf_nuc_cosmic_results .= $mafvar."|".$vcf_nuc_results.",";
           }
           chop($vcf_nuc_cosmic_results) if ($vcf_nuc_cosmic_results);
           ## Amino Acid match part
           my $strand_char = "+";
           unless($strand_op) { $strand_char = "="; }
           my %count_gene = ();
           foreach my $res_info_line (@{$resinfo}) {
               my ($maftrans, $mafaa, $gene, $chr_dg) = split(/\|/, $res_info_line);
               my @vcf_ps_omim = ($gene, $maftrans, $mafaa, $aa_ran, $strand_char);
               my $vcf_cosmicaa_omim_results = COSMIC_Compare_aa($cosmics, $omim, \@vcf_ps_omim);
               $vcf_cosmicaa_omim_results =~s/\t/\|/;
               $vcf_aa_cosmic_results .= $maftrans."|".$gene."|".$vcf_cosmicaa_omim_results.",";
               next if (exists($count_gene{$gene}));
               $count_gene{$gene} = 1;
               my $vcf_tmp_targets = GDSCtargetDrugs($g_d, $gene);
               my $vcf_tmp_manova  = GDSCMANOVADrugs($m_pout, $gene, $def_or, $def_p, $def_ic);
               $vcf_gdsc_results .= $gene."|".$vcf_tmp_targets."|".$vcf_tmp_manova.",";
               my $vcf_tmp_dgene = dGeneCompare($dg_pout, $gene, $chr_dg);
               $vcf_dgene_results .= $gene."|".$vcf_tmp_dgene.",";
           }
           chop($vcf_aa_cosmic_results) if ($vcf_aa_cosmic_results);
           chop($vcf_gdsc_results)      if ($vcf_gdsc_results);
           chop($vcf_dgene_results)     if ($vcf_dgene_results);
           # output
           my @vcf_content = split(/\t/, $ll);
           my $iter = 0;
           my $test_cont = "";
           my $vcf_each_line = "";
           foreach my $vcf_item (@vcf_content) {
               if ($iter == $info_col) {
                   $vcf_item  .= "COSMIC_NUC=".$vcf_nuc_cosmic_results.";" if ($vcf_nuc_cosmic_results);
                   $vcf_item  .= "COSMIC_OMIM_AA=".$vcf_aa_cosmic_results.";" if ($vcf_aa_cosmic_results);
                   $vcf_item  .= "GDSC=".$vcf_gdsc_results.";" if ($vcf_gdsc_results);
                   $vcf_item  .= "DGENE=".$vcf_dgene_results.";" if ($vcf_dgene_results);
               }
               $vcf_each_line .= $vcf_item."\t";
               $iter++;
           }
           chop($vcf_each_line);
           $vcf_each_line .= "\n";
           print $outfh $vcf_each_line;
        }
        else { # MAF file or WU annotation
            my ($gene, $mafchr, $mafstart, $mafstop, 
                 $mafref, $mafvar, $maftrans, $mafstrand,
                 $mafaa) =  ParseCols($ll, $total, $rcols, $wu_anno, $vcf_anno);
            ## Nucletide match part
            my @maf_ps = ($mafchr, $mafstart, $mafstop, $mafref, $mafvar);
            my $cosmic_nuc_results = COSMIC_Compare_nucleotide($cosmics, \@maf_ps, $dis, $max_sp);
            ## Amino Acid match part
            my $strand_char = ($mafstrand eq "-1" ? "-" : "+");
            unless($strand_op) {
                $strand_char = "=";
            } 
            my @maf_ps_omim = ($gene, $maftrans, $mafaa, $aa_ran, $strand_char);
            my $cosmicaa_omim_results = COSMIC_Compare_aa($cosmics, $omim, \@maf_ps_omim);
            # GDSC target drugs
            my $target_drugs = GDSCtargetDrugs($g_d, $gene);
            my $gdsc_manova  = GDSCMANOVADrugs($m_pout, $gene, $def_or, $def_p, $def_ic);
            # GDSC MAVOVA results drugs
            my $gdsc_output = $target_drugs."|".$gdsc_manova;
            # dGene output
            my $dgene_output = dGeneCompare($dg_pout, $gene, $mafchr);
            my $results_output = $cosmic_nuc_results."\t".$cosmicaa_omim_results."\t".$gdsc_output."\t".$dgene_output;
            # load match hits information
            print $outfh $ll."\t".$results_output."\n";
        }
    }
}
## GDSC MANOVA
sub GDSCMANOVADrugs {
    my ($m_pout, $gene, $def_or, $def_p, $def_ic) = @_;
    # GDSC MAVOVA results drugs
    my $gdsc_manova_con = "";
    if ($def_or) {
        my $temp_con = "";
        my $d6 = exists($m_pout->{def}->{$gene});
        if ($d6) {
            my @t9 = keys %{$m_pout->{def}->{$gene}};
            foreach my $item9 (@t9) { $temp_con .= $item9.","; }
        }
        chop($temp_con) if ($temp_con);
        $gdsc_manova_con .= $temp_con if ($temp_con); 
        $gdsc_manova_con = "NA" if ($gdsc_manova_con eq "");
    }
    else {
        my $temp_con = "";
        my $d7 = exists($m_pout->{psort}->{$gene});
        if ($d7) {
            my @t10 = @{$m_pout->{psort}->{$gene}};
            foreach my $item10 (@t10) {
                next if ($item10 > $def_p);
                my @t11 = keys %{$m_pout->{dlist}->{$gene}->{$item10}};
                foreach my $item11 (@t11) {
                    next if ($item11 > $def_ic);
                    my @t12 = keys %{$m_pout->{dlist}->{$gene}->{$item10}->{$item11}};
                    foreach my $item12 (@t12) { $temp_con .= $item12.","; } 
                }
            }
        }
        chop($temp_con) if ($temp_con);
        $gdsc_manova_con .= $temp_con if ($temp_con); 
        $gdsc_manova_con = "NA" if ($gdsc_manova_con eq "");
    }
    # drug list
    return($gdsc_manova_con);
}
## GDSC target drugs
sub GDSCtargetDrugs {
    my ($gene_drugs, $gene) = @_;
    # GDSC target drugs
    my $gdsc_con = "";
    my $temp_con = "";
    my $d5 = exists($gene_drugs->{$gene});
    if ($d5) {
        my @t = keys %{$gene_drugs->{$gene}};
        foreach my $item (@t) { $temp_con .= $item.","; }
    }
    # cut last ch
    chop($temp_con) if ($temp_con);
    $gdsc_con .= $temp_con if ($temp_con); 
    $gdsc_con = "NA" if ($gdsc_con eq "");

    # results
    return ($gdsc_con);
}

# dGene match
sub dGeneCompare {
    my ($dgout, $gene, $chr) = @_;
    my $gdsc_con = "";
    my $d5 = exists($dgout->{$gene}->{$chr});
    if ($d5) { $gdsc_con = $dgout->{$gene}->{$chr}; }
    # cut last ch
    $gdsc_con = "NA" if ($gdsc_con eq "");
    return($gdsc_con);
}

## COSMIC nucle comparing
sub COSMIC_Compare_nucleotide {
    my ($cosmics, $maf_ps, $dis, $max_sp) = @_;
    my($mafchr, $mafstart, $mafstop, $mafref, $mafvar) = @$maf_ps; 
    # Nucleotide match infor container
    my $perfect_con     = "NA";
    my $samlocation_con = "NA";
    my $neighbor_con    = "NA";
    ## Nucletide match part
    # perfect match
    my $d1 = exists($cosmics->{nues}->{$mafchr}->{$mafstart}->{$mafstop.$mafref.$mafvar});
    if ($d1) {
        $perfect_con = "";
        # Allel both
        my $href1 = $cosmics->{nues}->{$mafchr}->{$mafstart}->{$mafstop.$mafref.$mafvar};
        my @t = keys %{$href1};
        my $samples = @t;
        $perfect_con .= $samples.":";
        my $href2 = $cosmics->{nuet}->{$mafchr}->{$mafstart}->{$mafstop.$mafref.$mafvar};
        my @t0 = keys %{$href2};
        foreach my $item (@t0) { $perfect_con .= $item." ".$href2->{$item}.","; }
        chop($perfect_con);
    }
    # same locations hits
    my %te_samples = ();
    my %te_tissue  = ();
    my $d2 = exists($cosmics->{nues}->{$mafchr}->{$mafstart});
    my $ex_item = $mafstop.$mafref.$mafvar;
    if ($d2) {
        my $href3 = $cosmics->{nues}->{$mafchr}->{$mafstart}; 
        my @t = keys %{$href3};
        foreach my $item (@t) {
            next if ($item eq $ex_item);
            my @t0 = keys %{$href3->{$item}};
            foreach my $item0 (@t0) { $te_samples{$item0} = 1; }
            my $href4 = $cosmics->{nuet}->{$mafchr}->{$mafstart}->{$item};
            my @t2 = keys %{$href4};
            foreach my $item2 (@t2) {
                $te_tissue{$item2} += $href4->{$item2};
            }
        }
        # write results
        my @t3 = keys %te_samples;
        my $samples = @t3;
        if ($samples) {
            $samlocation_con = "";
            $samlocation_con .= $samples.":";
            foreach my $item3 (keys %te_tissue) {
                $samlocation_con .= $item3." ".$te_tissue{$item3}.",";
            }
            chop($samlocation_con);
        }
    }
    # neighbor hits
    %te_samples = ();
    %te_tissue  = ();
    my $nstart    = $mafstart - $dis;
    my $det_start = $nstart   - $max_sp;
    my $nstop     = $mafstop  + $dis;
    $nstart    = 1 if ($nstart    <= 0);
    $det_start = 1 if ($det_start <= 0);
    foreach my $k ($det_start..$nstop) {
        # skip the first two steps
        next if ($k == $mafstart);
        my $d3 = exists($cosmics->{nups}->{$mafchr}->{$k});
        if ($d3) {
            my @t4 = keys %{$cosmics->{nups}->{$mafchr}->{$k}};
            foreach my $item5 (@t4) {
                next if ($item5 < $nstart);
                my @t5 = keys %{$cosmics->{nups}->{$mafchr}->{$k}->{$item5}};
                foreach my $item6 (@t5) {
                    $te_samples{$item6} = 1;
                }
                my @t6 = keys %{$cosmics->{nupt}->{$mafchr}->{$k}->{$item5}};
                foreach my $item7 (@t6) {
                    $te_tissue{$item7} += $cosmics->{nupt}->{$mafchr}->{$k}->{$item5}->{$item7};
                }
            } 
        }
    }
    # write results
    my @t7 = keys %te_samples;
    my $samples = @t7;
    if ($samples) {
        $neighbor_con = "";
        $neighbor_con .= $samples.":";
        foreach my $item8 (keys %te_tissue) {
            $neighbor_con .= $item8." ".$te_tissue{$item8}.",";
        }
        chop($neighbor_con);
    }
    # results
    my $cosmic_nuc_results = $perfect_con."|".$samlocation_con."|".$neighbor_con;
    return($cosmic_nuc_results);
}
## COSMIC database comparing (AA)
sub COSMIC_Compare_aa {
    my ($cosmics, $omim, $maf_ps) = @_;
    my($gene, $maftrans, $mafaa, $aa_range, $strand_char) = @$maf_ps;
    # Amino Acid change match infor container
    my $aa_perfect_con     = "NA";
    my $aa_samlocation_con = "NA";
    my $aa_neighbor_con    = "NA";
    # Grab aa information
    my($mafres1, $mafrstart, $mafres2, $mafrstop, $mafnew_res) = AA_Check($mafaa);
    unless($mafrstart) {
        return("NA|NA|NA\tNA|NA|NA");
    }
    # Perfect match
    my $ad1 = exists($cosmics->{resp}->{$strand_char}->{$maftrans}->{$mafrstart}->{$mafres1.$mafrstop.$mafres2});
    if ($ad1) {
        my $ahref0 = $cosmics->{resp}->{$strand_char}->{$maftrans}->{$mafrstart}->{$mafres1.$mafrstop.$mafres2};
        $aa_perfect_con = "";
        my @at0 = keys %{$ahref0};
        foreach my $aitem0 (@at0) { $aa_perfect_con .= $aitem0.","; }
        chop($aa_perfect_con) if ($aa_perfect_con);
        $aa_perfect_con = "NA" if ($aa_perfect_con eq "");
    }
    # Position match
    my %ate_aac = ();
    my $ad2 = exists($cosmics->{resp}->{$strand_char}->{$maftrans}->{$mafrstart});
    my $aex_item = $mafres1.$mafrstop.$mafres2;
    if ($ad2) {
        my $ahref1 = $cosmics->{resp}->{$strand_char}->{$maftrans}->{$mafrstart};
        my @at1 = keys %{$ahref1};
        foreach my $aitem1 (@at1) {
            next if ($aitem1 eq $aex_item);
            my @at2 = keys %{$ahref1->{$aitem1}};
            foreach my $aitem2 (@at2) { $ate_aac{$aitem2} = 1; }
        }
        #write results
        my @at3 = keys %ate_aac;
        if (@at3) {
            $aa_samlocation_con = "";
            foreach my $aitem3 (@at3) { $aa_samlocation_con .= $aitem3.","; }
            chop($aa_samlocation_con) if ($aa_samlocation_con);
        }
        $aa_samlocation_con = "NA" if ($aa_samlocation_con eq "");
    }
    # Proximity match
    my $iter_start = $mafrstart - $aa_range;
    my $iter_stop  = $mafrstart + $aa_range;
    $iter_start = 1 if ($iter_start <= 0);
    %ate_aac = ();
    foreach my $aitem4 ($iter_start..$iter_stop) {
        next if ($aitem4 == $mafrstart);
        my $ad3 = exists($cosmics->{respp}->{$strand_char}->{$maftrans}->{$aitem4});
        if ($ad3) {
            my @at4 = keys %{$cosmics->{respp}->{$strand_char}->{$maftrans}->{$aitem4}};
            foreach my $aitem5 (@at4) { $ate_aac{$aitem5} = 1; }
        }
    }
    # write results
    my @at5 = keys %ate_aac;
    if (@at5) {
        $aa_neighbor_con = "";
        foreach my $aitem6 (@at5) { $aa_neighbor_con .= $aitem6.","; }
        chop($aa_neighbor_con) if ($aa_neighbor_con);
    }
    $aa_neighbor_con = "NA" if ($aa_neighbor_con eq "");

    my $cosmic_aa_results = $aa_perfect_con."|".$aa_samlocation_con."|".$aa_neighbor_con;
    ## OMIM match
    my $omim_match_con = OMIM_Compare($omim, $gene, $mafrstart, $mafres1, $mafres2, 2);
    my $cosmic_omim = $cosmic_aa_results."\t".$omim_match_con;
    # return both cosmic and omim results
    return($cosmic_omim);
}
## OMIM database comparing 
sub OMIM_Compare {
    my ($omim, $gene, $start, $aa0, $aa1, $aa_range) = @_;
    my $perfect_con     = "NA";
    my $samlocation_con = "NA";
    my $neighbor_con    = "NA";
    my $its = $aa0.$start.$aa1;
    # perfect match
    my $ad = exists($omim->{$gene}->{$start}->{$its});
    if ($ad) { $perfect_con = $its; }
    # position match
    my $ad0 = exists($omim->{$gene}->{$start});
    if ($ad0) {
        my @t = keys %{$omim->{$gene}->{$start}};
        $samlocation_con = "";
        foreach my $it0 (@t) {
            next if ($it0 eq $its);
            $samlocation_con .= $it0.",";
        }
        chop($samlocation_con) if ($samlocation_con);
        $samlocation_con = "NA" if ($samlocation_con eq "");
    }
    # proximity match
    my $iter_start = $start - $aa_range;
    my $iter_stop  = $start + $aa_range;
    $iter_start = 1 if ($iter_start <= 0);
    $neighbor_con = "";
    foreach my $it0 ($iter_start..$iter_stop) {
        next if ($it0 == $start);
        my $ad1 = exists($omim->{$gene}->{$it0});
        if ($ad1) {
            my @t = keys %{$omim->{$gene}->{$it0}};
            foreach my $it0 (@t) { $neighbor_con .= $it0.","; }
        }
    }
    chop($neighbor_con) if ($neighbor_con);
    $neighbor_con ="NA" if ($neighbor_con eq "");
    my $omim_results = $perfect_con."|".$samlocation_con."|".$neighbor_con;
    return($omim_results);
}
# Amino Acid change check
sub AA_Check {
    my ($aachange) = @_;
    my ($res1, $start, $res2, $stop, $new_res);
    unless (defined($aachange)) { return ($res1, $start, $res2, $stop, $new_res); }
    unless ($aachange =~ /^p\./) { return ($res1, $start, $res2, $stop, $new_res); }
    #__FORMULATE ERROR STRING JUST IN CASE
    #__VALIDATE
     $aachange =~ s/^p\.//x;
     if ($aachange =~ /^ (\D+) (\d+) _ (\D+) (\d+) (.*) $/x ) {
         ($res1, $start, $res2, $stop, $new_res) = ($1, $2, $3, $4, $5);
     }
     elsif ($aachange =~ /^ (\D+) (\d+) (\D+) (.*) $/x ) {
         ($res1, $start, $res2, $new_res) = ($1, $2, $3, $4);
         $stop = $start;
     }
     elsif ($aachange =~ /^ (\d+) (.*) $/x ) {
         ($start, $res2) = ($1, $2);
         $res1 = '*';
         $stop = $start;
         $new_res = $res2;
     } 
     elsif ($aachange =~ /^ (\D+) (\d+) $/x ) {
         ($res1, $start) = ($1, $2);
         $new_res = ' ';
         $res2 = ' ';
         $stop = $start;
     }
     if (defined($new_res)) { $new_res =~ s/^ > //x; }
     $new_res ||= '';
     return ($res1, $start, $res2, $stop, $new_res);
}
# Parse the head of input file
sub FileHeadParse {
    my ($fh, $wuanno, $vcf) = @_;
    # Parse input file and get head info
    my @cols;
    my $header;
    my $temp_header;
    my %wu_nah  = ();
    my $total   = 1;
    $temp_header = "\t"."COSMIC_Nuc\tCOSMIC_AA\tOMIM\tGDSC\tDGENE\n";
    # appended VCF header information
    my $c1_info = "##INFO=<ID=COS_NUC,Number=.,Type=String,Description=\"Cosmic nucleotide matches. Format: Exact_matches|Position_matches|Proximity_matches\">\n";
    my $c2_info = "##INFO=<ID=COS_AA,Number=.,Type=String,Description=\"Cosmic Amino Acid Change matches. Format: Exact_matches|Position_matches|Proximity_matches\">\n";
    my $om_info = "##INFO=<ID=OMIM,Number=.,Type=String,Description=\"Cosmic nucleotide matches. Format: Exact_matches|Position_matches|Proximity_matches\">\n";
    my $gd_info = "##INFO=<ID=GDSC,Number=.,Type=String,Description=\"Target drugs matches. Format: Target_drugs|Target_drugs_based_on_pvalue\">\n";
    my $dgene_info = "##INFO=<ID=DGENE,Number=.,Type=String,Description=\"Drugable gene matches. Format: Gene|Describtion\">\n";
    if ($vcf) {
        while (my $ll = <$fh>) {
            if ($ll =~ m/^#CHROM/) {
                # added vcf header
                $header .= $c1_info.$c2_info.$om_info.$gd_info.$dgene_info;
                chomp($ll);
                my $i = 0;
                %wu_nah = map {($_, $i++)} split(/\t/, $ll);
                $header .= $ll.$temp_header;
                last;
            }
            else {
                $header .= $ll;
            }
        }
        unless (    defined($wu_nah{"#CHROM"}) 
                and defined($wu_nah{"POS"}) 
                and defined($wu_nah{"ID"})                           
                and defined($wu_nah{"REF"}) 
                and defined($wu_nah{"ALT"})
                and defined($wu_nah{"QUAL"}) 
                and defined($wu_nah{"FILTER"})
                and defined($wu_nah{"INFO"})) {
            die "not a valid VCF annotation file!\n";
        }
        @cols = ($wu_nah{"#CHROM"}, 
                 $wu_nah{"POS"},
                 $wu_nah{"REF"}, 
                 $wu_nah{"ALT"}, 
                 $wu_nah{"INFO"});
        return($total, $header, \@cols, $wu_nah{"INFO"});
    }
    if ($wuanno) {
        while (my $ll = <$fh>) {
            if ($ll =~ m/^Hugo_Symbol/) {
                chomp($ll);
                my $i = 0;
                %wu_nah = map {($_, $i++)} split(/\t/, $ll);
                $header .= $ll.$temp_header;
                last;
           }
           else {
               $header .= $ll;
           }
       }
       unless(    defined($wu_nah{"gene_name"}) 
              and defined($wu_nah{"chromosome_name"}) 
              and defined($wu_nah{"start"})
              and defined($wu_nah{"stop"}) 
              and defined($wu_nah{"reference"})
              and defined($wu_nah{"variant"}) 
              and defined($wu_nah{"transcript_name"})
              and defined($wu_nah{"strand"}) 
              and defined($wu_nah{"amino_acid_change"})) {
           die "not a valid WU annotation file!\n";
       } 
       @cols = ($wu_nah{"gene_name"}, 
                $wu_nah{"chromosome_name"}, 
                $wu_nah{"start"},
                $wu_nah{"stop"}, 
                $wu_nah{"reference"}, 
                $wu_nah{"variant"},
                $wu_nah{"transcript_name"}, 
                $wu_nah{"strand"},
                $wu_nah{"amino_acid_change"});
       return($total, $header, \@cols);
    }
    # MAF file
    while (my $ll = <$fh>) {
       if ($ll =~ m/^Hugo_Symbol/) {
            chomp($ll);
            my $i = 0;
            %wu_nah = map {($_, $i++)} split(/\t/, $ll);
            $header .= $ll.$temp_header;
            last;
       }
       else {
           $header .= $ll;
       } 
    }
    unless (    defined($wu_nah{"Hugo_Symbol"}) 
            and defined($wu_nah{"Chromosome"}) 
            and defined($wu_nah{"Start_Position"})                                
            and defined($wu_nah{"End_Position"}) 
            and defined($wu_nah{"Reference_Allele"})                                
            and defined($wu_nah{"Tumor_Seq_Allele1"}) 
            and defined($wu_nah{"Tumor_Seq_Allele2"})) {
        die "not a valid MAF annotation file!\n";
    }
    @cols = ($wu_nah{"Hugo_Symbol"}, 
             $wu_nah{"Chromosome"}, 
             $wu_nah{"Start_Position"},                        
             $wu_nah{"End_Position"}, 
             $wu_nah{"Reference_Allele"}, 
             $wu_nah{"Tumor_Seq_Allele1"},                        
             $wu_nah{"Tumor_Seq_Allele2"}, 
             $wu_nah{"transcript_name"}, 
             $wu_nah{"strand"},
             $wu_nah{"amino_acid_change"});
    unless(     defined($wu_nah{"transcript_name"}) 
            and defined($wu_nah{"strand"}) 
            and defined($wu_nah{"amino_acid_change"})) {
        $total = 0;
    }
    return($total, $header, \@cols);
}
# Parse info of lines 
sub ParseCols {
    my ($line, $total, $rcols, $wuanno, $vcf) = @_;
    # For vcf 
    my ($pos, $vcfref, $vcfalt, $vcfinfo, @vcfnuc, @vcfres); 
    # For MAF
    my ($a1, $a2);
    my ($chr, $start, $stop, $ref, $var);
    my ($gene, $trans, $strand, $aa) = ("NA", "NA", "NA", "NA");
    my @cols = split(/\t/, $line);
    # VCF annotation
    if ($vcf) {
        ($chr, $pos, $vcfref, $vcfalt, $vcfinfo) = @cols[@$rcols];
        unless ($vcfinfo =~ /CSQ=/) {
            die "input is not valid Ensembl VEP annotation file\n";
        }
        VariantVCF($chr, $pos, $vcfref, $vcfalt, $vcfinfo, \@vcfnuc, \@vcfres);
        return(\@vcfnuc, \@vcfres);
    }
    # WU annotation
    if ($wuanno) {
        ($gene, $chr, $start, $stop, 
         $ref, $var, $trans, $strand, $aa) = @cols[@$rcols];
        return($gene, $chr, $start, $stop, $ref, $var, $trans, $strand, $aa);
    }
    # MAF annotation
    if ($total) {
        ($gene, $chr, $start, $stop, 
         $ref, $a1, $a2, 
         $trans, $strand, $aa) = @cols[@$rcols];
    }
    else {
        ($gene, $chr, $start, 
         $stop, $ref, $a1, $a2) = @cols[@$rcols[0..6]];
    }
    $var = ($ref eq $a1 ? $a2 : $a1);
    return($gene, $chr, $start, $stop, $ref, $var, $trans, $strand, $aa);
}


# Parse variants of VCF 
sub VariantVCF {
    my ($vcfchr, $pos, $vcfref, $vcfalt, $vcfinfo, $vcfnuc, $vcfres) = @_;
    my ($start, $stop, $ref, $var);
    my @t = split(/,/, $vcfalt);
    foreach my $item (@t) {
        if (length($vcfref) == length($item)) {
            $start = $pos;
            $stop  = $pos + length($item) - 1;
            $ref   = $vcfref;
            $var   = $item;
        } 
        elsif (length($vcfref) > length($item)) {
            $start = $pos + 1;
            $stop  = $pos + length($vcfref) - 1;
            $ref   = substr($vcfref, 1);
            $var   = "-";
        }
        else {
            $start = $pos;
            $stop  = $pos + 1;
            $ref   = "-";
            $var   = substr($item, 1);;
        }
        push(@$vcfnuc, $vcfchr."|".$start."|".$stop."|".$ref."|".$var);
    }
    # amino acid part
    $vcfinfo =~ /CSQ=(\S+)/;
    my @t0 = split(/,/, $1);
    foreach my $it0 (@t0) {
        $it0 .= " ";
        my @t1 = split(/\|/, $it0);
        my ($allel, $trans, $transt, $pos, $res, $gene) = @t1[0,2,3,7,8,13];
        unless ($transt eq "Transcript") {next;}
        unless ($trans) { $trans = "NA"; } 
        unless ($gene) { $gene = "NA"; }
        my $aa = "NA";
        if ($pos and $res) {
            if ($res =~ /(\S+)\/(\S+)/) {
                $aa = "p.$1$pos$2";
            }
            else {
                $aa = "p.$res$pos";
            }
        }
        push(@$vcfres, $trans."|".$aa."|".$gene."|".$vcfchr);
    }
    return 1;
}

1;

