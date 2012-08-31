package Genome::Model::Tools::Annotate::ParseVep;

use strict;
use warnings;
use Genome;
use IO::File;
use File::Basename;
use FileHandle;

class Genome::Model::Tools::Annotate::ParseVep {
    is => 'Command',
    has => [
	    output_file => {
	        is => 'Text',
	        is_optional => 0,
	        doc => "Outputs annotation-like file",
	    },
        vep_input => {
            is => 'Text',
            is_optional => 0,
            doc => "VEP annotation file",
        },
        strict_sites => {
            is => 'Text',
            is_optional => 1,
            doc => "A tab-delimited list of sites to limit the output to, \"chr\\tpos\"",
        },
	    reference_build	=> {
            is => 'Text',
            doc => "reference build -- \"NCBI-human-build36\" or \"GRCh37-lite-build37\"",
            is_optional => 1,
            default => 'GRCh37-lite-build37',
            is_input => 1
        },
	],
};


sub help_brief {                            # keep this to just a few words <---
    "This tool is designed to complement vcf-to-variant-matrix in how it parses pipeline files to a format better suited for the statistical analysis methods. Change the VEP default output into something more readable. Also, try to pick a \"top\" annotation."
}


sub help_synopsis {
<<'HELP';
This tool is designed to complement vcf-to-variant-matrix in how it parses pipeline files to a format better suited for the statistical analysis methods.
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
<<'HELP';
This tool is designed to complement vcf-to-variant-matrix in how it parses pipeline files to a format better suited for the statistical analysis methods.
HELP
}

###############

sub execute {                               # replace with real execution logic.
	my $self = shift;

    my $vep_annotation = $self->vep_input;
    my $output_annotation = $self->output_file;
    my $reference_build = $self->reference_build;

    my $inFh_vep_annot = Genome::Sys->open_file_for_reading($vep_annotation);
    my $output_annotation_fh = Genome::Sys->open_file_for_writing($output_annotation);

    my $only_strict = 0;
    my $strict_sites;
    my %strict_sites_hash;
    if ($self->strict_sites) {
        $only_strict = 1;
        $strict_sites = $self->strict_sites;
        my $inFh_strict_sites = Genome::Sys->open_file_for_reading($strict_sites);
        while(my $line = $inFh_strict_sites->getline) {
            chomp $line;
            my ($chr, $pos) = split "\t", $line;
            my $strict_matcher = "$chr.$pos";
            $strict_sites_hash{$strict_matcher}++;
        }
    }

    my %variant_hash;
    my %gene_name_hash;
    while(my $line = $inFh_vep_annot->getline) {
        chomp $line;
        if ($line =~ m/^##/) { #these are info lines
            next;
        }
        if ($line =~ m/^#/) { #this is the header line
            next;
        }
        #Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
        my ($Uploaded_variation,$Location,$Allele,$Gene,$Feature,$Feature_type,$Consequence,$cDNA_position,$CDS_position,$Protein_position,$Amino_acids,$Codons,$Existing_variation,$extra) = split "\t", $line;

        if ($only_strict) {
            my ($chr, $pos) = split ":", $Location;
            my $strict_matcher = "$chr.$pos";
            unless (defined $strict_sites_hash{$strict_matcher}) {
                next;
            }
        }

		my ($chrom, $chr_start, $alleles, $ref, $var);
        if ($Uploaded_variation =~ m/rs/) {
            ($chrom, $chr_start) = split(/:/,$Location);
            $var = $Allele;
            my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "$reference_build");
            my $reference_build_fasta = $reference_build_fasta_object->data_directory . "/all_sequences.fa";
            $ref = `samtools faidx $reference_build_fasta $chrom:$chr_start-$chr_start | grep -v ">"`;
            chomp($ref);
        }
        else {
			($chrom, $chr_start, $alleles) = split(/\_/, $Uploaded_variation);
			($ref, $var) = split(/\//, $alleles);
        }
        my $variant_name = $chrom."_".$chr_start."_".$ref."_".$var;

		## Reset extra variables
        my $gene_name = my $condel_prediction = my $polyphen_prediction = my $sift_prediction = "-";
        my @gene_info = split ";", $extra;
        if ($extra ne "-") {
            foreach my $entry (@gene_info) {
    			my ($key, $value) = split(/\=/, $entry);
                if ($entry =~ m/HGNC/) {
                    ($gene_name) = $entry =~ m/=(\S+)/;
                }
                else {
                    $entry =~ s/ /_/;
                    if ($entry =~ m/SIFT/) {
                        ($sift_prediction) = $entry =~ m/=(\D+)\(/;
                    }
                    elsif ($entry =~ m/PolyPhen/) {
                        ($polyphen_prediction) = $entry =~ m/=(\D+)\(/;
                    }
                    elsif ($entry =~ m/Condel/) {
                        ($condel_prediction) = $entry =~ m/=(\D+)\(/;
                    }
                    else {
                        print "Thing = $entry\n$line\n";
                    }
                }
            }
        }

        my ($aa1, $aa2) = split "/", $Amino_acids;
        my $aa_change;
        if ($Amino_acids eq '-') {
            $aa_change = "-";
        }
        else {
            if (defined $aa1 && defined $aa2) {
                if (defined $gene_name) {
                    $aa_change = $gene_name.".".$aa1.$Protein_position.$aa2;
                }
                else {
                    $aa_change = "-".".".$aa1.$Protein_position.$aa2;
                }
            }
            elsif (defined $aa1 && !defined $aa2) {
                if (defined $gene_name) {
                    $aa_change = $gene_name.".".$aa1.$Protein_position;
                }
                else {
                    $aa_change = "-".".".$aa1.$Protein_position;
                }
            }
            else {
                $aa_change = $gene_name.".".$aa1.$Protein_position.$aa2;
            }
        }

        if (defined($variant_hash{$variant_name})) {
            my $prediction1 = "$variant_hash{$variant_name}";
            my $prediction2 = "$gene_name\t$Consequence\t$aa_change\t$condel_prediction\t$polyphen_prediction\t$sift_prediction";
            my $choice = $self->rank_VEP($prediction1,$prediction2,\%gene_name_hash);
            if ($choice == 1) {
                $variant_hash{$variant_name} = $prediction1;
                my ($old_gene_name,$old_Consequence,$old_aa_change,$old_condel_prediction,$old_polyphen_prediction,$old_sift_prediction) = split "\t", $variant_hash{$variant_name};
                $gene_name_hash{$old_gene_name}++;
            }
            elsif ($choice == 2) {
                $variant_hash{$variant_name} = $prediction2;
                $gene_name_hash{$gene_name}++;
            }

        }
        else {
            $variant_hash{$variant_name} = "$gene_name\t$Consequence\t$aa_change\t$condel_prediction\t$polyphen_prediction\t$sift_prediction";
        }
    }
    print $output_annotation_fh "Variant_Name\tChromosome\tPosition\tReference\tVariant\tGene_Name\tTrv_Type\taa_change\tcondel_prediction\tpolyphen_prediction\tsift_prediction\n";
    foreach my $variant_name (sort keys %variant_hash) {
        my ($chrom,$position,$ref,$var) = split "_", $variant_name;
        print $output_annotation_fh "$variant_name\t$chrom\t$position\t$ref\t$var\t$variant_hash{$variant_name}\n";
    }

    return 1;
}

sub rank_VEP {
	my $self = shift;
    my $prediction1 = shift;
    my $prediction2 = shift;
    my $hash_ref = shift;
    my %gene_name_hash = %{$hash_ref};
    
    my ($old_gene_name,$old_Consequence,$old_aa_change,$old_condel_prediction,$old_polyphen_prediction,$old_sift_prediction) = split "\t", $prediction1;
    my ($new_gene_name,$new_Consequence,$new_aa_change,$new_condel_prediction,$new_polyphen_prediction,$new_sift_prediction) = split "\t", $prediction2;
    my $old_prediction = "$old_condel_prediction\t$old_polyphen_prediction\t$old_sift_prediction";
    my $new_prediction = "$new_condel_prediction\t$new_polyphen_prediction\t$new_sift_prediction";
    my ($prediction_winner) = $self->byPrediction($old_prediction,$new_prediction);
    if ($prediction_winner) {
        return $prediction_winner;
    }
    else {
	    my ($code_a,$code_b) = $self->byTrv($old_Consequence,$new_Consequence);
        if ($code_a > $code_b) {
            return 1;
        }
        elsif ($code_b > $code_a) {
            return 2;
        }
        else {
            my ($aa_ranking) = $self->byAminoAcid($old_aa_change,$new_aa_change);
            if ($aa_ranking) {
                return $aa_ranking;
            }
            else {
                my ($gene_ranking) = $self->byGene($old_gene_name,$new_gene_name,\%gene_name_hash);
                if ($gene_ranking) {
                    return $gene_ranking;
                }
                else {
                    print "No ranking possible :( (all ties!) so I just randomly picked first one listed here:\n\t$prediction1\n\t$prediction2\n";
                    return 1;
                }
            }
        }
    }
}

sub byPrediction {
	my $self = shift;
    my $old_prediction = shift;
    my $new_prediction = shift;

    my ($old_condel_prediction,$old_polyphen_prediction,$old_sift_prediction) = split(/\t/, $old_prediction);
    my ($new_condel_prediction,$new_polyphen_prediction,$new_sift_prediction) = split(/\t/, $new_prediction);
    if ($old_condel_prediction ne "-" && $new_condel_prediction ne "-") {
        my $prediction_sum1 = my $prediction_sum2 = "benign";
        if ($old_condel_prediction =~ m/damaging/ || $old_condel_prediction =~ m/deleterious/) {
            $prediction_sum1 = "bad";
        }
        if ($new_condel_prediction =~ m/damaging/ || $new_condel_prediction =~ m/deleterious/) {
            $prediction_sum2 = "bad";
        }

        if ($prediction_sum1 eq 'bad' && $prediction_sum2 eq 'bad') {
            return 0;
        }
        elsif ($prediction_sum1 eq 'bad') {
            return 1;
        }
        elsif ($prediction_sum2 eq 'bad') {
            return 2;
        }
        else {
            return 0;
        }
    }
    elsif ($old_condel_prediction ne "-") {
        return 1;
    }
    elsif ($new_condel_prediction ne "-") {
        return 2;
    }
    elsif ($old_polyphen_prediction ne "-" && $new_polyphen_prediction ne "-") {
        my $prediction_sum1 = my $prediction_sum2 = "benign";
        if ($old_polyphen_prediction =~ m/damaging/ || $old_polyphen_prediction =~ m/deleterious/) {
            $prediction_sum1 = "bad";
        }
        if ($new_polyphen_prediction =~ m/damaging/ || $new_polyphen_prediction =~ m/deleterious/) {
            $prediction_sum2 = "bad";
        }

        if ($prediction_sum1 eq 'bad' && $prediction_sum2 eq 'bad') {
            return 0;
        }
        elsif ($prediction_sum1 eq 'bad') {
            return 1;
        }
        elsif ($prediction_sum2 eq 'bad') {
            return 2;
        }
        else {
            return 0;
        }
    }
    elsif ($old_polyphen_prediction ne "-") {
        return 1;
    }
    elsif ($new_polyphen_prediction ne "-") {
        return 2;
    }
    elsif ($old_sift_prediction ne "-" && $new_sift_prediction ne "-") {
        my $prediction_sum1 = my $prediction_sum2 = "benign";
        if ($old_sift_prediction =~ m/damaging/ || $old_sift_prediction =~ m/deleterious/) {
            $prediction_sum1 = "bad";
        }
        if ($new_sift_prediction =~ m/damaging/ || $new_sift_prediction =~ m/deleterious/) {
            $prediction_sum2 = "bad";
        }

        if ($prediction_sum1 eq 'bad' && $prediction_sum2 eq 'bad') {
            return 0;
        }
        elsif ($prediction_sum1 eq 'bad') {
            return 1;
        }
        elsif ($prediction_sum2 eq 'bad') {
            return 2;
        }
        else {
            return 0;
        }
    }
    elsif ($old_sift_prediction ne "-") {
        return 1;
    }
    elsif ($new_sift_prediction ne "-") {
        return 2;
    }
    return 0;
}

sub byTrv {
	my $self = shift;
    my $choice_one = shift;
    my $choice_two = shift;

	my $code_a = my $code_b = 0;

    my %vep_class_rank = ();
    $vep_class_rank{'-'} = 0;
    $vep_class_rank{'NMD_TRANSCRIPT'} = 0;
    $vep_class_rank{'PARTIAL_CODON'} = 0;
    $vep_class_rank{'INTERGENIC'} = 0;
    $vep_class_rank{'UPSTREAM'} = 1;
    $vep_class_rank{'DOWNSTREAM'} = 2;
    $vep_class_rank{'INTRONIC'} = 3;
    $vep_class_rank{'5PRIME_UTR'} = 4;
    $vep_class_rank{'3PRIME_UTR'} = 5;
    $vep_class_rank{'WITHIN_NON_CODING_GENE'} = 6;
    $vep_class_rank{'WITHIN_MATURE_miRNA'} = 7;
    $vep_class_rank{'SYNONYMOUS_CODING'} = 8;
    $vep_class_rank{'SPLICE_SITE'} = 9;
    $vep_class_rank{'ESSENTIAL_SPLICE_SITE'} = 10;
    $vep_class_rank{'NON_SYNONYMOUS_CODING'} = 11;
    $vep_class_rank{'STOP_LOST'} = 12;
    $vep_class_rank{'STOP_GAINED'} = 13;


	foreach my $test_class (keys %vep_class_rank)
	{
        if ($choice_one =~ m/$test_class/) {
            if ($vep_class_rank{$test_class} > $code_a) {
    	        $code_a = $vep_class_rank{$test_class};
            }
        }
        if ($choice_two =~ m/$test_class/) {
            if ($vep_class_rank{$test_class} > $code_b) {
        		$code_b = $vep_class_rank{$test_class};
            }
        }
	}
	return ($code_a,$code_b);
}

sub byAminoAcid {
	my $self = shift;
    my $old_aa_change = shift;
    my $new_aa_change = shift;
    if ($old_aa_change ne "-" && $new_aa_change ne "-") {
        return 0;
    }
    elsif ($old_aa_change ne "-" && $old_aa_change =~ m/\S+/) {
        return 1;
    }
    elsif ($new_aa_change ne "-" && $new_aa_change =~ m/\S+/) {
        return 2;
    }
    return 0;
}
sub byGene {
	my $self = shift;
    my $old_gene_name = shift;
    my $new_gene_name = shift;
    my $hash_ref = shift;
    my %gene_name_hash = %{$hash_ref};
    if ($old_gene_name ne "-" && $new_gene_name ne "-") {
        if (defined($gene_name_hash{$old_gene_name}) && defined($gene_name_hash{$new_gene_name})) {
            return 0;
        }
        elsif (defined($gene_name_hash{$old_gene_name})) {
            return 1;
        }
        elsif (defined($gene_name_hash{$new_gene_name})) {
            return 2;
        }
        return 0;
    }
    elsif ($old_gene_name ne "-" && $old_gene_name =~ m/\S+/) {
        return 1;
    }
    elsif ($new_gene_name ne "-" && $new_gene_name =~ m/\S+/) {
        return 2;
    }
    return 0;
}

sub bin_consequence {
	my $self = shift;
	my $Consequence = shift;

    ###BINNING TO EXAMINE TRV CONSEQUENCE###
    $Consequence =~ s/STOP_GAINED,SPLICE_SITE/Truncation/;
    $Consequence =~ s/STOP_GAINED/Truncation/;
    $Consequence =~ s/ESSENTIAL_SPLICE_SITE/Truncation/;

    $Consequence =~ s/NON_SYNONYMOUS_CODING,SPLICE_SITE/Non-truncation/;
    $Consequence =~ s/NON_SYNONYMOUS_CODING/Non-truncation/;
    $Consequence =~ s/STOP_LOST/Non-truncation/;

    $Consequence =~ s/UPSTREAM/Regulatory_Change/;
    $Consequence =~ s/3PRIME_UTR,SPLICE_SITE/Regulatory_Change/;
    $Consequence =~ s/5PRIME_UTR,SPLICE_SITE/Regulatory_Change/;
    $Consequence =~ s/3PRIME_UTR/Regulatory_Change/;
    $Consequence =~ s/5PRIME_UTR/Regulatory_Change/;

    $Consequence =~ s/SPLICE_SITE,INTRONIC/Unknown_Protein_Effect/;
    $Consequence =~ s/SPLICE_SITE,WITHIN_NON_CODING_GENE/Unknown_Protein_Effect/;
    $Consequence =~ s/WITHIN_NON_CODING_GENE,INTRONIC/Unknown_Protein_Effect/;
    $Consequence =~ s/WITHIN_NON_CODING_GENE/Unknown_Protein_Effect/;
    $Consequence =~ s/INTRONIC/Unknown_Protein_Effect/;
    $Consequence =~ s/^SYNONYMOUS_CODING,SPLICE_SITE/Unknown_Protein_Effect/;
    $Consequence =~ s/^SYNONYMOUS_CODING/Unknown_Protein_Effect/;
    $Consequence =~ s/DOWNSTREAM/Unknown_Protein_Effect/;
    $Consequence =~ s/INTERGENIC/Unknown_Protein_Effect/;

    return $Consequence;
}

1;
