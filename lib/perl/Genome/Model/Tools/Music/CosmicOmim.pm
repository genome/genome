package Genome::Model::Tools::Music::CosmicOmim;

#__STANDARD PERL PACKAGES
   use warnings;
   use strict;
   use Genome;
   use FileHandle;
   use Text::CSV_XS;

class Genome::Model::Tools::Music::CosmicOmim {
    is => 'Genome::Model::Tools::Music::Base',
    has_input => [
       maf_file => {
           is => 'Path', file_format => 'maf',
           doc => 'list of annotated mutations in MAF format (or any file with MAF+annotation headers)',
       },
       output_file => {
           is => 'Path', is_output => 1,
           doc => 'Output file contains the input file with two columns appended to the end, corresponding to cosmic and omim mutation comparisons, respectively',
       },
       reference_build => {
           is => 'Text',
           example_values => ['Build36', 'Build37'],
           doc => 'Put either "Build36" or "Build37"',
       }
    ],
    has_optional_input => [
       omimaa_dir => {
           is => 'Path', default => Genome::Sys->dbpath( 'omim', 'latest' ),
           doc => 'omim amino acid mutation database folder',
       },
       cosmic_dir => {
           is => 'Path', default => Genome::Sys->dbpath( 'cosmic', 'latest' ),
           doc => 'cosmic amino acid mutation database folder',
       },
       verbose => {
           is => 'Boolean', default => 0,
           doc => 'Use this to display the larger working output',
       },
       wu_annotation_headers => {
           is => 'Boolean', default => 0,
           doc => "Use this if input MAF contains WUSTL annotation format headers",
       },
       aa_range => {
           is => 'Integer', default => 2,
           doc => "Set how close a 'near' match is when searching for amino acid near hits",
       },
       nuc_range => {
           is => 'Integer', default => 5,
           doc => "Set how close a 'near' match is when searching for nucleotide position near hits",
       },
       show_known_hits => {
           is => 'Boolean', default => 1,
           doc => "When a finding is novel, show known AA in that gene",
       },
    ],
    doc => 'Compare the amino acid changes of supplied mutations to COSMIC and OMIM databases.'
};

sub help_synopsis {
    return <<EOS;
 ... music cosmic-omim \\
        --maf-file input_dir/myMAF.tsv \\
        --output-file output_dir/myMAF_output.tsv \\
        --no-verbose

 ... music cosmic-omim \\
        --maf-file input_dir/myMAF.tsv \\
        --output-file output_dir/myMAF_output.tsv \\
        --omimaa-dir omim_dir/ \\
        --cosmic-dir cosmic_dir/ \\
        --no-verbose
EOS
}

sub help_detail {
    return <<EOS;
This tool looks at the amino acid changes for the given set of mutations and compares the genomic coordinates as well as the affected amino acid to the coordinates and amino acids of all cancer-specific mutations listed in the Cosmic and OMIM databases. The database files are specially prepared for this task and provided with the MuSiC suite. The tool reports various types of matches, including matches within "near proximity", where "near proximity" is currently defined as a linear DNA distance of 5 bases or 2 amino acids. (This type of matching helps to account for the possibility of subtle differences in reported positions for variants due to differences in transcript definitions or other things of this nature.) Any site without a match in a particular databases is reported as "novel" with respect to that database.

The output of this script returns each row the original input MAF file with two columns appended to the end of each, one column for each of the databases. Also included is a STDOUT printout of a summary of what was found in the input MAF. Neither output can be suppressed in the current version. The --verbose option is used to display working notes that are useful for various purposes in debugging potential MAF problems. The Omim and Cosmic directories must point to the output of the downloader, named appropriately, as they don't recognize the OMIM database in the raw download format.

This tool only compares build 36 or build 37 coordinates that are specified in Cosmic to the coordinates in your maf file. This is a weakness of the Cosmic database (not all position entries are currently available for both builds) that is out of our control.

In addition to the standard version 2.3 MAF headers, there needs to be 3 columns appended. These
column headers in the MAF must have these names in the header in order for the tool to find them:
   transcript_name - the transcript name, such as NM_000028
 amino_acid_change - the amino acid change, such as p.R290H
        c_position - the nucleotide position changed, such as c.869

EOS
}

sub _doc_credits {
    return <<EOS
This tool depends on copies of data from the following databases, packaged in a form useable for quick analysis:

 * COSMIC - http://www.sanger.ac.uk/genetics/CGP/cosmic/
 * OMIM - http://www.ncbi.nlm.nih.gov/omim
EOS
}

sub _doc_authors {
    return <<EOS
 Brian Dunford-Shore
 David Larson, Ph.D.
 Michael C. Wendl, Ph.D.
 William Schierding, M.S.
EOS
}

sub doc_copyright_years {
    return (2007,2011);
}

sub execute {
    my $self = shift;

####################
#  PRE-PROCESSING  #
####################

    my $mut_file = $self->maf_file;
    my $cosmic_dir = $self->cosmic_dir;

    if (-e "$cosmic_dir/music") {
        # new versions of COSMIC data contain multiple data sets besides the MuSiC files.
        $cosmic_dir = "$cosmic_dir/music";
    }

    my $basename = $self->output_file;
    my $ref_build = $self->reference_build;
    unless($ref_build =~ m/build36/i xor $ref_build =~ m/build37/i){
        $self->error_message("You must either specify reference_build as either \"Build36\" or \"Build37\"");
        die $self->error_message;
    }
    my $omimaa_dir = $self->omimaa_dir;
    my $omimaa;
    if (-d $omimaa_dir){
        $omimaa = "$omimaa_dir/omim_aa.csv";
        unless (-e $omimaa) {
            $omimaa = "$omimaa_dir/OMIM_aa_will.csv";
        }
    }
    my $cosmic_database_file = 'Cosmic_Database.tsv';
    my $cosmic_database;
    if (-d $cosmic_dir){
        $cosmic_database = "$cosmic_dir/$cosmic_database_file";
    }
    my $verbose = $self->verbose;
    my $wuheaders = $self->wu_annotation_headers;

    $self->debug_message("Using $omimaa as omima db file");
    $self->debug_message("Using $cosmic_dir as cosmic db folder");

    my $aa_range = $self->aa_range;
    my $nuc_range = $self->nuc_range;

#####################
#  MAIN PROCESSING  #
#####################

#Set Stats hash that counts everything
    my %stats;
    $stats{'OMIMDB'} = 0;
    $stats{'OMIM'}{'doublematch'} = 0;
    $stats{'OMIM'}{'ntmatch'} = 0;
    $stats{'OMIM'}{'aamatch'} = 0;
    $stats{'OMIM'}{'posmatch'} = 0;
    $stats{'OMIM'}{'ntposmatch'} = 0;
    $stats{'OMIM'}{'aaposmatch'} = 0;
    $stats{'OMIM'}{'nearmatch'} = 0;
    $stats{'OMIM'}{'ntnearmatch'} = 0;
    $stats{'OMIM'}{'aanearmatch'} = 0;
    $stats{'OMIM'}{'novel'} = 0;
    $stats{'OMIM'}{'ntnovel'} = 0;
    $stats{'OMIM'}{'aanovel'} = 0;
    $stats{'OMIM'}{'silent'} = 0;
    $stats{'OMIM'}{'nomatch'} = 0;

    $stats{'COSMICDB'} = 0;
    $stats{'COSMIC'}{'doublematch'} = 0;
    $stats{'COSMIC'}{'ntmatch'} = 0;
    $stats{'COSMIC'}{'aamatch'} = 0;
    $stats{'COSMIC'}{'posmatch'} = 0;
    $stats{'COSMIC'}{'ntposmatch'} = 0;
    $stats{'COSMIC'}{'aaposmatch'} = 0;
    $stats{'COSMIC'}{'nearmatch'} = 0;
    $stats{'COSMIC'}{'ntnearmatch'} = 0;
    $stats{'COSMIC'}{'aanearmatch'} = 0;
    $stats{'COSMIC'}{'novel'} = 0;
    $stats{'COSMIC'}{'ntnovel'} = 0;
    $stats{'COSMIC'}{'aanovel'} = 0;
    $stats{'COSMIC'}{'silent'} = 0;
    $stats{'COSMIC'}{'nomatch'} = 0;

#__PARSE MUTATION DATA
    my $fh = new FileHandle;
    unless ($fh->open (qq{$mut_file})) { die "Could not open mutation project file '$mut_file' for reading"; }

    if ($verbose) {print "Parsing mutation file...\n";}
    my $mutation = ParseMutationFile($fh, $mut_file, $wuheaders);
    $fh->close;
    if ($verbose) {print "Done Parsing Mutation File! Yippee!\n";}

#__READ IN OMIM FILE
    my %omimaa;
    if (defined($omimaa) && -e $omimaa) {
        open(OMIMAA,$omimaa) || die "Could not open omim file '$omimaa'";
        if ($verbose) {print "Loading OMIM Database\n";}
        my $omimaa_header = <OMIMAA>;
        while(<OMIMAA>) {
            $stats{'OMIMDB'}++;
            chomp;
            my ($gene, $omim_entry, $position, $aa_ori, $aa_mut, $description, $diseases) = split("\t");
            $omimaa{$gene}{$omim_entry}{$position}{residue1} = $aa_ori;
            $omimaa{$gene}{$omim_entry}{$position}{residue2} = $aa_mut;
            $omimaa{$gene}{$omim_entry}{$position}{description} = $description;
        }
        if ($verbose) {print "Finished Loading OMIM Database\n";}
        close(OMIMAA);
    }

#__READ IN COSMIC FILE
    my %cosmic_gene;
    my %cosmic_position;
    my %cosmic_position_only;
    my %cosmic_tissue;
    my %cosmic_histology;
    my %aa_count;
    my %residue_match;
    if (defined($cosmic_database) && -e $cosmic_database) {
        my ($gene_col, $chr_col, $start_col, $stop_col, $chr_col_37, $start_col_37, $stop_col_37, $amino_col, $nucleo_col, $somatic_col, $primary_tissue_col, $tissue_sub_1_col, $tissue_sub_2_col, $histology_col,
        $histology_sub_1_col, $histology_sub_2_col, $gene, $chr, $start, $stop, $chr_37, $start_37, $stop_37, $amino, $nucleo, $somatic, $primary_tissue, $tissue_sub_1, $tissue_sub_2, $histology, $histology_sub_1,
        $histology_sub_2);
        open(COSMIC,$cosmic_database) || die "Could not open omim file '$cosmic_database'";
        if ($verbose) {print "Loading COSMIC Database\n";}
        my $cosmic_header = <COSMIC>;
        chomp($cosmic_header);
        my @parser = split(/\t/, $cosmic_header);
        my $parsecount = 0;
        my %parsehash;

        foreach my $item (@parser) {
           $parsehash{$item} = $parsecount;
           $parsecount++;
        }

        $gene_col = $parsehash{'Gene'};
        $chr_col = $parsehash{'Chromosome'};
        $start_col = $parsehash{'Genome Start'};
        $stop_col = $parsehash{'Genome Stop'};
        $chr_col_37 = $parsehash{'Chromosome Build37'};
        $start_col_37 = $parsehash{'Genome Start Build37'};
        $stop_col_37 = $parsehash{'Genome Stop Build37'};
        $amino_col = $parsehash{'Amino Acid'};
        $nucleo_col = $parsehash{'Nucleotide'};
        $somatic_col = $parsehash{'Somatic Status'};
        $primary_tissue_col = $parsehash{'Primary_Tissue'};
        $tissue_sub_1_col = $parsehash{'Tissue_subtype_1'};
        $tissue_sub_2_col = $parsehash{'Tissue_subtype_2'};
        $histology_col = $parsehash{'Histology'};
        $histology_sub_1_col = $parsehash{'Histology_subtype_1'};
        $histology_sub_2_col = $parsehash{'Histology_subtype_2'};

       while(my $line = <COSMIC>) {
            $stats{'COSMICDB'}++;
            chomp($line);
            my @parser = split(/\t/, $line);
            $gene = $parser[$gene_col];
            if ($ref_build =~ m/build36/i) {
                $chr = $parser[$chr_col];
                $start = $parser[$start_col];
                $stop = $parser[$stop_col];
            }
            else {
                $chr = $parser[$chr_col_37];
                $start = $parser[$start_col_37];
                $stop = $parser[$stop_col_37];
            }
            $chr =~ s/NotListed/ /;
            $start =~ s/NotListed/ /;
            $stop =~ s/NotListed/ /;

            $chr =~ s/23/X/;
            $chr =~ s/24/Y/;

            $amino = $parser[$amino_col];
            $nucleo = $parser[$nucleo_col];
            $somatic = $parser[$somatic_col];
            $primary_tissue = $parser[$primary_tissue_col];
            $tissue_sub_1 = $parser[$tissue_sub_1_col];
            $tissue_sub_2 = $parser[$tissue_sub_2_col];
            $histology = $parser[$histology_col];
            $histology_sub_1 = $parser[$histology_sub_1_col];
            $histology_sub_2 = $parser[$histology_sub_2_col];
            my ($residue1, $res_start, $residue2, $res_stop, $new_residue);
            if ($amino eq 'E597A') {
                $amino = 'p.E597A';
            }

            if ($amino =~ m/^p/) {
                unless ($amino =~ m/^p\./) {
                    $amino =~ s/^p/p\./;
                }
            }

            if ($amino =~ m/^P\./) {
                $amino =~ s/^P\./p\./;
            }

            $amino =~ s/\s+$//;
            $amino =~ s/^\s+//;

            ($residue1, $res_start, $residue2, $res_stop, $new_residue) = AA_Check($amino);

            if (defined $res_start){
                if ($res_start == $res_stop){
                    $residue_match{$gene}{$res_start}{$res_stop}{$residue1}{$residue2}++;
                    my $addition = $residue1.$res_start.$residue2;
                    $aa_count{$gene}{$addition}++;
                }
                else {
                    my $addition = $residue1.$res_start."-".$res_stop.$residue2;
                    $aa_count{$gene}{$addition}++;
                }
            }
            else {
                if ($amino =~ m/^p\.\?/ || $amino =~ m/^p\.\>/ || $amino eq 'p.INS' || $amino eq 'p.DEL' || $amino eq 'p.fs*?' || $amino eq 'p.fs') {
                    my $addition = $amino;
                    $addition =~ s/p\.//;
                    $aa_count{$gene}{$addition}++;
                }
                else {
                    print "$residue1, $res_start, $residue2, $res_stop, $new_residue\n";
                    die "$amino not found\n";
                }
            }


            if (defined($nucleo)) {
                $cosmic_gene{$gene}{$amino}{$nucleo}++;
            }
            else {
                $cosmic_gene{$gene}{$amino}++;
            }

            if (defined($chr) && defined($start) && defined($stop) && $chr ne ' ' && $start ne ' ' && $stop ne ' ') {
                $cosmic_position{$chr}{$start}{$stop}{$gene}++;
                if (defined($nucleo)) {
                    $cosmic_position_only{$chr}{$start}{$stop}{$nucleo}++;
                }
                else {
                    $cosmic_position_only{$chr}{$start}{$stop}++;
                }
            }

            $cosmic_tissue{$gene}{$amino} = "$primary_tissue\t$tissue_sub_1\t$tissue_sub_2";
            $cosmic_histology{$gene}{$amino} = "$histology\t$histology_sub_1\t$histology_sub_2";
        }
        close(COSMIC);
        if ($verbose) {print "Finished Loading COSMIC Database! Hooray!\n";}
    }

    my %cosmic_results;
    my %omim_results;

    my $summary_file = $basename;
    unless (open(SUMMARY,">$summary_file")) {
    die "Could not open output file '$summary_file' for writing";
    }

    unless ($fh->open (qq{$mut_file})) { die "Could not open mutation project file '$mut_file' for reading"; }
    my %fileline;
    my $i = 1;
    while (my $filehandleline = <$fh>) {
        chomp $filehandleline;
        while ($filehandleline =~ /^#/) {
            print SUMMARY "$filehandleline\n";
            $filehandleline = <$fh>;
            chomp $filehandleline;
        }
        $fileline{$i} = $filehandleline;
        $i++;
    }
    $fh->close;

    print SUMMARY "$fileline{'1'}\tInput_MAF_Line_Number\tCosmic_Results\tOMIM_Results\n";
    if ($verbose) {print "Starting COSMIC/OMIM to Mutation File Comparisons\n";}
    foreach my $hugo (sort keys %{$mutation}) {
    foreach my $sample (keys %{$mutation->{$hugo}}) {
    foreach my $line_num (keys %{$mutation->{$hugo}->{$sample}}) {
        if ($verbose) {print ".";}   #report that we are starting a sample

        #read in the alleles. The keys may change with future file formats. If so, a new version should be added to
#       my ($entrez_gene_id, $line, $aa_change,$transcript,$mstatus,$Variant_Type,$Chromosome,$Start_position,$End_position,$Reference_Allele,$Tumor_Seq_Allele1,$gene) =
        my ($entrez_gene_id, $line, $aa_change,$transcript,$Variant_Type,$Chromosome,$Start_position,$End_position,$Reference_Allele,$Tumor_Seq_Allele1,$Tumor_Seq_Allele2,$gene) =
        (
        $mutation->{$hugo}->{$sample}->{$line_num}->{ENTREZ_GENE_ID},
        $mutation->{$hugo}->{$sample}->{$line_num}->{file_line},
        $mutation->{$hugo}->{$sample}->{$line_num}->{AA_CHANGE},
        $mutation->{$hugo}->{$sample}->{$line_num}->{TRANSCRIPT},
#        $mutation->{$hugo}->{$sample}->{$line_num}->{MUTATION_STATUS},
        $mutation->{$hugo}->{$sample}->{$line_num}->{VARIANT_TYPE},
        $mutation->{$hugo}->{$sample}->{$line_num}->{CHROMOSOME},
        $mutation->{$hugo}->{$sample}->{$line_num}->{START_POSITION},
        $mutation->{$hugo}->{$sample}->{$line_num}->{END_POSITION},
        $mutation->{$hugo}->{$sample}->{$line_num}->{REFERENCE_ALLELE},
        $mutation->{$hugo}->{$sample}->{$line_num}->{TUMOR_SEQ_ALLELE1},
        $mutation->{$hugo}->{$sample}->{$line_num}->{TUMOR_SEQ_ALLELE2},
        $mutation->{$hugo}->{$sample}->{$line_num}->{HUGO_SYMBOL},
        );

#        if ($mstatus){
        #Annotate the allele's effect on all known (ie transcript without the 'unknown' status) transcripts
        my $proper_allele;
        if ($Reference_Allele ne $Tumor_Seq_Allele1) {
            $proper_allele = $Tumor_Seq_Allele1;
        }
        elsif($Reference_Allele ne $Tumor_Seq_Allele2) {
            $proper_allele = $Tumor_Seq_Allele2;
        }
        else {
            warn "line $line_num has both Tumor_Seq_Allele1 ($Tumor_Seq_Allele1) and Tumor_Seq_Allele2 ($Tumor_Seq_Allele2) as the Reference ($Reference_Allele)\n";
            $proper_allele = $Tumor_Seq_Allele2;
        }

#LOOK FOR ONLY SINGLE CHARACTER PROPER ALLELE TYPES - A, C, T, G, 0, or -
        unless($Reference_Allele ne $proper_allele) {
            die "Ref allele: $Reference_Allele same as mutation allele: $proper_allele ('line num' $line_num)";
        }
        unless($Reference_Allele =~ /[ACTG0\-]/ && $proper_allele =~ /[ACTG0\-]/) {
            die "Read in improper alleles from mutation file ref: $Reference_Allele var: $proper_allele ('line num' $line_num)";
        }

## SHOULD NEVER BE '--' BECAUSE THE CODE WILL THINK THIS IS A DINUCLEOTIDE POLYMORPHISM INSTEAD OF AN INDEL, SO CHANGE THESE
        $Reference_Allele = '-' if $Reference_Allele eq '--';
        $proper_allele = '-' if $proper_allele eq '--';
        chomp($line_num);
        my %results_hash;
        #parse the amino acid string
        my ($residue1, $res_start, $residue2, $res_stop, $new_residue) = AA_Check($aa_change);
        if(defined($Start_position) && defined($End_position) && (!$residue2 || $residue2 eq ' ')){
            if ($verbose) {print "Skipping Non-coding (Silent) Mutation";}
            my $find_type = $self->CheckPositionMatch($Start_position, $End_position, $Reference_Allele, $proper_allele, \%cosmic_position_only, \%cosmic_position);
            my $createspreadsheet;
            if ($find_type =~ m/no_match/) {
                $createspreadsheet = "$fileline{$line_num}\t$line_num\tSkipped - Silent Mutation\tSkipped - Silent Mutation";
                $stats{'COSMIC'}{'silent'}++;
                $stats{'OMIM'}{'silent'}++;
            }
            elsif ($find_type =~ m/position/) {
                my ($findtype,$cosmic_genes,$chr,$gen_start,$gen_stop) = split(":",$find_type);
                $results_hash{NT}{POSITION}{COSMIC}{$transcript}=": Nucleotide -> Cosmic Gene(s):$cosmic_genes Position:$chr:$gen_start-$gen_stop";
                my $matchtype_cosmic;
                ($cosmic_results{$line_num}, $matchtype_cosmic) = score_results(\%results_hash, "COSMIC");
                $stats{'COSMIC'}{$matchtype_cosmic}++;
                $stats{'OMIM'}{'silent'}++;
                $createspreadsheet = "$fileline{$line_num}\t$line_num\t$cosmic_results{$line_num}\tSkipped - Silent Mutation";
            }
            elsif ($find_type =~ m/position_nucleotide/) {
                my ($findtype,$cosmic_genes,$chr,$gen_start,$gen_stop,$reference,$mutant) = split(":",$find_type);
                $results_hash{NT}{MATCH}{COSMIC}{$transcript} = ": Nucleotide -> Cosmic Gene(s): $cosmic_genes Position:$chr:$gen_start-$gen_stop,ref:$reference,mut:$mutant";
                my $matchtype_cosmic;
                ($cosmic_results{$line_num}, $matchtype_cosmic) = score_results(\%results_hash, "COSMIC");
                $stats{'COSMIC'}{$matchtype_cosmic}++;
                $stats{'OMIM'}{'silent'}++;
                $createspreadsheet = "$fileline{$line_num}\t$line_num\t$cosmic_results{$line_num}\tSkipped - Silent Mutation";
            }
            elsif ($find_type =~ m/near_match/) {
                my ($findtype,$diff_position) = split(":",$find_type);
                $results_hash{NT}{ALMOST}{COSMIC}{$transcript}=": Nucleotide $diff_position bases away";
                my $matchtype_cosmic;
                ($cosmic_results{$line_num}, $matchtype_cosmic) = score_results(\%results_hash, "COSMIC");
                $stats{'COSMIC'}{$matchtype_cosmic}++;
                $stats{'OMIM'}{'silent'}++;
                $createspreadsheet = "$fileline{$line_num}\t$line_num\t$cosmic_results{$line_num}\tSkipped - Silent Mutation";
            }
            print SUMMARY "$createspreadsheet\n";
            next; #skip silent mutations
        }

        #look for mutation file gene ($hugo) to match cosmic file gene (%cosmic_gene{gene}), note matched name in $cosmic_hugo
        my $cosmic_hugo;
        my $uc_hugo = uc($hugo);
        if (exists($cosmic_gene{$hugo})) {
            $cosmic_hugo = $hugo;
        }
        elsif (exists($cosmic_gene{$uc_hugo})) {       # check for UPPERCASE hugo match
            $cosmic_hugo = $uc_hugo;
        }
        else {
            #if cosmic key needs to be uppercase to match mutation file (for example, maf default is all uppercase)
            foreach my $key (keys %cosmic_gene) {       # check for UPPERCASE keys match
                if ($uc_hugo eq uc($key)) {
                    $cosmic_hugo = $key;
                }
            }
        }

        #genes that didn't find a match go here, will check for position matches later
        unless (defined($cosmic_hugo)) {
            if (-e "$cosmic_dir/$hugo\.csv") { #database flatfile has only genes with AA changes, check source files for gene existance
                $results_hash{NT}{NOVEL}{COSMIC}{$transcript}=": Gene $hugo in Cosmic but No Amino Acid Results for Gene";
            }
            else {
                $results_hash{NT}{NOVEL}{COSMIC}{$transcript}=": Gene $hugo not in Cosmic Database";
            }
        }

        #retrieve COSMIC match
        if (!defined($aa_change) || $aa_change eq 'NULL') {
            warn "We skipped silent mutations, how do we have undel or NULL amino acids at non-silent sites? Gene:$hugo AA:$aa_change\n";
        }

        #Start checks,First check position (then check amino acid)

        my $cosmic_find_type;
        my @aa_holder;

        if(defined($Start_position) && defined($End_position) && $Start_position ne ' ' && $End_position ne ' ') {
            my $find_type = $self->CheckPositionMatch($Start_position, $End_position, $Reference_Allele, $proper_allele, \%cosmic_position_only, \%cosmic_position);
            if ($find_type =~ m/no_match/) {
                $results_hash{NT}{NOVEL}{COSMIC}{$transcript}=": Nucleotide";
            }
            elsif ($find_type =~ m/position/) {
                my ($findtype,$cosmic_genes,$chr,$gen_start,$gen_stop) = split(":",$find_type);
                $results_hash{NT}{POSITION}{COSMIC}{$transcript}=": Nucleotide -> Cosmic Gene(s):$cosmic_genes Position:$chr:$gen_start-$gen_stop";
            }
            elsif ($find_type =~ m/position_nucleotide/) {
                my ($findtype,$cosmic_genes,$chr,$gen_start,$gen_stop,$reference,$mutant) = split(":",$find_type);
                $results_hash{NT}{MATCH}{COSMIC}{$transcript} = ": Nucleotide -> Cosmic Gene(s): $cosmic_genes Position:$chr:$gen_start-$gen_stop,ref:$reference,mut:$mutant";
            }
            elsif ($find_type =~ m/near_match/) {
                my ($findtype,$diff_position) = split(":",$find_type);
                $results_hash{NT}{ALMOST}{COSMIC}{$transcript}=": Nucleotide $diff_position bases away";
            }
        }

        #check amino acid here
        if ($cosmic_hugo && exists($cosmic_gene{$cosmic_hugo})) {
            $cosmic_find_type = 'no_match';
            foreach my $key (keys %{$aa_count{$cosmic_hugo}}) {
                if ($key =~ m/\S+/) {
                    push(@aa_holder,"$key ($aa_count{$cosmic_hugo}{$key})");
                }
            }
            unless (@aa_holder) {
                @aa_holder = "AA?Unknown?";
            }
            if ($res_start && $res_stop && exists($residue_match{$cosmic_hugo}) && exists($residue_match{$cosmic_hugo}{$res_start}) && exists($residue_match{$cosmic_hugo}{$res_start}{$res_stop})) { #match amino acid
                $cosmic_find_type = 'position';
                $results_hash{AA}{POSITION}{COSMIC}{$transcript}=": Amino Acid -> Matched $cosmic_hugo, $res_start, $res_stop";
                if ($residue1 && $residue2 && ($residue_match{$cosmic_hugo}{$res_start}{$res_stop}{$residue1}{$residue2} || $residue_match{$cosmic_hugo}{$res_start}{$res_stop}{uc($residue1)}{$residue2} || $residue_match{$cosmic_hugo}{$res_start}{$res_stop}{$residue1}{uc($residue2)} || $residue_match{$cosmic_hugo}{$res_start}{$res_stop}{uc($residue1)}{uc($residue2)})) { # matches both amino acid and amino acid position
                    $cosmic_find_type = 'position_aminoacid';
                    my $addition;
                    if ($res_start == $res_stop){
                        $addition = $residue1.$res_start.$residue2;
                    }
                    else {
                        $addition = $residue1.$res_start."-".$res_stop.$residue2;
                    }
                    $results_hash{AA}{MATCH}{COSMIC}{$transcript} = ": Amino Acid -> Matched $cosmic_hugo, $addition";
                }
            }
            if($cosmic_find_type && $cosmic_find_type eq 'no_match') {
                if ($self->show_known_hits) {
                    $results_hash{AA}{NOVEL}{COSMIC}{$transcript}=": Amino Acid -> Known AA = @aa_holder";
                }
                else {
                    $results_hash{AA}{NOVEL}{COSMIC}{$transcript}=": Amino Acid -> Known AA Not Shown";
                }
                my $iter_start = $res_start - $aa_range;
                my $iter_stop = $res_stop + $aa_range;
                my $iter;
                for($iter = $iter_start; $iter <= $iter_stop; $iter++) {
                    if ($res_start && $res_stop && exists($residue_match{$cosmic_hugo}) && exists($residue_match{$cosmic_hugo}{$iter})) {
                        if ($self->show_known_hits) {
                            $results_hash{AA}{ALMOST}{COSMIC}{$transcript} = ": Amino Acid -> Known AA for Gene = @aa_holder";
                        }
                        else {
                            $results_hash{AA}{ALMOST}{COSMIC}{$transcript} = ": Amino Acid -> Known AA Not Shown";
                        }
                    }
                }
            }
        }

        #retrieve OMIM match
        my $omim_find_type;
        my $omim = \%omimaa;
        if (exists($omim->{$hugo})) {
            $omim_find_type = FindOMIM(\%omimaa,$hugo,$res_start,$res_stop,$residue1,$residue2,$aa_range);
        }
        #Add OMIM result to the results hash
        if(defined($omim_find_type)) {
            if ($omim_find_type eq 'position_aminoacid') {
                $results_hash{AA}{MATCH}{OMIM}{$transcript} =": Amino Acid";
            }
            elsif ($omim_find_type eq 'position') {
                $results_hash{AA}{POSITION}{OMIM}->{$transcript}=": Amino Acid";
            }
            elsif ($omim_find_type eq 'almost') {
                $results_hash{AA}{ALMOST}{OMIM}->{$transcript}=": Amino Acid";
            }
            else {
                $results_hash{AA}{NOVEL}{OMIM}->{$transcript}=": Amino Acid";
            }
        }
        else {
            $results_hash{AA}{NOVEL}{OMIM}->{$transcript}=": Amino Acid - OMIM Gene Name Not Found";
        }

        #now check to see what the 'best' cosmic score was
        my $matchtype_cosmic;
        ($cosmic_results{$line_num}, $matchtype_cosmic) = score_results(\%results_hash, "COSMIC");
        $stats{'COSMIC'}{$matchtype_cosmic}++;
        #now check to see what the 'best' omim score was
        my $matchtype_omim;
        ($omim_results{$line_num}, $matchtype_omim) = score_results(\%results_hash, "OMIM");
        $stats{'OMIM'}{$matchtype_omim}++;

        my $createspreadsheet = "$fileline{$line_num}\t$line_num\t$cosmic_results{$line_num}\t$omim_results{$line_num}";
        print SUMMARY "$createspreadsheet\n";
#        }
    }
    }
    }
    close(SUMMARY);
    if ($verbose) {print "Finished COSMIC/OMIM to Mutation File Comparisons! HAPPY HAPPY JOY JOY!\n";}

    print "\n";
    print "Number of Genes in OMIM: $stats{'OMIMDB'}\n";
    print "Number of AA and NT Matches: $stats{'OMIM'}{'doublematch'}\n";
    print "Number of NT only Matches: $stats{'OMIM'}{'ntmatch'}\n";
    print "Number of AA only Matches: $stats{'OMIM'}{'aamatch'}\n";
    print "Number of AA and NT Position Matches: $stats{'OMIM'}{'posmatch'}\n";
    print "Number of NT only Position Matches: $stats{'OMIM'}{'ntposmatch'}\n";
    print "Number of AA only Position Matches: $stats{'OMIM'}{'aaposmatch'}\n";
    print "Number of NT and AA Novel Sites with Matches in Near Proximity: $stats{'OMIM'}{'nearmatch'}\n";
    print "Number of NT only Novel Sites with Matches in Near Proximity: $stats{'OMIM'}{'ntnearmatch'}\n";
    print "Number of AA only Novel Sites with Matches in Near Proximity: $stats{'OMIM'}{'aanearmatch'}\n";
    print "Number of NT and AA Novel Sites with Nothing in Near Proximity: $stats{'OMIM'}{'novel'}\n";
    print "Number of NT Novel Sites with Nothing in Near Proximity: $stats{'OMIM'}{'ntnovel'}\n";
    print "Number of AA Novel Sites with Nothing in Near Proximity: $stats{'OMIM'}{'aanovel'}\n";
    print "Number of Silent Mutations Skipped: $stats{'OMIM'}{'silent'}\n";
    print "Number of Lines that Exited with No Hit: $stats{'OMIM'}{'nomatch'}\n";
    print "\n";
    print "Number of Genes in COSMIC: $stats{'COSMICDB'}\n";
    print "Number of AA and NT Matches: $stats{'COSMIC'}{'doublematch'}\n";
    print "Number of NT only Matches: $stats{'COSMIC'}{'ntmatch'}\n";
    print "Number of AA only Matches: $stats{'COSMIC'}{'aamatch'}\n";
    print "Number of AA and NT Position Matches: $stats{'COSMIC'}{'posmatch'}\n";
    print "Number of NT only Position Matches: $stats{'COSMIC'}{'ntposmatch'}\n";
    print "Number of AA only Position Matches: $stats{'COSMIC'}{'aaposmatch'}\n";
    print "Number of NT and AA Novel Sites with Matches in Near Proximity: $stats{'COSMIC'}{'nearmatch'}\n";
    print "Number of NT only Novel Sites with Matches in Near Proximity: $stats{'COSMIC'}{'ntnearmatch'}\n";
    print "Number of AA only Novel Sites with Matches in Near Proximity: $stats{'COSMIC'}{'aanearmatch'}\n";
    print "Number of NT and AA Novel Sites with Nothing in Near Proximity: $stats{'COSMIC'}{'novel'}\n";
    print "Number of NT Novel Sites with Nothing in Near Proximity: $stats{'COSMIC'}{'ntnovel'}\n";
    print "Number of AA Novel Sites with Nothing in Near Proximity: $stats{'COSMIC'}{'aanovel'}\n";
    print "Number of Silent Mutations Skipped: $stats{'COSMIC'}{'silent'}\n";
    print "Number of Lines that Exited with No Hit: $stats{'COSMIC'}{'nomatch'}\n";

    return 1;
}




################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################


sub CheckPositionMatch {
    my $self = shift;
    my $genomic_start = shift;
    my $genomic_stop = shift;
    my $nt1 = shift;
    my $nt2 = shift;
    my $cp_only = shift;
    my %cosmic_position_only = %{$cp_only};
    my $cp_all = shift;
    my %cosmic_position = %{$cp_all};
    my $verbose = $self->verbose;
    my $nuc_range = $self->nuc_range;

    my $find_type = 'no_match';
    foreach my $chr (sort keys %cosmic_position_only) {
        # Test that it at least matches position
        foreach my $gen_start (keys %{$cosmic_position_only{$chr}}) {
            my $diff_start = $gen_start - $genomic_start;
            if ($gen_start == $genomic_start) {
                foreach my $gen_stop (keys %{$cosmic_position_only{$chr}{$gen_start}}) {
                    my $diff_stop = $gen_stop - $genomic_stop;
                    if ($gen_stop == $genomic_stop) {
                        my $cosmic_genes;
                        if (keys %{$cosmic_position{$chr}{$gen_start}{$gen_stop}}) {
                            my @cosmic_genes = keys %{$cosmic_position{$chr}{$gen_start}{$gen_stop}};
                            $cosmic_genes = join(",",@cosmic_genes);
                        }
                        if ($find_type eq 'no_match' || $find_type eq 'near_match') {
                            $find_type = "position:$cosmic_genes:$chr:$gen_start:$gen_stop";
                        }
                        # Test that it matches both
                        foreach my $nucleo (keys %{$cosmic_position_only{$chr}{$gen_start}{$gen_stop}}) {
                            my ($start,$stop,$type_length,$type,$reference,$mutant) = parse_nucleotide($nucleo,$verbose);
                            if($reference && $mutant && $reference eq $nt1 && $mutant eq $nt2) {
                                $find_type = "position_nucleotide:$cosmic_genes:$chr:$gen_start:$gen_stop:$reference:$mutant";
                                return $find_type;
                            }
                        }
                    }
                    elsif ($diff_stop <= $nuc_range && $diff_stop >= -$nuc_range) {
                        if ($find_type eq 'no_match') {
                            $find_type = "near_match:$diff_stop";
                        }
                    }
                }
            }
            elsif ($diff_start <= $nuc_range && $diff_start >= -$nuc_range) {
                if ($find_type eq 'no_match') {
                    $find_type = "near_match:$diff_start";
                }
            }
        }
    }
    return $find_type;
}

sub FindOMIM {
    my ($omim, $hugo,$res_start, $res_stop, $residue1, $residue2, $aa_range) = @_;

    my $return_value = 'no_match';
    unless (exists($omim->{$hugo})) {
        warn "No omim entry for: $hugo";
        return $return_value;
    }
    foreach my $sample (keys %{$omim->{$hugo}}) {
        # Test that it at least matches position
        if (exists($omim->{$hugo}{$sample}{$res_start})) {
            $return_value = 'position';
            # Test that it matches both
            if (exists($omim->{$hugo}{$sample}{$res_start}{residue1}) &&
                    exists($omim->{$hugo}{$sample}{$res_start}{residue2}) &&
                    defined($omim->{$hugo}{$sample}{$res_start}{residue1}) &&
                    defined($omim->{$hugo}{$sample}{$res_start}{residue2}) &&
                    uc($omim->{$hugo}{$sample}{$res_start}{residue1}) eq uc($residue1) &&
                    uc($omim->{$hugo}{$sample}{$res_start}{residue2}) eq uc($residue2)) {
                return 'position_aminoacid';
            }
        }
        elsif ($return_value eq 'no_match') {
            my $iter_start = $res_start - $aa_range;
            my $iter_stop = $res_stop + $aa_range;
            my $iter;
            for($iter = $iter_start; $iter <= $iter_stop; $iter++) {
                if (exists($omim->{$hugo}{$sample}{$iter})) {
                    $return_value = 'almost';
                }
            }
        }
    }
    return $return_value;
}

sub score_results {
    my ($results, $database) = @_;
    my $matchtype;
    if(exists($results->{NT}{MATCH}->{$database}) && exists($results->{AA}{MATCH}->{$database})) {
#best hit was a DOUBLE MATCH. Huzzah!
        $matchtype = 'doublematch';
        my ($transcript) = keys %{$results->{NT}{MATCH}{$database}};
        my ($transcript2) = keys %{$results->{AA}{MATCH}{$database}};
        my $ret_value = "NT and AA Match".$results->{AA}{MATCH}->{$database}{$transcript2}." and ".$results->{NT}{MATCH}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{NT}{MATCH}->{$database})) {
#best hit was a MATCH. Huzzah!
        $matchtype = 'ntmatch';
        my ($transcript) = keys %{$results->{NT}{MATCH}{$database}};
        my $ret_value = "Match".$results->{NT}{MATCH}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{AA}{MATCH}->{$database})) {
#best hit was a MATCH. Huzzah!
        $matchtype = 'aamatch';
        my ($transcript) = keys %{$results->{AA}{MATCH}{$database}};
        my $ret_value = "Match".$results->{AA}{MATCH}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{NT}{POSITION}->{$database}) && exists($results->{AA}{POSITION}->{$database})) {
#best hit was a position match
        $matchtype = 'posmatch';
        my ($transcript) = keys %{$results->{NT}{POSITION}{$database}};
        my ($transcript2) = keys %{$results->{AA}{POSITION}{$database}};
        my $ret_value = "NT and AA Position Match".$results->{AA}{POSITION}->{$database}{$transcript2}." and ".$results->{NT}{POSITION}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{NT}{POSITION}->{$database})) {
#best hit was a position match
        $matchtype = 'ntposmatch';
        my ($transcript) = keys %{$results->{NT}{POSITION}{$database}};
        my $ret_value = "Position Match".$results->{NT}{POSITION}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{AA}{POSITION}->{$database})) {
#best hit was a position match
        $matchtype = 'aaposmatch';
        my ($transcript) = keys %{$results->{AA}{POSITION}{$database}};
        my $ret_value = "Position Match".$results->{AA}{POSITION}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{NT}{ALMOST}->{$database}) && exists($results->{AA}{ALMOST}->{$database})) {
#best hit was near a position match
        $matchtype = 'nearmatch';
        my ($transcript) = keys %{$results->{NT}{ALMOST}{$database}};
        my ($transcript2) = keys %{$results->{AA}{ALMOST}{$database}};
        my $ret_value = "NT and AA Novel, but near match".$results->{AA}{ALMOST}->{$database}{$transcript2}." and ".$results->{NT}{ALMOST}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{NT}{ALMOST}->{$database})) {
#best hit was near a position match
        $matchtype = 'ntnearmatch';
        my ($transcript) = keys %{$results->{NT}{ALMOST}{$database}};
        my $ret_value = "Novel, but near match".$results->{NT}{ALMOST}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{AA}{ALMOST}->{$database})) {
#best hit was near a position match
        $matchtype = 'aanearmatch';
        my ($transcript) = keys %{$results->{AA}{ALMOST}{$database}};
        my $ret_value = "Novel, but near match".$results->{AA}{ALMOST}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{AA}{NOVEL}->{$database}) && exists($results->{NT}{NOVEL}->{$database})) {
#no hits, novel
        $matchtype = 'novel';
        my ($transcript) = keys %{$results->{AA}{NOVEL}{$database}};
        my ($transcript2) = keys %{$results->{NT}{NOVEL}{$database}};
        my $ret_value = "Novel".$results->{AA}{NOVEL}->{$database}{$transcript}. " and ".$results->{NT}{NOVEL}->{$database}{$transcript2};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{AA}{NOVEL}->{$database})) {
#no hits, novel
        $matchtype = 'aanovel';
        my ($transcript) = keys %{$results->{AA}{NOVEL}{$database}};
        my $ret_value = "Novel".$results->{AA}{NOVEL}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    elsif(exists($results->{NT}{NOVEL}->{$database})) {
#no hits, novel
        $matchtype = 'ntnovel';
        my ($transcript) = keys %{$results->{NT}{NOVEL}{$database}};
        my $ret_value = "Novel".$results->{NT}{NOVEL}->{$database}{$transcript};
        return ($ret_value, $matchtype);
    }
    else {
#it was a nomatch! this shouldn't happen.
        $matchtype = 'nomatch';
        my ($transcript) = keys %{$results->{NOMATCH}};
        my $ret_value = (defined($transcript) && $results->{NOMATCH}{$transcript}) ?  $results->{NOMATCH}{$transcript} : "Unknown/NULL";
        return ($ret_value, $matchtype);
    }
}

sub parse_nucleotide {
    my ($string, $verbose) = @_;
    my ($change, $modifier);
    my ($start,$stop,$type_length,$type,$reference,$mutant);
    ($start,$stop,$change) = $string =~ /^c\. (\d+) _ (\d+) (.*) $/x;
    unless (defined $change) {
        ($start,$change) = $string =~ /^c\. (\d+) (.*) $/x;
    }
    if ($string =~ /^c\.\w* \Q?\E.* /x) {
#ambiguous entry
        return;
    }
    if( defined $change) {
#first check to make sure its not intronic
        if ($change =~ /^ (-|\+) (\d+) (.*) $/x) {
            if ($verbose) {print "Mutation Appears Intronic";}
            $change = $3;
        }
#This could be one of several possiblities
#First distinguish between > notation and del18 notation
        ($reference,$mutant) = $change =~ /^ (\D*) > (\D+) $/x;
        if(defined($reference) && $reference eq '') {
            $reference = undef;
        }
        if(defined $mutant) {
#Then we expect that this format is correct.
            if((defined($stop) &&  $stop eq $start) || !defined($stop))  {
#We have a snp
                $type = 'SNP';
                $type_length = 1;
            }
            else {
#assuming that if it is listed explicitly it is an indel
                $type = 'indel';
                if(defined $reference) {
                    $type_length = length $reference;
                }
                else {
                    $type_length = abs($start-$stop)+1;
                }
            }
            return ($start, $stop, $type_length, $type, $reference, $mutant);
        }
        else {
#did not guess right. Should be either del15 or insAAT type of format
            $type = substr $change, 0,3;
            $modifier = substr $change, 3, (length($change) - 1);
            if($type =~ /^ (del|ins|delins) $/xi) {
                if($type eq 'delins') {
                    $type = 'indel';
                }
#then insertion
                if($modifier =~ /^ (\d+) $/x) {
#it is a digit
                    $type_length = $1;
                    return ($start, $stop, $type_length, $type, $reference, $mutant);
                }
                else {

#it is a sequence
                    $type_length = length $modifier;
                    return ($start, $stop, $type_length, $type, $reference, $modifier);
                }
            }
            else {
#unrecognized format
                warn "Unable to parse nucleotide format in: $string\n";
                return;
            }
        }
    }
    else {
        warn "Unable to parse nucleotide format in: $string\n";
        return;
    }
}

sub ParseMutationFile {
    my ($fh,$file,$wuheaders) = @_;
    my $source = 'mutation_csv'; # 'CSV'
    my $keyfields = 'HUGO_SYMBOL:TUMOR_SAMPLE_ID:file_line_num';
    my $field_subset_array = undef;
    my $header_translation;
    if ($wuheaders) {
        $header_translation = {
            'chromosome_name' => 'CHROMOSOME',
            'start' => 'START_POSITION',
            'stop' => 'END_POSITION',
            'reference' => 'REFERENCE_ALLELE',
            'variant' => 'TUMOR_SEQ_ALLELE1',
#                'variant' => 'TUMOR_SEQ_ALLELE2',
            'type' => 'VARIANT_TYPE',
            'gene_name' => 'HUGO_SYMBOL',
            'transcript_name' => 'TRANSCRIPT',
            'strand' => 'TUMOR_SAMPLE_ID', #meaningless proxy
            'amino_acid_change' => 'AA_CHANGE',
        };
    }
    else {
        $header_translation = {
			'Chromosome' => 'CHROMOSOME',
			'Start_position' => 'START_POSITION',
			'End_position' => 'END_POSITION',
			'Reference_Allele' => 'REFERENCE_ALLELE',
			'Tumor_Seq_Allele1' => 'TUMOR_SEQ_ALLELE1',
			'Tumor_Seq_Allele2' => 'TUMOR_SEQ_ALLELE2',
			'Variant_Type' => 'VARIANT_TYPE',
			'Hugo_Symbol' => 'HUGO_SYMBOL',
			'transcript_name' => 'TRANSCRIPT',
			'Strand' => 'TUMOR_SAMPLE_ID', #meaningless proxy
			'amino_acid_change' => 'AA_CHANGE',
        };
    }
    my $header_skip = 0;
    my $no_header = 0;
    my $header_fields = undef;
    my $ucheader_fields = undef;
    my $separator = "\t";
    my $no_spaces = 0;
    my $line_number_field = 'file_line_num';
    my $line_field = 'file_line';
    my $line_num = 1;
    my $record = {};
    my $csv = Text::CSV_XS->new({'sep_char' => $separator});
    #parse MAF header
    my $header = <$fh>;
    while ($header =~ /^#/) {
        $header = <$fh>;
    }
    unless ($header =~ m/Hugo_Symbol/) {
        die "Header field \"Hugo_Symbol\" not found. Therefore, your file is assumed to not have a header. If this is in error, please contact the authors to fix this line in the code. Otherwise, please include a header in your input file. Thank you.\n";
    }
    $line_num++;
    $csv->parse($header);
    my @header_fields = $csv->fields();
    if ($ucheader_fields) {
        @header_fields = map { uc($_) } @header_fields;
    }
    if(defined($no_spaces) && $no_spaces) {
        @header_fields = map { $_ =~ s/ /_/g } @header_fields;
    }
    unless(defined($separator)) {
        $header =~ s/\t/,/gx;
    }
    unshift( @header_fields, $line_field);        # Add 'extra' fields of the input line
    unshift( @header_fields, $line_number_field);    # and the line number
    # Translate the header names, if a translation is given
    for (my $h = 0; $h <= $#header_fields; $h++) {
        $header_fields[$h] = (exists($header_translation->{$header_fields[$h]})) ?
            $header_translation->{$header_fields[$h]} : $header_fields[$h];
    }

    my (@key_fields) = split(':',$keyfields);
    my %key_fields;
    @key_fields{ @key_fields } = ( 0 .. $#key_fields );

    # Construct field name to position lookup
    my %header_fields;
    @header_fields{ @header_fields } = ( 0 .. $#header_fields );
    # Construct a subset of the fields--the default is the complete set of fields
    my %field_subset;
    if (defined($field_subset_array)) {
        @field_subset{ @{$field_subset_array} } = @header_fields{ @{$field_subset_array } };
    } else {
        @field_subset{ @header_fields } = ( 0 .. $#header_fields );
    }
    # Construct a list of fields that are not key (are values only)
    my @value_fields;
    foreach my $field (@header_fields[ (values %field_subset ) ]) {
        unless (exists($key_fields{$field})) {
            push @value_fields, ($field);
        }
    }
    #__PARSE FILE
    my $line;
    while ($line = <$fh>) {
        chomp ($line);
        my $temp = $line;

        #edit for punctuation characters in annotation section of MAF
        $temp =~ s/"//g;
        $temp =~ s/'//g;
        unless(defined($separator)) {
            #maintain original default behavior of handling both tabs and
            #commas
            $temp =~ s/\t/,/gx;
        }
            $csv->parse($temp);
        my @values = $csv->fields();

            if(defined($no_spaces) && $no_spaces) {
                @values = map { $_ =~ s/ /_/g } @values;
            }

        unshift (@values, $line);        # Add 'extra' fields of the input line
        unshift (@values, $line_num++);    # and the line number

        my $sub_record;
        # Construct the hierarchical key structure
        $sub_record = $record;
        foreach my $sub_key  (@key_fields) {
            unless (exists($sub_record->{ $values[ $header_fields{ $sub_key } ] } )) {
                $sub_record->{ $values[ $header_fields{ $sub_key } ] } = {};
            }
            $sub_record = $sub_record->{ $values[ $header_fields{ $sub_key } ] };
        }
        # Get the hash array slice of the values
        @{$sub_record}{ @value_fields } = @values[ @header_fields{@value_fields} ];
    }

    #__DUMP PARSED RESULTS AND STOP IF JUST CHECKING
        Process($record);

    #__RETURN DATA STRUCTS
    return ($record);
}

sub Process {
   my ($input) = @_;
   my ($output) = {};

   foreach my $hugo (keys (%{$input})) {
      foreach my $sample (keys (%{$input->{$hugo}})) {
         foreach my $line_num (keys (%{$input->{$hugo}->{$sample}})) {
            $output->{$hugo}->{$sample}->{$line_num} =
                  $input->{$hugo}->{$sample}->{$line_num};
            my ($residue1, $res_start, $residue2, $res_stop, $new_residue) =
            AA_Check(
               $input->{$hugo}->{$sample}->{$line_num}->{AA_CHANGE} , $line_num
            );

            $output->{$hugo}->{$sample}->{$line_num}->{res_start} = $res_start;
            $output->{$hugo}->{$sample}->{$line_num}->{res_stop} = $res_stop;
            $output->{$hugo}->{$sample}->{$line_num}->{residue1} = $residue1;
            $output->{$hugo}->{$sample}->{$line_num}->{residue2} = $residue2;
            $output->{$hugo}->{$sample}->{$line_num}->{new_residue} = $new_residue;
         }
      }
   }

   #__DUMP PARSED RESULTS AND STOP IF JUST CHECKING
   return ($output);
}

sub AA_Check {
   my ($AminoAcidChange_string, $line_num) = @_;
     my ($residue1, $res_start, $residue2, $res_stop, $new_residue);

     unless (defined($AminoAcidChange_string)) {
         return ($residue1, $res_start, $residue2, $res_stop, $new_residue);
     }
     unless ($AminoAcidChange_string =~ m/^p\./) {
         return ($residue1, $res_start, $residue2, $res_stop, $new_residue);
     }
   #__FORMULATE ERROR STRING JUST IN CASE
   my $string = "'$AminoAcidChange_string' is not a valid AminoAcidChange";
   $string .= " on line $line_num" if defined $line_num;

   #__VALIDATE
     $AminoAcidChange_string =~ s/^p\.//x;
     if ($AminoAcidChange_string =~ /^ (\D+) (\d+) _ (\D+) (\d+) (.*) $/x ) {
         ($residue1, $res_start, $residue2, $res_stop, $new_residue) =
             ($1, $2, $3, $4, $5);
     } elsif ($AminoAcidChange_string =~ /^ (\D+) (\d+) (\D+) (.*) $/x ) {
         ($residue1, $res_start, $residue2, $new_residue) =
             ($1, $2, $3, $4);
         $res_stop = $res_start;
     } elsif ($AminoAcidChange_string =~ /^ (\d+) (.*) $/x ) {
         ($res_start, $residue2, $res_stop, $new_residue) =
             ($1, $2);
         $residue1 = '*';
         $res_stop = $res_start;
         $new_residue = $residue2;
     } elsif ($AminoAcidChange_string =~ /^ (\D+) (\d+) $/x ) {
         ($residue1, $res_start) =
             ($1, $2);
         $new_residue = ' ';
         $residue2 = ' ';
         $res_stop = $res_start;
     }

     if (defined($new_residue)) {
         $new_residue =~ s/^ > //x;
     }
     $new_residue ||= '';

     return ($residue1, $res_start, $residue2, $res_stop, $new_residue);
}



