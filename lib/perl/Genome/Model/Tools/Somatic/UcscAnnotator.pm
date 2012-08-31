#review tmooney
#Quite the mix of variable naming conventions and generally inscrutable activities


package Genome::Model::Tools::Somatic::UcscAnnotator;

use strict;
use warnings;
use Genome;
use Genome::Utility::HugoGene::HugoGeneMethods;
use Carp;
use IO::File;

class Genome::Model::Tools::Somatic::UcscAnnotator{
    is => 'Command',
    has => [
    input_file => {
        is  => 'String',
        is_input => 1,
        doc => 'The input file of variants to be annotated',
    },
    output_file => {
        is => 'Text',
        is_input => 1,
        is_output => 1,
        doc => "Store annotation in the specified file"
    },
    unannotated_file => {
        is => 'Text',
        is_input => 1,
        is_optional => 1,
        doc => "File of sites unable to be annotated",
        default => 'ucsc_unannotated_variants',
    },
    skip => {
        is => 'Boolean',
        default => '0',
        is_input => 1,
        is_optional => 1,
        doc => "If set to true... this will do nothing! Fairly useless, except this is necessary for workflow.",
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
    "runs ucsc annotation on some variants",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic ucsc-annotator...    
EOS
}

sub help_detail {                           
    return <<EOS 
runs ucsc annotation on some variants
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip) {
        $self->status_message("Skipping execution: Skip flag set");
        return 1;
    }
    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $input_file = $self->input_file;
    my $output_file = $self->output_file;

    # Default to the same name as output file with an extension if not specifically provided
    my $unannotated_file = $self->unannotated_file;
    $unannotated_file ||= $self->output_file . ".unannotated";

    open(UNANNOTATED, "> $unannotated_file") || die "Could not open '$unannotated_file': $!";

    open(OUT, "> $output_file") || die "Could not open '$output_file': $!";

    # Make look-up hash for UCSC gene name to Hugo gene name
    my $hugo_gene_ref = Genome::Utility::HugoGene::HugoGeneMethods::makeHugoGeneObjects();
    my %ucsc_to_hugo;
    foreach my $hugo ( keys %{$hugo_gene_ref} ) {
        my $ucsc_name = $$hugo_gene_ref{$hugo}->ucsc();
        if ( !defined $ucsc_name || $ucsc_name eq "" ) { next; }
        if ( defined $ucsc_to_hugo{$ucsc_name} ) {
            print "WARNING: $ucsc_name => $hugo and $ucsc_name => $ucsc_to_hugo{$ucsc_name} \n";
        }
        $ucsc_to_hugo{$ucsc_name} = $hugo;
    }


    use DBI;
    my $db = "ucsc";
    my $user = "mgg_admin";
    my $password = "c\@nc3r"; 
    my $dataBase = "DBI:mysql:$db:mysql2";
    my $dbh = DBI->connect($dataBase, $user, $password) ||
    die "ERROR: Could not connect to database: $! \n";


    # The tables have slightly different entries.  The repeatmask tables 'chr$Chr_rmsk', 
    # knowGene, tfbsConsSites and recombRate tables are special cases.

    # use 'varType' for table 'dgv'  
    
    my ( @tables_with_scores, @tables_with_names, @tables_with_neither, $gene_table_query, @tables_with_variation_type,
        %chr_to_repeatmask, %queries_with_description, %queries_without_description,
        $table, $repeat_mask_query, $scores_table_query, $start, $stop, $description, $chr_start, $chr_stop,
        $names_table_query, $no_score_name_table_query, $exon_count, $exon_starts, $exon_ends, %chr_to_self_chain);


    my %table_with_no_bin = ( cnpSebat2 => 1, cpgIslandExt => 1, gad => 1, knownGene => 1, recombRate => 1);   #keep track of which tables we cannot use binning on to speed up queries

    @tables_with_scores = qw ( delConrad2  eponine firstEF genomicSuperDups phastConsElements17way phastConsElements28way polyaDb polyaPredict simpleRepeat switchDbTss targetScanS  vistaEnhancers wgEncodeGisChipPet wgEncodeGisChipPetHes3H3K4me3 wgEncodeGisChipPetMycP493 wgEncodeGisChipPetStat1Gif wgEncodeGisChipPetStat1NoGif);

    @tables_with_variation_type = qw ( cnpLocke cnpSharp2 );  

    @tables_with_names = qw (  cnpSebat2   cpgIslandExt gad  microsat );

    @tables_with_neither = qw(cnpTuzun cnpRedon cnpIafrate2 encodeUViennaRnaz exaptedRepeats delHinds2 delMccarroll uppsalaChipH3acSignal uppsalaChipUsf1Signal uppsalaChipUsf2Signal wgEncodeUcsdNgTaf1Signal wgEncodeUcsdNgTaf1ValidH3K4me wgEncodeUcsdNgTaf1ValidH3ac wgEncodeUcsdNgTaf1ValidRnap wgEncodeUcsdNgTaf1ValidTaf oreganno regPotential7X laminB1 );


    # Query for gene table.  Call with $chr, $start, $end
    $gene_table_query = "SELECT name, txStart, txEnd, exonCount, exonStarts, exonEnds
    FROM knownGene
    WHERE chrom = ? && txEnd >= ? && txStart <= ?
    ORDER BY txStart";
    my $gene_statement = $dbh->prepare($gene_table_query) ||
    die "Could not prepare statement '$gene_table_query': $DBI::errstr \n";
    $DB::single =1 ;
    
    # Define available chromosomes for repeatmasker to use
    my @available_chromosomes = (1..22, "X","Y");

    #  Query for repeatmask tables.  Call with start, end of region to check
    #  Call with $chr_to_repeatmaskStatements{$chr}->execute($start, $stop)
    foreach my $chr (@available_chromosomes) {
        $table = "chr$chr"."_rmsk";
        $repeat_mask_query = "SELECT repFamily, genoStart, genoEnd 
        FROM $table
        WHERE genoEnd >= ? && genoStart <= ?   
        AND %s
        ORDER BY genoStart";

        $chr_to_repeatmask{$chr} = $repeat_mask_query;
    }
    #   Query for selfChain tables
    foreach my $chr (@available_chromosomes) {
        my $table = "chr$chr"."_chainSelf";
        my $self_chain_query = "SELECT normScore, tStart, tEnd 
        FROM $table
        WHERE tEnd >= ? && tStart <= ?   
        AND %s
        ORDER BY tStart";

        $chr_to_self_chain{$chr} = $self_chain_query;    
    }
    
    # Query with tables that have scores; Want to display the score in output
    # To handle binning we will not include it in the string
    foreach $table (@tables_with_scores) {
        # Query for tables that have score
        $scores_table_query = "SELECT score, chromStart, chromEnd
        FROM $table
        WHERE chrom = ? && chromEnd >= ? && chromStart <= ? 
        AND %s
        ORDER BY chromStart";

        $queries_with_description{$table} = $scores_table_query;
    }


    # Query for tfbsConsSites table, which uses zScore rather than score
    my $tfbs_cons_query = "SELECT zScore, chromStart, chromEnd
    FROM tfbsConsSites 
    WHERE chrom = ? && chromEnd >= ? && chromStart <= ? 
    AND %s
    ORDER BY chromStart";
    $queries_with_description{"tfbsConsSites"} = $tfbs_cons_query;  


    # Query for tables that do not have a score but have a 'name' field
    # Want to display the name instead of the score
    foreach $table (@tables_with_names) {
        $names_table_query = "SELECT name, chromStart, chromEnd
        FROM $table
        WHERE chrom = ? && chromEnd >= ? && chromStart <= ? 
        AND %s
        ORDER BY chromStart";

        $queries_with_description{$table} = $names_table_query;
    }

    # These tables have a variation type which is displayed
    foreach $table (@tables_with_variation_type) {
        my $query = "SELECT variationType, chromStart, chromEnd
        FROM $table
        WHERE chrom = ? && chromEnd >= ? && chromStart <= ? 
        AND %s
        ORDER BY chromStart"; 

        $queries_with_description{$table} = $query;
    }

    # The table 'dgv' has 'varType' as column name rather than variationType
    my $query = "SELECT varType, chromStart, chromEnd
    FROM dgv
    WHERE chrom = ? && chromEnd >= ? && chromStart <= ?
    AND %s
    ORDER BY chromStart"; 
    $queries_with_description{"dgv"} = $query;


    # The other tables do not have anything to display.  Just annotate with a '+'
    # if appropriate
    foreach $table (@tables_with_neither) {
        #Query for tables that do not have score or name fields
        $no_score_name_table_query = "SELECT chromStart, chromEnd
        FROM $table
        WHERE chrom = ? && chromEnd >= ? && chromStart <= ? 
        AND %s
        ORDER BY chromStart";

        $queries_without_description{$table} = $no_score_name_table_query;
    } 


    # Need special sub for recombRate.  Genome is divided into 1000000 bp windows.  Only get
    # rate in one window.  There are three values (avg, male, female) for three maps (Decode, Marshfield, 
    # Genethon).  Not all maps have a rate.  Ignore values of '0'
    #
    # This does not support bins. Leaving alonge
    my $recombination_query = "SELECT decodeAvg, marshfieldAvg, genethonAvg
    FROM recombRate
    WHERE chrom = ? && chromEnd >= ? && chromStart <= ?";
    my $recombination_statement = $dbh->prepare($recombination_query) ||
    die "Could not prepare statement '$recombination_query': $DBI::errstr \n";

    # First print out the headers. This order MUST match the order in which the
    # statements are executed.
    print OUT "chr\tstart\tstop\tdecode,marshfield,genethon\trepeatMasker\tselfChain";
    #print OUT "chr\tstart\tstop\trepeatMasker";
    foreach $table ( sort keys %queries_with_description ) { print OUT "\t$table"; }
    foreach $table (sort keys %queries_without_description) { print OUT "\t$table"; }
    print OUT "\tknownGenes\tHUGO symbol\n";


    open(IN, "< $input_file") ||
    die "Could not open '$input_file': $!";
    my @entire_file = <IN>;
    close IN;
    my ($got_entry, %description_list, );
    
    foreach my $line (@entire_file) {
        chomp $line;
        my ($Chr, $start, $stop) = split /\s+/, $line;

        #check if we can annotate this site
        unless(grep {$Chr eq $_} @available_chromosomes) {
            warn "Unable to annotate $Chr\t$start\t$stop. Chromosome unavailable for annotation\n";
            print UNANNOTATED "$Chr\t$start\t$stop\n"; 
            next;
        }

        print OUT "$Chr\t$start\t$stop\t"; 
        $start = $start - 1; #change to 0 based
        $stop = $stop - 1; #change to 0 based
        $got_entry = 0; %description_list = ();


        # Recombination query
        $recombination_statement->execute("chr$Chr", $start, $start) ||
        die "Could not execute statement for repeat masker table with ($start, $stop): $DBI::errstr \n";
        while ( my ($decode, $marsh, $genethon) = $recombination_statement->fetchrow_array() ) {
            if ( $decode == 0 ) { $decode = "-"; }
            if ( $marsh == 0 ) { $marsh = "-"; }
            if ( $genethon == 0 ) { $genethon = "-"; }
            print OUT "$decode $marsh $genethon ";
            $got_entry = 1;
        }
        if ( !$got_entry ) { print OUT "- - -"; }
        print OUT "\t";
        my $bin_string = $self->bin_query_string($start,$stop);

        $got_entry = 0; 
        # Repeatmasker query
        my $repeat_mask_query_string = sprintf($chr_to_repeatmask{$Chr},$bin_string); 
        my $repeat_statement = $dbh->prepare_cached($repeat_mask_query_string) ||
        die "Could not prepare statement '$repeat_mask_query_string': $DBI::errstr \n";

        $repeat_statement->execute($start,$stop);
       


        #$chr_to_repeatmaskStatements{$Chr}->execute($start, $stop) ||
        #die("Could not execute statement for repeat masker table with ($start, $stop) for chromosome $Chr : $DBI::errstr \n");
#        while ( ($description, $chr_start, $chr_stop) =  $chr_to_repeatmaskStatements{$Chr}->fetchrow_array() ) {
        while ( ($description, $chr_start, $chr_stop) =  $repeat_statement->fetchrow_array() ) {
        
            $description_list{$description} = 1;
            $got_entry = 1;
        }
        if ( $got_entry ) { 
            foreach (keys %description_list) { print OUT "$_ "; }
        } else {
            print OUT  "-"; 
        }
        print OUT  "\t";

        # selfChain query
        %description_list = ();
        $got_entry = 0;

        
        my $self_chain_query_string = sprintf($chr_to_self_chain{$Chr}, $bin_string); 
        
        my $self_chain_statement = $dbh->prepare_cached($self_chain_query_string) ||
        die "Could not prepare statement '$self_chain_query_string': $DBI::errstr \n";

        $self_chain_statement->execute($start,$stop);

        #$chr_to_self_chainStatements{$Chr}->execute($start, $stop) || die "Could not execute statement for selfChain table with ($start, $stop): $DBI::errstr \n";
        while ( ($description, $chr_start, $chr_stop) =  $self_chain_statement->fetchrow_array() ) {
            $description_list{$description} = 1;
            $got_entry = 1;
        }
        if ( $got_entry ) { 
            foreach (keys %description_list) { print OUT "$_ "; }
        } else {
            print OUT  "-"; 
        }
        print OUT  "\t";

        # Tables that have a description or score
        foreach $table (sort keys %queries_with_description) {
            $got_entry = 0; %description_list = ();
            my $query_string = sprintf($queries_with_description{$table}, exists($table_with_no_bin{$table}) ? "1" : $bin_string);
            my $statement = $dbh->prepare_cached($query_string) ||
            die "Could not prepare statement '$query_string': $DBI::errstr \n";

            $statement->execute("chr$Chr", $start, $stop) ||
            die "Could not execute statement for table '$table' with ($Chr, $start, $stop): $DBI::errstr \n";
            while ( ($description, $chr_start, $chr_stop) =  $statement->fetchrow_array() ) {
                if ( $description eq "" ) { next; }
                $description_list{$description} = 1;
                $got_entry = 1;
            }
            if ( $got_entry ) { 
                foreach (keys %description_list) { print OUT "$_ "; }
            } else {
                print OUT "-"; 
            }
            print OUT "\t";
        }


        # Tables that are annotated with either a '+' or '-'
        foreach $table (sort keys %queries_without_description) {
            $got_entry = 0; 
            my $query_string = sprintf($queries_without_description{$table}, exists($table_with_no_bin{$table}) ? "1" : $bin_string);
            my $statement = $dbh->prepare_cached($query_string) ||
            die "Could not prepare statement '$query_string': $DBI::errstr \n";
            $statement->execute("chr$Chr", $start, $stop) ||
            die "Could not execute statement for table '$table' with ($Chr, $start, $stop): $DBI::errstr \n";
            while ( ($chr_start, $chr_stop) =  $statement->fetchrow_array() ) {
                # Only print out one '+' even if there are several entries that overlap
                if ( !$got_entry ) { print OUT "+"; }
                $got_entry = 1;
            }
            if ( !$got_entry ) { print OUT  "-"; }
            print OUT "\t";
        }   

        # Genes in regions
        $got_entry = 0;
        my (@hugo_names, $exon);
        $gene_statement->execute("chr$Chr", $start, $stop) ||
        die "Could not execute statement for table '$table' with ($Chr, $start, $stop): $DBI::errstr \n";
        while ( ($description, $chr_start, $chr_stop, $exon_count, $exon_starts, $exon_ends) = $gene_statement->fetchrow_array() ) {
            print OUT "$description ";
            if ( defined $ucsc_to_hugo{$description} ) { push @hugo_names, $ucsc_to_hugo{$description}; }
            $got_entry = 1;
            $exon = $self->regionOverlapsExons($start, $stop, $exon_count, $exon_starts, $exon_ends);
            if ( $exon ) { print OUT " $exon "
            }
        }
        if ( !$got_entry ) { print OUT "-\t-"; }
        print OUT "\t";
        if ( scalar(@hugo_names) >= 1 ) { print OUT "@hugo_names"; } else { print OUT "-"; }
        print OUT "\n";

    }
    
    $dbh->disconnect();
    close UNANNOTATED;
    close OUT;
}

# Returns "exonNumber start stop"  if one of the region overlaps one of the exons
# 0 if not
sub regionOverlapsExons {
    my $self = shift;
    my ($start, $end, $exon_count, $exon_starts, $exon_ends) = @_;

    my ( @starts, @ends, );
    @starts =  split /,/, $exon_starts;
    @ends = split /,/, $exon_ends;
    ( scalar(@starts) == $exon_count && scalar(@ends) == $exon_count ) ||
    confess "Did not get expected number of exons";
    for ( my $i = 0; $i <= $#starts; $i++ ) {
        if ( $ends[$i] >= $start && $starts[$i] <= $end ) {
            my $num = $i + 1; 
            my $exon = " Exon $num $starts[$i] $ends[$i] ";
            return $exon;
        }
    }
    return 0;
}

sub calculate_bin_from_range {
    my ($self, $start, $end) = @_;
    #This code derived from C code from http://genomewiki.ucsc.edu/index.php/Bin_indexing_system

    #This file is copyright 2002 Jim Kent, but license is hereby
    #granted for all use - public, private or commercial. */

    # add one new level to get coverage past chrom sizes of 512 Mb
    #      effective limit is now the size of an integer since chrom start
    #      and end coordinates are always being used in int's == 2Gb-1
    my @bin_offsets_extended = (4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0);
    my $_bin_first_shift =  17;       # How much to shift to get to finest bin.
    my $_bin_next_shift = 3;          # How much to shift to get to next larger bin.
    my $_bin_offset_old_to_extended = 4681;

    # Given start,end in chromosome coordinates assign it
    # a bin.   There's a bin for each 128k segment, for each
    # 1M segment, for each 8M segment, for each 64M segment,
    # for each 512M segment, and one top level bin for 4Gb.
    #      Note, since start and end are int's, the practical limit
    #      is up to 2Gb-1, and thus, only four result bins on the second
    #      level.
    # A range goes into the smallest bin it will fit in. */
    my ($start_bin, $end_bin) = ($start, $end-1);
    $start_bin >>= $_bin_first_shift;
    $end_bin >>= $_bin_first_shift;
    for (my $i = 0; $i < scalar(@bin_offsets_extended); $i++) {
        if ($start_bin == $end_bin) {
            return $_bin_offset_old_to_extended + $bin_offsets_extended[$i] + $start_bin;
            $start_bin >>= $_bin_next_shift;
            $end_bin >>= $_bin_next_shift;
        }
        $self->error_message(sprintf("start %d, end %d out of range in calculate_bin_from_range (max is 2Gb)", $start, $end));
        return 0;
    }
}

sub bin_query_string {
    my ($self, $start, $end) = @_;  #start and end should probably be 0 based as UCSC is 0 based
    #This code taken from function in kent src tree of UCSC called static void hAddBinToQueryStandard(char *binField, int start, int end, struct dyString *query, boolean selfContained)
    #Found this online at http://code.google.com/p/genomancer/source/browse/trunk/poka/src/genomancer/ucsc/das2/BinRange.java?spec=svn66&r=66 and haven't bothered to look in the actual source
    my ($b_first_shift, $_b_next_shift) = (17,3);
    my $start_bin = ($start>>$b_first_shift);
    my $end_bin = (($end)>>$b_first_shift);#TODO figure out if the -1 is necessary
    my @binOffsets = ( 512+64+8+1, 64+8+1, 8+1, 1, 0); #Not using the extended binning scheme...
    my $_bin_offset_old_to_extended = 4681;

    my $bin_query_string = "(";
    for (my $i = 0; $i < scalar(@binOffsets); ++$i) {
        my $offset = $binOffsets[$i];
        if ($i != 0) {
            $bin_query_string .=  " or ";
        }
        if ($start_bin == $end_bin) {
            #assuming the binField is actually bin in all cases (may not be true?)
            $bin_query_string .= sprintf("%s=%u", "bin", $start_bin + $offset);
        }
        else {
            $bin_query_string .= sprintf("( %s>=%u and %s<=%u )", "bin", $start_bin + $offset, "bin", $end_bin + $offset);
        }
        $start_bin >>= $_b_next_shift;
        $end_bin >>= $_b_next_shift;
    }
    $bin_query_string .= sprintf(" or %s=%u )", "bin", $_bin_offset_old_to_extended);
    return $bin_query_string;
}


1;
