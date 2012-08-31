package Genome::Model::Tools::Snp::Evaluation;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bio::DB::Fasta;

class Genome::Model::Tools::Snp::Evaluation {
    is => 'Command',
    has => [
    snp_file => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "input file of snps in the experimental metrics format",
    },
    original_snp_file => 
    { 
        type => 'String',
        is_optional => 1,
        doc => "Input file of snps in the experimental metrics format. Used to compare the metrics of a filtered snp file to it's parent. Required for calculation of database specificity and sensitivity",
    },
    gold_snp_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "input file of snp locations and calls from the intersection of the affy and illumina platforms",
    },
    dbsnp_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "file of dbsnp/watson/venter status for all locations in the model",
    },
    print =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Print the results in human readable format",
        default => 0,
    },        
    exclude_y =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Don't consider SNPs present on the Y chromosome",
        default => 0,
    },        
    total_snp_positions =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_concordant_positions =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_gold_heterozygote_snps =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_gold_homozygous_ref_positions =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_gold_het_concordant_snps =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_gold_hom_ref_concordant_snps =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_database_germline =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_database_wildtype =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_database_somatic =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_database_present_wildtype =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_database_present_germline =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    total_database_present_somatic =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 0,
    },        
    REFDIR =>
    {
        #This is the path to the reference sequence used for aligning the model
        type => 'String',
        is_optional => 0,
        default => "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c",
    },        
    #The following are bitmasks to determine novel vs non-novel status
    #They are stored as strings to hopefully reduce memory slightly since we are hashing the whole shebang
    DBSNP_MASK =>
    {
        type => 'Bitmask',
        is_optional => 0,
        default => "100",
    },        
    WATSON_MASK =>
    {
        type => 'Bitmask',
        is_optional => 0,
        default => "010",
    },        
    VENTER_MASK =>
    {
        type => 'Bitmask',
        is_optional => 0,
        default => "001",
    },        
    ]
};


sub execute {
    my $self=shift;

    #Check on the file names
    unless(-f $self->snp_file) {
        $self->error_message("Snps file is not a file: " . $self->snp_file);
        return;
    }
    if(defined($self->original_snp_file) && !(-f $self->original_snp_file)) {
        $self->error_message("Original SNP file is not a file: " . $self->original_snp_file);
        return;
    }
    unless(-f $self->gold_snp_file) {
        $self->error_message("Gold snp file is not a file: " . $self->gold_snp_file);
        return;
    }
    unless(-f $self->dbsnp_file) {
        $self->error_message("dbSNP file is not a file: " . $self->dbsnp_file);
        return;
    }

    if(defined($self->original_snp_file) && !(-f $self->dbsnp_file)) {
        $self->error_message("dbSNP file is not a file: " . $self->dbsnp_file);
        return;
    }

    #Check and open filehandles
    my $snp_fh=IO::File->new($self->snp_file);
    unless($snp_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->snp_file );
        return;
    }
    my $original_fh;
    if(defined($self->original_snp_file)) {
        $original_fh=IO::File->new($self->original_snp_file);
        unless($original_fh) {
            $self->error_message("Failed to open filehandle for: " .  $self->original_snp_file );
            return;
        }
    }
    my $gold_fh=IO::File->new($self->gold_snp_file);
    unless($gold_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->gold_snp_file );
        return;
    }
    my $dbsnp_fh=IO::File->new($self->dbsnp_file);
    unless($dbsnp_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->dbsnp_file);
        return;
    }

    my $dbsnp_hash_ref = $self->create_dbSNP_hash($dbsnp_fh);
    close($dbsnp_fh);

    unless(defined($dbsnp_hash_ref)) {
        $self->error_message("Fatal error creating dbSNP hash");
        return;
    }


    my ($gold_het_hash_ref, $gold_ref_hash_ref) = $self->create_gold_snp_hashes($gold_fh);
    close($gold_fh);
    
    unless(defined($gold_het_hash_ref)) {
        $self->error_message("Fatal error creating Gold SNP hash");
        return;
    }

    my ($validation_db_hash_ref) = $self->create_db_validation_hash();
    #calculate the metrics on the primary file

    #set up temporary values
    my ($total_snp_positions, 
        $total_concordant_positions, 
        $total_gold_het_concordant_snps, 
        $total_gold_hom_ref_concordant_snps, 
        $total_database_present_wildtype, 
        $total_database_present_germline, 
        $total_database_present_somatic) = 
    $self->calculate_metrics($snp_fh,$dbsnp_hash_ref,$gold_het_hash_ref,$gold_ref_hash_ref,$validation_db_hash_ref);

    #initialize if undefined
    $total_snp_positions ||= 0; 
    $total_concordant_positions ||= 0; 
    $total_gold_het_concordant_snps ||= 0; 
    $total_gold_hom_ref_concordant_snps ||= 0; 
    $total_database_present_wildtype ||= 0; 
    $total_database_present_germline ||= 0; 
    $total_database_present_somatic ||= 0; 

    
    #Set all the class variables for the primary file
    $self->total_concordant_positions($total_concordant_positions);
    $self->total_snp_positions($total_snp_positions);
    $self->total_gold_hom_ref_concordant_snps($total_gold_hom_ref_concordant_snps);
    $self->total_gold_het_concordant_snps($total_gold_het_concordant_snps);

    $self->total_database_present_somatic($total_database_present_somatic);
    $self->total_database_present_wildtype($total_database_present_wildtype);
    $self->total_database_present_germline($total_database_present_germline);

    if($self->print) {
        printf("Evaluation of the SNP File\n");
        $self->print_gold_and_dbsnp($total_concordant_positions,$total_snp_positions,$total_gold_het_concordant_snps, $self->total_gold_heterozygote_snps, $total_gold_hom_ref_concordant_snps, $self->total_gold_homozygous_ref_positions);
        printf("There were %i validated somatic SNVs detected out of %i total validated somatic SNVs in the database\n", $self->total_database_present_somatic,$self->total_database_somatic);
    }
    #Calculate metrics on the original file
    if(defined($original_fh)) {

        my ($org_total_snp_positions, 
            $org_total_concordant_positions, 
            $org_total_gold_het_concordant_snps, 
            $org_total_gold_hom_ref_concordant_snps, 
            $org_total_database_present_wildtype, 
            $org_total_database_present_germline, 
            $org_total_database_present_somatic) = 
        $self->calculate_metrics($original_fh, $dbsnp_hash_ref,$gold_het_hash_ref,$gold_ref_hash_ref,$validation_db_hash_ref); 

        #set class variables for sensitivity and specificity
        $self->total_database_somatic( $org_total_database_present_somatic );
        $self->total_database_germline( $org_total_database_present_germline );
        $self->total_database_wildtype( $org_total_database_present_wildtype );


        if($self->print) {
            printf("\nEvaluation of the Original SNP File\n");
            $self->print_gold_and_dbsnp($org_total_concordant_positions,$org_total_snp_positions,$org_total_gold_het_concordant_snps, $self->total_gold_heterozygote_snps, $org_total_gold_hom_ref_concordant_snps, $self->total_gold_homozygous_ref_positions);

            #print sensitivity and specificity
            printf("\nEvaluation of the SNP File based on database SNPs present in the Original SNP file\n");

            printf("There were %i validated SNVs detected out of %i total validated SNVs in the original file for a sensitivity of %0.2f%%\n", $self->total_database_present_somatic + $self->total_database_present_germline, $org_total_database_present_somatic + $self->total_database_germline,(100 * ($self->total_database_present_somatic+$self->total_database_present_germline) / ($self->total_database_somatic+$self->total_database_germline)));

            printf("There were %i validated false positive SNVs detected out of %i total validated false positive SNVs in the database for a specificity of %0.2f%%\n", $self->total_database_present_wildtype,$self->total_database_wildtype,100 * (1 - $self->total_database_present_wildtype / $self->total_database_wildtype));
        }
    }

    return 1;
}

    


1;

sub help_brief {
    "Calculates the dbSNP concordance and FNR/FPR of a snp file"
}

sub calculate_metrics {
    my ($self,$snp_fh,$dbsnp_hash_ref, $gold_het_hash_ref, $gold_ref_hash_ref, $validation_db_hash_ref) = @_;

    #set up temporary values
    my $total_gold_het_concordant_snps;
    my $total_snp_positions;
    my $total_gold_hom_ref_concordant_snps;
    my $total_concordant_positions;

    my $total_database_present_germline;
    my $total_database_present_wildtype;
    my $total_database_present_somatic;

    #skip header
    $snp_fh->getline;

    my %seen_at; #hash to store the locations we have already evaluated. Necessary because experimental metrics splits bi-allelic het sites.
    while(my $line = $snp_fh->getline) {
        #$DB::single = 1;
        next if($line =~ /chromosome/); #Skip experimental metrics headers if included because of say cat
        my ($chr,$pos,$ref,$var) = split /,\s*|\s+/, $line; #handle both types of experimental metrics formats,and SNP files

        if($self->exclude_y) {
            next if($chr eq 'Y'); #female patient these are BS
        }

        next if($ref eq ' ' || $ref eq '' || $ref eq 'N'); #skip 'SNPs' where the reference is N or non-existent

        #Check db status
        if(exists($validation_db_hash_ref->{$chr}{$pos})) { 
            if( $validation_db_hash_ref->{$chr}{$pos} eq 'G') {
                $total_database_present_germline++;
            }
            if( $validation_db_hash_ref->{$chr}{$pos} eq 'WT') {
                $total_database_present_wildtype++;
            }
            if( $validation_db_hash_ref->{$chr}{$pos} eq 'S') {
                $total_database_present_somatic++;
            }

        }


        #If there is a het SNP site, it is possible that there are two lines in the metrics file
        #corresponding to that position. We want to check both but only count the position once
        #Check that the genotypes match
        if(exists($gold_het_hash_ref->{$chr}{$pos})) {
            if( (@{$gold_het_hash_ref->{$chr}{$pos}}[0] eq $ref && 
                    @{$gold_het_hash_ref->{$chr}{$pos}}[1] eq $var)
                || (@{$gold_het_hash_ref->{$chr}{$pos}}[0] eq $var &&
                    @{$gold_het_hash_ref->{$chr}{$pos}}[1] eq $ref)
            ) {
                #Any het site should be counted. Bi-allelic het sites should be counted if either matches
                $total_gold_het_concordant_snps++;
            }

        }
        next if(exists($seen_at{$chr}{$pos})); #only count each position once
        $seen_at{$chr}{$pos} = 1;
        $total_snp_positions++;
        if(exists($gold_ref_hash_ref->{$chr}{$pos}) &&
            @{$gold_ref_hash_ref->{$chr}{$pos}}[0] eq $ref) {
            #){
            $total_gold_hom_ref_concordant_snps++;
        }
        if(exists($dbsnp_hash_ref->{$chr}{$pos})) {
            if(int($dbsnp_hash_ref->{$chr}{$pos} & $self->DBSNP_MASK)) {
                $total_concordant_positions++;
            }
        }
    }
    return ($total_snp_positions, $total_concordant_positions, $total_gold_het_concordant_snps, $total_gold_hom_ref_concordant_snps, $total_database_present_wildtype, $total_database_present_germline, $total_database_present_somatic);

}

sub print_gold_and_dbsnp {
    my ($self, $concordant_positions, $snp_positions,  $gold_het_concordant, $gold_het_total, $gold_ref_concordant, $gold_ref_total) = @_; 
    printf("There were %i positions in dbSNP out of %i total SNP positions for a dbSNP concordance of %0.2f%%\n", $concordant_positions, $snp_positions, ($concordant_positions/$snp_positions)*100);

    printf("There were %i heterozygous Gold SNPs detected out of %i total for a false negative rate of %0.2f%%\n", $gold_het_concordant, $gold_het_total, (1 - ($gold_het_concordant/$gold_het_total)) * 100);

    printf("There were %i homozygous reference Gold SNPs detected out of %i total for a false positive rate of %0.2f%%\n", $gold_ref_concordant, $gold_ref_total, ($gold_ref_concordant/$gold_ref_total) * 100);
}

#Functions to create hashes of data

#Not sure about the memory requirements of this
#It would be more efficient, memory-wise to stream these directly from the file on the fly
#But more complicated

#Create hash of strings indicating dbSNP status of each position
sub create_dbSNP_hash {
    my ($self, $dbsnp_fh) = @_;

    #read in header
    my $header = $dbsnp_fh->getline;

    my %dbsnp_at;
    while(my $line = $dbsnp_fh->getline) {
        chomp $line;
        my ($chr, $start, $end, $dbsnp, $watson, $venter) = split /\t/, $line;
        #add in compatibility for non-watson/venter including files
        #TODO Fix this when we get watson and venter in OLAP
        $watson ||= 0;
        $venter ||= 0;
        if(!exists($dbsnp_at{$chr}{$start})) {
            $dbsnp_at{$chr}{$start} = $dbsnp.$watson.$venter; #represent as a binary string
        }
        else {
            $self->error_message( "Duplicate entry in dbSNP/Watson/Venter file");
            return;
        }
    }
    return \%dbsnp_at;
}

#Create hashes of gold SNPs

sub create_gold_snp_hashes {
    my ($self,$gold_fh) = @_;

    #create db instance to check the ref seq
    my $refdb = Bio::DB::Fasta->new($self->REFDIR);

    #create temporary counter variables
    my $total_gold_homozygous_ref_positions;
    my $total_gold_heterozygote_snps;
    
    my %heterozygous_snp_at;
    my %reference_snp_at;
    #Read through the gold snp file. Any heterozygous SNPs that are found can be used for false negative rate calculation
    #Any homozygous ref that can be found can be used for false positive calculation

    #there is no header on this file
    #Format is tab separated
    #Chr\tPos\tPos\tAllele1\tAllele2\tPlatform1_Allele1\tPlatform1_Allele2\tPlatform2_Allele1\tPlatform2_Allele2
    while(my $line = $gold_fh->getline) {
        chomp $line;
        my ($chr, $pos, $pos2, $allele1, $allele2, $allele1_type1,$allele2_type1, $allele1_type2, $allele2_type2) = split /\t/, $line;

        my $ref_a= $refdb->seq($chr, $pos => $pos2);  
        chomp($ref_a);
        $ref_a=uc($ref_a);

        if($allele1 eq $allele2) {
            #homozygous
            if($allele1_type1 eq $allele1_type2) {
                #Check that the platform agree is internally consistent
                if($allele1_type1 ne $allele2_type1 || $allele1_type2 ne $allele2_type2) {
                    $self->error_message("Inconsistent types within a platform on a homozygous SNP at " . $gold_fh->input_line_number);
                    next;
                }
                #            $total_gold++;
                if($allele1_type1 eq 'ref') {
                    if($allele1 eq $ref_a) {
                        #it was in our reference as ref
                        $reference_snp_at{$chr}{$pos} = [$allele1,$allele2];
                        $total_gold_homozygous_ref_positions++;
                    }
                    else {
                        #Gold SNP reference base is altered
                        #So it is not actually a reference base
                        $self->error_message("Gold SNP reference doesn't match B36 reference sequence");
                    }
                }
                else {
                    #assuming assuming homozygous SNP at this position
                    #we aren't tracking these though
                    #if($allele1 ne $ref_a) {
                    #    $homozygous_snp_at{$chr}{$pos} = [$allele1,$allele2];
                    #    $total_hom++;
                    #}
                    #else {
                    #    $self->error_message("Gold SNP is listed as reference in B36 sequence");
                    #}
                }

            }
            else {
                #platforms disagree
                if($allele1_type1 ne $allele2_type1 || $allele1_type2 ne $allele2_type2) {
                    $self->error_message("Inconsistent types within a platform on a homozygous SNP at ".$gold_fh->input_line_number);
                    next;
                }
                #$total_gold++;
                #Check if the allele matches the reference base in B36
                if($allele1 ne $ref_a) {
                    #It's a SNP!
                    #$homozygous_snp_at{$chr}{$pos} = [$allele1,$allele2];
                    #$total_hom++;
                }
                else {
                    #it was in our reference as ref
                    $reference_snp_at{$chr}{$pos} = [$allele1,$allele2];
                    $total_gold_homozygous_ref_positions++;
                }

            }
        }
        else {
            #heterozygous site
            #
            #Check that the platforms agree
            if($allele1_type1 eq $allele1_type2 && $allele2_type1 eq $allele2_type2) {
                #het site
                #check that the allele is actually a snp
                if($allele1_type1 eq 'SNP' ) {
                    #$total_gold++;
                    $total_gold_heterozygote_snps++;
                    $heterozygous_snp_at{$chr}{$pos} = [$allele1,$allele2];
                }
                elsif($allele2_type1 eq 'SNP' ) {
                    #$total_gold++;
                    $total_gold_heterozygote_snps++;
                    $heterozygous_snp_at{$chr}{$pos} = [$allele2,$allele1];
                }
                else {
                    #something is up, neither alleles is labeled as SNP
                    $self->error_message("Supposedly heterozygous SNP not labeled as such at " . $gold_fh->input_line_number);
                }
                    
            }
            else {
                $self->error_message("Platforms disagree on reference at line ".$gold_fh->input_line_number);
            }
        }
    }
    #set the class variables
    $self->total_gold_heterozygote_snps($total_gold_heterozygote_snps);
    $self->total_gold_homozygous_ref_positions($total_gold_homozygous_ref_positions);

    return (\%heterozygous_snp_at,\%reference_snp_at);
}

#In my neverending quest to use as many hashes as possible...
#This creates a hash of the db results
#This shouldn't suck up too much memory unless we validate a whole lotta crap
sub create_db_validation_hash {
    my ($self) = @_;
    my @variants = Genome::VariantReviewDetail->get(variant_type => 'S');
    my %validated_status_for;

    my $total_database_wildtype;
    my $total_database_germline;
    my $total_database_somatics;

    my %rank = ('S' => 0,
        'G' => 1,
        'WT' => 2,
        'A' => 3,
        'V' => 3,
        'O' => 3,
        'LQ' => 4,
        'X' => 4,
        'NC' => 5,
    );

    foreach my $variant (@variants) {
        my $chr = $variant->chromosome;
        my $pos = $variant->begin_position;
        my $status = $variant->somatic_status;
        next if(!defined($chr) || !defined($pos));

        next if(!defined($status));

        #The following handles two levels of discrepancy
        #TODO make it recursive so it can handle arbitrary levels
        if(!exists($validated_status_for{$chr}{$pos})) {
            my ($level1_discrep) = $status =~ /^DISCREPANCY\((.*)/;
            my $level2_discrep;
            if(defined($level1_discrep)) {
                ($level2_discrep) = $level1_discrep =~ /^DISCREPANCY\((.*)/;
            }
            if(defined($level2_discrep)) {
                my @status = split /:/, $level2_discrep;
                map {s/^\s*(\S*)\s*$/$1/} @status;
                my $new_status = 'NC';
                foreach my $discr_status (@status) {
                    if(exists($rank{$discr_status}) && $rank{$discr_status} < $rank{$new_status} ) {
                        $new_status = $discr_status;
                    }
                }

                $validated_status_for{$chr}{$pos} = $new_status;

            }
            elsif(defined($level1_discrep)) {
                my @status = split /:/, $level1_discrep;
                map {s/^\s*(\S*)\s*$/$1/} @status;
                my $new_status = 'NC';
                foreach my $discr_status (@status) {
                    unless(defined($rank{$discr_status})) {
                        $self->error_message("Unexpected status: '$discr_status'"); #quotes for whitespace visualization
                    }
                    if($rank{$discr_status} < $rank{$new_status} ) {
                        $new_status = $discr_status;
                    }
                }
                $validated_status_for{$chr}{$pos} = $new_status;
            }
            else {
                $validated_status_for{$chr}{$pos} = $status;
            }
        }
        else {
            #check and see if it is an improvement over NC or LQ
            if(uc($status) eq 'G' || uc($status) eq 'WT' || uc($status) eq 'S') {
                $validated_status_for{$chr}{$pos} = $status;
            }
        }
        #At this point, assuming the database has only 1 entry per position, we have the final status
        if($validated_status_for{$chr}{$pos} eq 'G') {
            $total_database_germline++;
        }
        if($validated_status_for{$chr}{$pos} eq 'S') {
            $total_database_somatics++;
        }
        if($validated_status_for{$chr}{$pos} eq 'WT') {
            $total_database_wildtype++;
        }
    }
    #Set class variables
    #$self->total_database_wildtype($total_database_wildtype);
    $self->total_database_somatic($total_database_somatics);
    #$self->total_database_germline($total_database_germline);
    return \%validated_status_for;
}
