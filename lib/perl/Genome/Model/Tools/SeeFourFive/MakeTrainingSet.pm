package Genome::Model::Tools::SeeFourFive::MakeTrainingSet;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use List::Util 'shuffle';
use Bio::DB::Fasta;

class Genome::Model::Tools::SeeFourFive::MakeTrainingSet {
    is => 'Command',
    has => [
    snp_file => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "input file of known real snps with
        lots of data about them in columns",
    },
    gold_snp_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "input file of known false snps with
        lots of data about them in columns",
    },
    gold_set_size =>
    {
        type => 'Integer',
        is_optional => 1,
        doc => 'Number of Germline Gold SNPs to include in the training set. The rest will be in the test set',
        default => '200'
    },        
    name_file =>
    {
        type => 'String',
        is_optional => 1,
        doc => "config file detailing name and range 
        of columns in the two snp files for C4.5",
    },
    data_file =>
    {
        type => 'String',
        is_optional => 1,
        doc => "file of data which C4.5 will use to 
        attempt to construct good rules. If not specified
        a default subset will be created at runtime",
    },
    test_file =>
    {  
        type => 'String',
        is_optional =>1,
        doc => "file of data which C4.5 will use to
        test out the rules it decides on. if not specified,
        good/bad concat will be used",
    },
    exclude_y =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Don't consider SNPs present on the Y chromosome",
        default => 0,
    },        
    REFDIR =>
    {
        #This is the path to the reference sequence used for aligning the model
        type => 'String',
        is_optional => 0,
        default => "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c",
    },        
    use_database =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => 'Use values from the variant_review database as a training set',
        default => 1,
    },
    data_format =>
    {
        type => 'String',
        is_optional => 0,
        doc => 'Which Maq::Metrics::Dtr module to use',
        default => 'MaqOSixThree',

    },
    ]
};

#Pseudocode
#This module should now read in the metrics file, and a Gold SNP file
#It should then determine which metrics refer to Germline SNPs, and which refer to wild-type snps
#Germline = anything in the metrics file that corresponds to a known Gold SNP genotype, or anything in the db as Germline or Somatic
#WT = anything in the metrics file that corresponds to a known reference Gold SNP, or anything in the db as WT
#Since we expect there to be many more G SNPs than WT in the Gold Set (Maq is pretty good at what it does)
#then our final step is to pad a random set of G SNPs from the Gold Set into our training set
#We can leave the rest as the test set

sub execute {
    my $self=shift;
    unless(-f $self->snp_file) {
        $self->error_message("Snp file is not a file: " . $self->snp_file);
        return;
    }
    unless(-f $self->gold_snp_file) {
        $self->error_message("Gold snp file is not a file: " . $self->gold_snp_file);
        return;
    }
    my $snp_fh=IO::File->new($self->snp_file);
    my $gold_fh=IO::File->new($self->gold_snp_file);
    unless($snp_fh && $gold_fh) {
        $self->error_message("Failed to open filehandles for: " .  $self->snp_file . " and/or " . $self->gold_snp_file);
        return;
    }

    #lots of processing goes here... to both bad and good files
    #independify various columsn, add WT/G to end of each line
    #they should probably be written either to two new files or just dumped right into test and data files. 

    #Get object for making, and reading tree data
    my $type = $self->data_format;
    my $dtr = eval "Genome::Model::Tools::Maq::Metrics::Dtr::$type->create()";
    unless(defined($dtr)) {
        $self->error_message($@);
        return;
    }
    my $names = $dtr->names_file_string;

    my ($gold_het_href,$gold_hom_href,$gold_ref_href) = $self->create_gold_snp_hashes($gold_fh);
    my ($validation_href) = $self->use_database ? $self->make_validation_hash() : undef;
    my ($db_wildtype_aref,
        $db_germline_aref,
        $gold_wildtype_aref,
        $gold_germline_aref) = $self->make_data_arrays($snp_fh,$dtr, $gold_het_href,$gold_hom_href, $gold_ref_href, $validation_href);

    map {undef $_} $gold_het_href,$gold_hom_href,$gold_ref_href,$validation_href; #hopefully free up some mem

    my @data = shuffle @$gold_germline_aref;
            
    my $data_file_handle;
    my $test_file_handle;
    my $name_file_handle;
    #if we didn't specify a data file. probably the normal case.
    unless($self->data_file) {
        $self->data_file("/tmp/C45.data");
        $data_file_handle=IO::File->new(">" . $self->data_file);
        ###use some procedure to fill this with a random subset of shit(stuff).
    } 
    unless($self->test_file) {
        $self->test_file("/tmp/C45.test");
        $test_file_handle=IO::File->new(">" . $self->test_file);
        #concatenate both files and shove them in here.
    } 
    $self->name_file("/tmp/C45.names");
    $name_file_handle=IO::File->new(">" . $self->name_file);
    ###use some procedure to fill this with a random subset of shit(stuff)


#write the names file
    print $name_file_handle $names;


#Use splice to hack apart the array into a test set and a training set.
#Then join together with commas etc
#The following line does both together
    my @data_lines;
    @data_lines = map {join ", ", @$_} splice(@data,0,$self->gold_set_size),@$db_germline_aref,@$db_wildtype_aref,@$gold_wildtype_aref;

    print $data_file_handle join "\n", @data_lines;
    print $data_file_handle "\n";

    #put the rest of the gold set into the test set
    my @test_lines = map {join ", ", @$_} @data;
    print $test_file_handle join "\n", @test_lines;
    print $test_file_handle "\n";

    ###Run C4.5 here



}
1;

sub help_detail {
    "This module is intended to be a front end for running C4.5 to generate decision trees"
}

#Create hashes of gold SNPs

sub create_gold_snp_hashes {
    my ($self,$gold_fh) = @_;

    #create db instance to check the ref seq
    my $refdb = Bio::DB::Fasta->new($self->REFDIR);

    #create temporary counter variables
    
    my %heterozygous_snp_at;
    my %reference_snp_at;
    my %homozygous_snp_at;
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
                    }
                    else {
                        #Gold SNP reference base is altered
                        #So it is not actually a reference base
                        $self->error_message("Gold SNP reference doesn't match B36 reference sequence");
                    }
                }
                else {
                    if($allele1 ne $ref_a) {
                        $homozygous_snp_at{$chr}{$pos} = [$allele1,$allele2];
                    }
                    else {
                        $self->error_message("Gold SNP is listed as reference in B36 sequence");
                    }
                }

            }
            else {
                #platforms disagree
                if($allele1_type1 ne $allele2_type1 || $allele1_type2 ne $allele2_type2) {
                    $self->error_message("Inconsistent types within a platform on a homozygous SNP at ".$gold_fh->input_line_number);
                    next;
                }
                #Check if the allele matches the reference base in B36
                if($allele1 ne $ref_a) {
                    #It's a SNP!
                    $homozygous_snp_at{$chr}{$pos} = [$allele1,$allele2];
                }
                else {
                    #it was in our reference as ref
                    $reference_snp_at{$chr}{$pos} = [$allele1,$allele2];
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
                    $heterozygous_snp_at{$chr}{$pos} = [$allele1,$allele2];
                }
                elsif($allele2_type1 eq 'SNP' ) {
                    $heterozygous_snp_at{$chr}{$pos} = [$allele2,$allele1];
                }
                else {
                    #something is up, neither allele is labeled as SNP
                    $self->error_message("Supposedly heterozygous SNP not labeled as such at " . $gold_fh->input_line_number);
                }
                    
            }
            else {
                $self->error_message("Platforms disagree on reference at line ".$gold_fh->input_line_number);
            }
        }
    }
    return (\%heterozygous_snp_at,\%homozygous_snp_at,\%reference_snp_at);
}

sub make_validation_hash {
    my $self = shift;
    my @variants = Genome::VariantReviewDetail->get();

    my %validated_status_for;
    
    #Use the following ranking system to resolve discrepancies
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
        my $variant_type = $variant->variant_type;
        next if(!defined($chr) || !defined($pos) || !defined($variant_type) || ($variant_type ne 'S'));

        if(!defined($status)) {
            #assume it failed manual review
            $status = 'FAILED MANUAL REVIEW';
            if(!exists($validated_status_for{$chr}{$pos})) {
                $validated_status_for{$chr}{$pos} = $status;
            }
            next;
        }
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
                    if(exists($rank{$discr_status}) && exists($rank{$new_status}) && $rank{$discr_status} < $rank{$new_status} ) {
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
                        $self->error_message("Unexpected status |$discr_status|");
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
            #This seems sort of dumb, why not use the ranking system
            #check and see if it is an improvement over NC or LQ
            if(uc($status) eq 'G' || uc($status) eq 'WT' || uc($status) eq 'S') {
                $validated_status_for{$chr}{$pos} = $status;
            }
        }
    }
    return \%validated_status_for;
}


sub make_data_arrays {
    my ($self, $handle, $dtr, $gold_het_href, $gold_hom_href, $gold_ref_href, $db_status_href) = @_;
    #useful array function
    #@shuffled = shuffle(@list);

    my @gold_germline;
    my @gold_wildtype;
    my @db_germline;
    my @db_wildtype;

    my %seen_at;

    while(my $line = $handle->getline) {
        #skip if header line
        next if $line =~ /^chromosome/;
        chomp $line;
        my ($chr,
            $pos,
            $al1,
            $al2,
        ) = split ",", $line;
        #create dependent variable ratios
        #everything is dependent on # of reads so make a ratio
        my @attributes = $dtr->make_attribute_array($line);

        if($self->exclude_y) {
            next if($chr eq 'Y'); #female patient these are BS
        }
        next if($al1 eq ' ' || $al1 eq '' || $al1 eq 'N'); #skip 'SNPs' where the reference is N or non-existent
        next if $seen_at{$chr}{$pos}{$al2};
        next if !@attributes;
        $seen_at{$chr}{$pos}{$al2} = 1;
        
        #First check the db
        if($self->use_database && exists($db_status_href->{$chr}{$pos})) {
            #Then there is a status in the db
            if($db_status_href->{$chr}{$pos} eq 'G' ||
                $db_status_href->{$chr}{$pos} eq 'S') {
                #then it is a true SNP
                push @attributes, 'G';
                push @db_germline, \@attributes;
                next;
            }
            if($db_status_href->{$chr}{$pos} eq 'WT') {
                push @attributes, 'WT';
                push @db_wildtype, \@attributes;
                next;
            }
        }

        if(exists($gold_het_href->{$chr}{$pos})) {
            if( (@{$gold_het_href->{$chr}{$pos}}[0] eq $al1 && 
                    @{$gold_het_href->{$chr}{$pos}}[1] eq $al2)
                || (@{$gold_het_href->{$chr}{$pos}}[0] eq $al2 &&
                    @{$gold_het_href->{$chr}{$pos}}[1] eq $al1)
            ) {
                #add to germline
                push @attributes, 'G';
                push @gold_germline, \@attributes;
                next;
            }

        }

        if(exists($gold_hom_href->{$chr}{$pos}) &&
            @{$gold_hom_href->{$chr}{$pos}}[0] eq $al2) {
            push @attributes, 'G';
            push @gold_germline, \@attributes;
            next;
        }

        if(exists($gold_ref_href->{$chr}{$pos})) { 
            push @attributes, 'WT';
            push @gold_wildtype, \@attributes;
            next;
        }
    }
    return (\@db_wildtype,\@db_germline,\@gold_wildtype,\@gold_germline);
}


####PSEUDOCODE#####
#INPUT: LIST OF GOOD SNPS with experimental appended
#INPUT: LIST OF BAD SNPS with experimental appended



####FIRST THING TO DO
####MODIFY COLUMNS TO MAKE ANY RELATED COLUMNS INDEPENDENT
####APPEND STATUS OF ,G to GOOD SNPS file. 
### APPEND STATUS OF ,WT TO BAD SNPS file.

##LIST ALTERATION/PROCESSING DONE

###CONCATENATE LISTS INTO NEW FILE

##GENERATE ANCILLARY CONFIG FILES FOR C4.5
##.NAMES - names of columns and how they vary. "continous"  "discrete: good, bad"
##.DATA - training set
##.TEST - test set
