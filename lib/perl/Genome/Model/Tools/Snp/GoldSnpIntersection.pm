package Genome::Model::Tools::Snp::GoldSnpIntersection;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;
use Command;
use IO::File;
use Bio::DB::Fasta;
use XML::LibXML;

class Genome::Model::Tools::Snp::GoldSnpIntersection {
    is => 'Command',
    has => [
    snp_file => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "maq cns2snp output or samtools pileup SNP output",
    },
    gold_snp_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "input file of snp locations and calls from the intersection of the affy and illumina platforms",
    },
    missed_snp_file => {
        type => 'String',
        is_optional => 1,
        doc => "file name of missed gold snv list",
    },
    exclude_y =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Don't consider SNPs present on the Y chromosome",
        default => 0,
    },        
    total_gold_heterozygote_snps => {
        type => 'Integer',
        is_optional => 1,
        doc => "Instance variable",
        default => 0,
    },
    total_gold_homozygous_ref_positions => {
        type => 'Integer',
        is_optional => 1,
        doc => "Instance variable",
        default => 0,
    },
    total_gold_homozygous_snps => {
        type => 'Integer',
        is_optional => 1,
        doc => "Instance variable",
        default => 0,
    },
    REFDIR =>
    {
        #This is the path to the reference sequence used for aligning the model
        type => 'String',
        is_optional => 0,
        default => "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c",
    },        
    snp_format => {
        type => 'String',
        is_optional => 1,
        doc => 'It should be either maq or sam. default is maq.',
        default => 'MAQ',
    },
    report_format => {
        type => 'String',
        is_optional => 1,
        doc => 'Format of the report as sent to STDOUT. txt or xml. Defaults to txt.',
        default => 'txt'
    },
    _report_txt => {
        type => 'String',
        is_optional => 1,
        doc => 'Report results stored here as text after execute.'
    },
    _report_xml => {
        type => 'String',
        is_optional => 1,
        doc => 'Report results stored here as xml after execute.'
    },                
    report_file => {
        is_optional => 1,
        doc => 'A file path to write the report to.',
    },
    _report_fh => {
        is_optional => 1,
        doc => 'A filehanlde to write the report to.',
    },
    ]
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    if ($self->report_file) {
        $self->_report_fh(Genome::Sys->open_file_for_writing($self->report_file));
    } else {
        open( my $fh,'>&STDOUT');
        unless ($fh) {
            die('Failed to open a filehandle to redirect output to STDOUT');
        }
        $self->_report_fh($fh);
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $snp_format = uc $self->snp_format;
    unless ($snp_format =~ /^(MAQ|SAM)$/) {
        $self->error_message("Invalid snp format: $snp_format");
        return;
    }
    
    #Check on the file names
    unless (-f $self->snp_file) {
        $self->error_message("Snps file is not a file: " . $self->snp_file);
        return;
    }
    unless (-f $self->gold_snp_file) {
        $self->error_message("Gold snp file is not a file: " . $self->gold_snp_file);
        return;
    }

    #Check and open filehandles
    my $snp_fh = IO::File->new($self->snp_file);
    unless ($snp_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->snp_file );
        return;
    }
    my $gold_fh = IO::File->new($self->gold_snp_file);
    unless ($gold_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->gold_snp_file );
        return;
    }

    my ($gold_het_hash_ref, $gold_hom_hash_ref, $gold_ref_hash_ref) = $self->create_gold_snp_hashes($gold_fh);
        
    unless (defined($gold_het_hash_ref)) {
        $self->error_message("Fatal error creating Gold SNP hash");
        return;
    }

    #Grab metrics
    my ($total_snp_positions,$ref_breakdown_ref, $het_breakdown_ref, $hom_breakdown_ref, $missed_snps) 
        = $self->calculate_metrics($snp_fh,$gold_het_hash_ref,$gold_hom_hash_ref,$gold_ref_hash_ref);
    
    $self->store_report_txt($ref_breakdown_ref,$self->total_gold_homozygous_ref_positions,$het_breakdown_ref, $self->total_gold_heterozygote_snps, $hom_breakdown_ref, $self->total_gold_homozygous_snps);

    $self->store_report_xml($ref_breakdown_ref,$self->total_gold_homozygous_ref_positions,$het_breakdown_ref, $self->total_gold_heterozygote_snps, $hom_breakdown_ref, $self->total_gold_homozygous_snps);

    if ($self->missed_snp_file) {
        my $miss_fh = Genome::Sys->open_file_for_writing($self->missed_snp_file);
        unless ($miss_fh) {
            $self->error_message('Failed to open missed_snp_file for writing: '.$self->missed_snp_file);
            return;
        }
    
        for my $type (keys %$missed_snps) {
            $miss_fh->print($type.":\n");
            for my $chr (sort {$a cmp $b} keys %{$missed_snps->{$type}}) {
                for my $pos (sort{$a <=> $b} keys %{$missed_snps->{$type}->{$chr}}) {
                    $miss_fh->print("$chr\t$pos\t".$missed_snps->{$type}->{$chr}->{$pos}."\n");
                }
            }
        }
    }
                
    # print report in requested format
    my $_report_fh = $self->_report_fh;
    if ($self->report_format eq "xml") {
        print $_report_fh $self->_report_xml;
    } elsif ($self->report_format eq "txt") {
        print $_report_fh $self->_report_txt;
    } else {
        $self->error_message("Report format '" . $self->report_format . "' not supported, please specify 'txt' or 'xml'.");
    }
    $_report_fh->close;
    $gold_fh->close;
    $snp_fh->close;

    return 1;
}

sub help_brief {
    "Performs a by-genotype comparison of a cns2snp file and a Gold SNP File";
}

sub help_detail {
    "This script performs a comparison of a maq cns2snp or samtools pileup SNP output file with a Gold SNP file. The comparisons are made on a by-genotype basis. Matches are reported only on an exact genotype match. Each type of array call is reported with the maq/sam calls broken down by match vs. partial match vs. mismatch and then further by type. In addition, the number of each is reported along with percentage of total calls and the average depth for those calls. Currently, no distinction is made between heterozygous Gold calls where one of the alleles is the reference and heterozygous Gold calls where neither allele is the reference. These are unlikely to occur, but this should be improved upon on some point. For maq/sam calls, the following types are reported:
'homozygous reference' - two alleles are reported, both are identical to the reference allele
'homozygous variant' - two alleles are reported, both are identical, but not the reference
'heterozygous - 1 allele variant' - two different alleles are reported, one is reference, the other is variant
'heterozygous - 2 alleles variant' - two different alleles are reported, neither are reference
'tri-allelic - with reference' - three different alleles are reported, one is reference
'tri-allelic - no reference' - three different allelic are reported, none are reference
'ambiguous call' - variant call met none of the above criteria. Should be N"
}


sub calculate_metrics {
    my ($self,$snp_fh,$gold_het_hash_ref_in, $gold_hom_hash_ref_in, $gold_ref_hash_ref_in ) = @_;

    my $gold_het_hash_ref = { %$gold_het_hash_ref_in };
    my $gold_hom_hash_ref = { %$gold_hom_hash_ref_in };
    my $gold_ref_hash_ref = { %$gold_ref_hash_ref_in };

    my %het_breakdown;
    my %hom_breakdown;
    my %ref_breakdown;
    my $total_snp_positions = 0;
    my $snp_format = uc $self->snp_format;

    my $exclude_y = $self->exclude_y;
    
    #no header in cns2snp
    while(my $line = $snp_fh->getline) {
        chomp $line;
        my ($chr,$pos,$ref,$call,$quality,@metrics) = split /\t/, $line; 
        my $rd_depth = $snp_format eq 'SAM' ? $metrics[2] : $metrics[0];
        
        if($exclude_y) {
            next if($chr eq 'Y'); #female patient these are BS
        }

        next if($ref eq ' ' || $ref eq '' || $ref eq 'N'); #skip 'SNPs' where the reference is N or non-existent
        $total_snp_positions++;
        
        my $snp_type = $self->define_snp_call($ref,$call); #get string describing snp type
        
        if(exists($gold_het_hash_ref->{$chr}{$pos})) {
            #Gold standard het call
            my $comparison = $self->compare_gold_to_call($gold_het_hash_ref->{$chr}{$pos},$call);
            $het_breakdown{$comparison}{$snp_type}{n} += 1;
            #print STDERR $line, "\n" if($self->compare_gold_to_maq($gold_het_hash_ref->{$chr}{$pos},$call) eq 'mismatch' && $maq_type eq 'mono-allelic variant');
            $het_breakdown{$comparison}{$snp_type}{depth} += $rd_depth;
            #if($maq_type eq 'homozygous variant') {
            #    print STDERR qq{"$comparison",},$metrics[0],"\n";
            #}
            delete $gold_het_hash_ref->{$chr}{$pos};
        }
        elsif(exists($gold_ref_hash_ref->{$chr}{$pos})) { 
            #Gold standard ref call
            my $comparison = $self->compare_gold_to_call($gold_ref_hash_ref->{$chr}{$pos},$call);
            $ref_breakdown{$comparison}{$snp_type}{n} += 1;
            $ref_breakdown{$comparison}{$snp_type}{depth} += $rd_depth;
            delete $gold_ref_hash_ref->{$chr}{$pos};
        }
        elsif(exists($gold_hom_hash_ref->{$chr}{$pos})) {
            #gold standard homozygous call at this site
            my $comparison = $self->compare_gold_to_call($gold_hom_hash_ref->{$chr}{$pos},$call);
            $hom_breakdown{$comparison}{$snp_type}{n} += 1;
            $hom_breakdown{$comparison}{$snp_type}{depth} += $rd_depth;
            delete $gold_hom_hash_ref->{$chr}{$pos};
        }
    }
    my $gold_missed = { 'homozygous' => $gold_hom_hash_ref, 'heterozygous' => $gold_het_hash_ref, 'reference' => $gold_ref_hash_ref };
    return ($total_snp_positions, \%ref_breakdown,\%het_breakdown,\%hom_breakdown, $gold_missed);

}

sub print_report {
    my ($self, $ref_breakdown_ref, $gold_ref_total, $het_breakdown_ref, $gold_het_total, $hom_breakdown_ref, $gold_hom_total) = @_;
    my $_report_fh = $self->_report_fh;
    print $_report_fh "There were $gold_ref_total homozygous ref sites\n";
    $self->print_breakdown($gold_ref_total,$ref_breakdown_ref);
    print $_report_fh "There were $gold_het_total heterozygous calls (could include bi-allelic calls)\n";
    $self->print_breakdown($gold_het_total,$het_breakdown_ref);
    print $_report_fh "There were $gold_hom_total homozygous calls\n";
    $self->print_breakdown($gold_hom_total, $hom_breakdown_ref);
}

sub print_breakdown {
    my ($self, $total, $hash) = @_;
    #first print predominant class (match)
    my $_report_fh = $self->_report_fh;
    if(exists($hash->{'match'})) {
        print $_report_fh "\tMatching Gold Genotype\n";
        foreach my $type (keys %{$hash->{'match'}}) {
            printf $_report_fh (
                           "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                           $type,
                           $hash->{'match'}{$type}{'n'},
                           $hash->{'match'}{$type}{'n'}/$total*100,
                           $hash->{'match'}{$type}{'depth'}/$hash->{'match'}{$type}{'n'}
                          );
        }
    }
    if(exists($hash->{'partial match'})) {
        print $_report_fh "\tPartially Matching Gold Genotype\n";
        foreach my $type (keys %{$hash->{'partial match'}}) {
            printf $_report_fh (
                           "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                           $type,
                           $hash->{'partial match'}{$type}{'n'},
                           $hash->{'partial match'}{$type}{'n'}/$total*100,
                           $hash->{'partial match'}{$type}{'depth'}/$hash->{'partial match'}{$type}{'n'}
                          );
        }
    }
    #next print un-matching classes
    if(exists($hash->{'mismatch'})) {
        print $_report_fh "\tMismatching Gold Genotype\n";
        foreach my $type (keys %{$hash->{'mismatch'}}) {
            printf $_report_fh (
                           "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                           $type,
                           $hash->{'mismatch'}{$type}{'n'},
                           $hash->{'mismatch'}{$type}{'n'}/$total*100,
                           $hash->{'mismatch'}{$type}{'depth'}/$hash->{'mismatch'}{$type}{'n'}
                          );
        }
    }
}

sub store_report_txt {
    my ($self, $ref_breakdown_ref, $gold_ref_total, $het_breakdown_ref, $gold_het_total, $hom_breakdown_ref, $gold_hom_total) = @_; 
    my $report_txt;
    # store txt report
    $report_txt .= "There were $gold_ref_total homozygous ref sites\n";
    $report_txt .= $self->get_breakdown($gold_ref_total,$ref_breakdown_ref);
    $report_txt .= "There were $gold_het_total heterozygous calls (could include bi-allelic calls)\n";
    $report_txt .= $self->get_breakdown($gold_het_total,$het_breakdown_ref);
    $report_txt .=  "There were $gold_hom_total homozygous calls\n";
    $report_txt .= $self->get_breakdown($gold_hom_total, $hom_breakdown_ref);

    $self->_report_txt($report_txt);

}

sub get_breakdown {
    my ($self, $total, $hash) = @_;
    my $breakdown;
    #first print predominant class (match)
    if (exists($hash->{'match'})) {
        $breakdown .= "\tMatching Gold Genotype\n";
        foreach my $type (keys %{$hash->{'match'}}) {
            $breakdown .= sprintf(
                                  "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                                  $type,
                                  $hash->{'match'}{$type}{'n'},
                                  $hash->{'match'}{$type}{'n'}/$total*100,
                                  $hash->{'match'}{$type}{'depth'}/$hash->{'match'}{$type}{'n'}
                                 );
        }
    }
    if (exists($hash->{'partial match'})) {
        $breakdown .= "\tPartially Matching Gold Genotype\n";
        foreach my $type (keys %{$hash->{'partial match'}}) {
            $breakdown .= sprintf (
                                   "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                                   $type,
                                   $hash->{'partial match'}{$type}{'n'},
                                   $hash->{'partial match'}{$type}{'n'}/$total*100,
                                   $hash->{'partial match'}{$type}{'depth'}/$hash->{'partial match'}{$type}{'n'}
                                  );
        }
    }
    #next print un-matching classes
    if (exists($hash->{'mismatch'})) {
        $breakdown .= "\tMismatching Gold Genotype\n";
        foreach my $type (keys %{$hash->{'mismatch'}}) {
            $breakdown .= sprintf(
                                  "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                                  $type,
                                  $hash->{'mismatch'}{$type}{'n'},
                                  $hash->{'mismatch'}{$type}{'n'}/$total*100,
                                  $hash->{'mismatch'}{$type}{'depth'}/$hash->{'mismatch'}{$type}{'n'});
        }
    }
    return $breakdown;
}



#Functions to create hashes of data

#Create hashes of gold SNPs

sub create_gold_snp_hashes {
    my ($self,$gold_fh) = @_;

    #create db instance to check the ref seq
    my $refdb = Bio::DB::Fasta->new($self->REFDIR);

    #create temporary counter variables
    my $total_gold_homozygous_ref_positions;
    my $total_gold_homozygous_snp_positions;
    my $total_gold_heterozygote_snps;
    
    my %heterozygous_snp_at;
    my %homozygous_snp_at;
    my %reference_snp_at;

    #there is no header on this file
    #Format is tab separated
    #Chr\tPos\tPos\tAllele1\tAllele2\tPlatform1_Allele1\tPlatform1_Allele2\tPlatform2_Allele1\tPlatform2_Allele2
    while(my $line = $gold_fh->getline) {
        chomp $line;
        my ($chr, $pos, $pos2, $allele1, $allele2, $allele1_type1,$allele2_type1, $allele1_type2, $allele2_type2) = split /\t/, $line;

        my $ref_a= $refdb->seq($chr, $pos => $pos2);  
        unless($ref_a) {
            $self->error_message("No reference base for position $chr $pos $pos2");
            next;
        }
        chomp($ref_a);
        $ref_a=uc($ref_a);

        if($allele1 eq $allele2) {
            #homozygous call
            if($allele1_type1 eq $allele1_type2) {
                #Check that the platform agree is internally consistent
                if($allele1_type1 ne $allele2_type1 || $allele1_type2 ne $allele2_type2) {
                    $self->error_message("Inconsistent types within a platform on a homozygous SNP at " . $gold_fh->input_line_number);
                    next;
                }
                if($allele1_type1 eq 'ref') {
                    if($allele1 eq $ref_a) {
                        #it was in our reference as ref
                        $reference_snp_at{$chr}{$pos} = $allele1.$allele2;
                        $total_gold_homozygous_ref_positions++;
                    }
                    else {
                        #Gold SNP reference base is altered
                        #So it is not actually a reference base
                        $self->error_message("Gold SNP reference doesn't match B36 reference sequence");
                    }
                }
                else {
                    #assuming homozygous SNP at this position
                    if($allele1 ne $ref_a) {
                        $homozygous_snp_at{$chr}{$pos} = $allele1.$allele2;
                        $total_gold_homozygous_snp_positions++;
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
                    $homozygous_snp_at{$chr}{$pos} = $allele1.$allele2;
                    $total_gold_homozygous_snp_positions++;
                }
                else {
                    #it was in our reference as ref
                    $reference_snp_at{$chr}{$pos} = $allele1.$allele2;
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
                if($allele1_type1 eq $allele2_type1) {
                    #non-ref bi-allelic SNP unlikely and unhandled. Let the user know
                    $self->error_message("Heterozygous snp where both alleles are non-reference detected at $chr\t$pos. Not added to hash");
                    #ignore
                }
                elsif($allele1_type1 eq 'SNP' ) {
                    #$total_gold++;
                    $total_gold_heterozygote_snps++;
                    $heterozygous_snp_at{$chr}{$pos} = $allele1.$allele2;
                }
                elsif($allele2_type1 eq 'SNP' ) {
                    #$total_gold++;
                    $total_gold_heterozygote_snps++;
                    $heterozygous_snp_at{$chr}{$pos} = $allele2.$allele1;
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
    $self->total_gold_homozygous_snps($total_gold_homozygous_snp_positions);
    

    return (\%heterozygous_snp_at,\%homozygous_snp_at,\%reference_snp_at);
}

sub define_snp_call {
    my ($self, $ref, $call) = @_;
    my $iub_string = Genome::Info::IUB->iub_to_string($call);
    
    if($self->is_homozygous_IUB($call)) {
        #homozygous call
        if($ref eq $call) {
            #homozygous ref call, will not happen
            return 'homozygous reference';
        }
        else {
            return 'homozygous variant';
        }
    }
    elsif(length $iub_string == 2) {
        #het call
        if($iub_string =~ qr{$ref}) {
            return 'heterozygous - 1 allele variant';
        }
        else {
            return 'heterozygous - 2 alleles variant';
        }
    }
    elsif(length $iub_string == 3) {
        #tri-allelic
        if($iub_string =~ qr{$ref}) {
            return 'tri-allelic - with reference';
        }
        else {
            return 'tri-allelic - no reference';
        }
    }
    else {
        #N
        return 'ambiguous call';
    }
}

sub is_homozygous_IUB {
    my ($self, $call) = @_;
    if($call =~ /[ACGT]/) {
        return 1;
    }
    else {
        return 0;
    }
}

sub compare_gold_to_call {
    my ($self, $gold_alleles, $call) = @_;
    my $alleles = Genome::Info::IUB->iub_to_string($call);
    
    if($gold_alleles eq $alleles || scalar(reverse($gold_alleles)) eq $alleles) {
        return 'match';
    }
    elsif( $self->overlap($gold_alleles, $call)) {
        return 'partial match';
    }
    else {
        return 'mismatch';
    }
}
    
#idea borrowed from ssmith's gmt snp intersect
sub overlap {

    my ($self, $gold_alleles, $call) = @_;
    my $iub_string = Genome::Info::IUB->iub_to_string($call);
    
    my $num_bases_overlapping = 0;
    my %gold_bases = map {$_ => 1} (split //, $gold_alleles); #make a list of unique bases in the gold call
    
    for my $base (keys %gold_bases) {
       if (index($iub_string, $base) != -1) {  
           $num_bases_overlapping++;
       }
    }
    
    return $num_bases_overlapping;
}

#
# XML reporting subs
#

sub store_report_xml {
    my ($self, $ref_breakdown_ref, $gold_ref_total, $het_breakdown_ref, $gold_het_total, $hom_breakdown_ref, $gold_hom_total) = @_; 

    my $report_xml;
    $self->{_xml} = XML::LibXML->createDocument;
    unless ( $self->{_xml} ) {
        $self->error_message("Can't create XML object");
        return;
    }

    # main node
    $self->{_main_node} = $self->_xml->createElement('gold-snp-intersection-report')
      or Carp::confess('Create main node to XML');
    $self->_xml->addChild( $self->_main_node )
      or Carp::confess('Add main node to XML');
    $self->_xml->setDocumentElement( $self->_main_node );
    
    # assemble gold-homozygous-ref node
    $self->{_gold_hom_ref_node} =  $self->_xml->createElement("gold-homozygous-ref");
    $self->_gold_hom_ref_node->addChild( $self->_xml->createAttribute("total", $gold_ref_total) );
    $self->_add_variants_node($gold_ref_total, $ref_breakdown_ref, $self->_gold_hom_ref_node);

    # assemble gold-heterozygous-snp node
    $self->{_gold_het_snp_node} =  $self->_xml->createElement("gold-heterozygous-snp");
    $self->_gold_het_snp_node->addChild( $self->_xml->createAttribute("total", $gold_het_total) );
    $self->_add_variants_node($gold_het_total, $het_breakdown_ref, $self->_gold_het_snp_node);

    # assemble gold-homozygous-snp node   
    $self->{_gold_hom_snp_node} =  $self->_xml->createElement("gold-homozygous-snp");
    $self->_gold_hom_snp_node->addChild( $self->_xml->createAttribute("total", $gold_hom_total) );
    $self->_add_variants_node($gold_hom_total, $hom_breakdown_ref, $self->_gold_hom_snp_node);

    # add nodes
    $self->_main_node->addChild( $self->_gold_hom_ref_node );
    $self->_main_node->addChild( $self->_gold_het_snp_node );
    $self->_main_node->addChild( $self->_gold_hom_snp_node );
    
    # store xml report
    $self->_report_xml($self->_xml->toString(1));

}

sub _add_variants_node {
    my ($self, $total, $hash, $node) = @_;
    
    #first print predominant class (match)
    if (exists($hash->{'match'})) {
        my $match_node = $node->addChild($self->_xml->createElement("match"));
        foreach my $type (keys %{$hash->{'match'}}) {
            my $variant_node = $match_node->addChild($self->_xml->createElement("variant"));
            $variant_node->appendTextChild("type", $type);
            $variant_node->appendTextChild("reads", $hash->{'match'}{$type}{'n'});
            $variant_node->appendTextChild("intersection", $hash->{'match'}{$type}{'n'}/$total*100);
            $variant_node->appendTextChild("depth", $hash->{'match'}{$type}{'depth'}/$hash->{'match'}{$type}{'n'});
        }
    }

    if (exists($hash->{'partial match'})) {
        my $partial_match_node = $node->addChild($self->_xml->createElement("partial-match"));
        foreach my $type (keys %{$hash->{'partial match'}}) {
            my $variant_node = $partial_match_node->addChild($self->_xml->createElement("variant"));
            $variant_node->appendTextChild("type", $type);
            $variant_node->appendTextChild("reads", $hash->{'partial match'}{$type}{'n'});
            $variant_node->appendTextChild("intersection", $hash->{'partial match'}{$type}{'n'}/$total*100);
            $variant_node->appendTextChild("depth", $hash->{'partial match'}{$type}{'depth'}/$hash->{'partial match'}{$type}{'n'});
        }
    }

    #next print un-matching classes
    if (exists($hash->{'mismatch'})) {
        my $mismatch_node = $node->addChild($self->_xml->createElement("mismatch"));
        foreach my $type (keys %{$hash->{'mismatch'}}) {
            my $variant_node = $mismatch_node->addChild($self->_xml->createElement("variant"));
            $variant_node->appendTextChild("type", $type);
            $variant_node->appendTextChild("reads", $hash->{'mismatch'}{$type}{'n'});
            $variant_node->appendTextChild("intersection", $hash->{'mismatch'}{$type}{'n'}/$total*100);
            $variant_node->appendTextChild("depth", $hash->{'mismatch'}{$type}{'depth'}/$hash->{'mismatch'}{$type}{'n'});
        }
    }
}

sub _xml {
    return $_[0]->{_xml};
}

sub _main_node {
    return $_[0]->{_main_node};
}

sub _gold_hom_ref_node {
    return $_[0]->{_gold_hom_ref_node};
}

sub _gold_het_snp_node {
    return $_[0]->{_gold_het_snp_node};    
}

sub _gold_hom_snp_node {
    return $_[0]->{_gold_hom_snp_node};
}

1;
