package Genome::Model::Tools::Snp::ExpectedSnvIntersection;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;
use Command;
use IO::File;
use Bio::DB::Fasta;
use XML::LibXML;

class Genome::Model::Tools::Snp::ExpectedSnvIntersection {
    is => 'Command',
    has => [
    test_snv_file => 
    { 
        type => 'String',
        is_optional => 0,
        doc => "maq cns2snv output or samtools pileup SNV output",
    },
    expected_snv_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "input file of snv locations and calls from the intersection of the affy and illumina platforms",
    },
    exclude_y =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Don't consider SNVs present on the Y chromosome",
        default => 0,
    },        
    total_expected_heterozygote_snvs => {
        type => 'Integer',
        is_optional => 1,
        doc => "Instance variable",
        default => 0,
    },
    total_expected_homozygous_ref_positions => {
        type => 'Integer',
        is_optional => 1,
        doc => "Instance variable",
        default => 0,
    },
    total_expected_homozygous_snvs => {
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
    snv_format => {
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
    }                
            
    ]
};


sub execute {
    my $self = shift;

    my $snv_format = uc $self->snv_format;
    unless ($snv_format =~ /^(MAQ|SAM)$/) {
        $self->error_message("Invalid snv format: $snv_format");
        return;
    }
    
    #Check on the file names
    unless (-f $self->snv_file) {
        $self->error_message("Snps file is not a file: " . $self->snv_file);
        return;
    }
    unless (-f $self->expected_snv_file) {
        $self->error_message("expected snv file is not a file: " . $self->expected_snv_file);
        return;
    }

    #Check and open filehandles
    my $snv_fh = IO::File->new($self->snv_file);
    unless ($snv_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->snv_file );
        return;
    }
    my $expected_fh = IO::File->new($self->expected_snv_file);
    unless ($expected_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->expected_snv_file );
        return;
    }

    my ($expected_het_hash_ref, $expected_hom_hash_ref, $expected_ref_hash_ref) = $self->create_expected_snv_hashes($expected_fh);
        
    unless (defined($expected_het_hash_ref)) {
        $self->error_message("Fatal error creating expected SNV hash");
        return;
    }

    #Grab metrics
    my ($total_snv_positions,$ref_breakdown_ref, $het_breakdown_ref, $hom_breakdown_ref) 
        = $self->calculate_metrics($snv_fh,$expected_het_hash_ref,$expected_hom_hash_ref,$expected_ref_hash_ref);
    
    $self->store_report_txt($ref_breakdown_ref,$self->total_expected_homozygous_ref_positions,$het_breakdown_ref, $self->total_expected_heterozygote_snvs, $hom_breakdown_ref, $self->total_expected_homozygous_snvs);

    $self->store_report_xml($ref_breakdown_ref,$self->total_expected_homozygous_ref_positions,$het_breakdown_ref, $self->total_expected_heterozygote_snvs, $hom_breakdown_ref, $self->total_expected_homozygous_snvs);

    # print report in requested format
    
    if ($self->report_format eq "xml") {
        print STDOUT $self->_report_xml;
    } elsif ($self->report_format eq "txt") {
        print STDOUT $self->_report_txt;
    } else {
        $self->error_message("Report format '" . $self->report_format . "' not supported, please specify 'txt' or 'xml'.");
    }
                              
    $expected_fh->close;
    $snv_fh->close;

    return 1;
}

sub help_brief {
    "Performs a by-genotype comparison of a cns2snv file and a expected SNV File";
}

sub help_detail {
    "This script performs a comparison of a maq cns2snv or samtools pileup SNV output file with a expected SNV file. The comparisons are made on a by-genotype basis. Matches are reported only on an exact genotype match. Each type of array call is reported with the maq/sam calls broken down by match vs. partial match vs. mismatch and then further by type. In addition, the number of each is reported along with percentage of total calls and the average depth for those calls. Currently, no distinction is made between heterozygous expected calls where one of the alleles is the reference and heterozygous expected calls where neither allele is the reference. These are unlikely to occur, but this should be improved upon on some point. For maq/sam calls, the following types are reported:
'homozygous reference' - two alleles are reported, both are identical to the reference allele
'homozygous variant' - two alleles are reported, both are identical, but not the reference
'heterozygous - 1 allele variant' - two different alleles are reported, one is reference, the other is variant
'heterozygous - 2 alleles variant' - two different alleles are reported, neither are reference
'tri-allelic - with reference' - three different alleles are reported, one is reference
'tri-allelic - no reference' - three different allelic are reported, none are reference
'ambiguous call' - variant call met none of the above criteria. Should be N"
}


sub calculate_metrics {
    my ($self,$snv_fh,$expected_het_hash_ref, $expected_hom_hash_ref, $expected_ref_hash_ref, ) = @_;

    my %het_breakdown;
    my %hom_breakdown;
    my %ref_breakdown;
    my $total_snv_positions = 0;
    my $snv_format = uc $self->snv_format;

    my $exclude_y = $self->exclude_y;
    
    #no header in cns2snv
    while(my $line = $snv_fh->getline) {
        chomp $line;
        my ($chr,$pos,$ref,$call,$quality,@metrics) = split /\t/, $line; 
        my $rd_depth = $snv_format eq 'SAM' ? $metrics[2] : $metrics[0];
        
        if($exclude_y) {
            next if($chr eq 'Y'); #female patient these are BS
        }

        next if($ref eq ' ' || $ref eq '' || $ref eq 'N'); #skip 'SNVs' where the reference is N or non-existent
        $total_snv_positions++;
        
        my $snv_type = $self->define_snv_call($ref,$call); #get string describing snv type
        
        if(exists($expected_het_hash_ref->{$chr}{$pos})) {
            #expected standard het call
            my $comparison = $self->compare_expected_to_call($expected_het_hash_ref->{$chr}{$pos},$call);
            $het_breakdown{$comparison}{$snv_type}{n} += 1;
            #print STDERR $line, "\n" if($self->compare_expected_to_maq($expected_het_hash_ref->{$chr}{$pos},$call) eq 'mismatch' && $maq_type eq 'mono-allelic variant');
            $het_breakdown{$comparison}{$snv_type}{depth} += $rd_depth;
            #if($maq_type eq 'homozygous variant') {
            #    print STDERR qq{"$comparison",},$metrics[0],"\n";
            #}
        }
        elsif(exists($expected_ref_hash_ref->{$chr}{$pos})) { 
            #expected standard ref call
            my $comparison = $self->compare_expected_to_call($expected_ref_hash_ref->{$chr}{$pos},$call);
            $ref_breakdown{$comparison}{$snv_type}{n} += 1;
            $ref_breakdown{$comparison}{$snv_type}{depth} += $rd_depth;
        }
        elsif(exists($expected_hom_hash_ref->{$chr}{$pos})) {
            #expected standard homozygous call at this site
            my $comparison = $self->compare_expected_to_call($expected_hom_hash_ref->{$chr}{$pos},$call);
            $hom_breakdown{$comparison}{$snv_type}{n} += 1;
            $hom_breakdown{$comparison}{$snv_type}{depth} += $rd_depth;
        }
    }
    return ($total_snv_positions, \%ref_breakdown,\%het_breakdown,\%hom_breakdown);

}

sub print_report {
    my ($self, $ref_breakdown_ref, $expected_ref_total, $het_breakdown_ref, $expected_het_total, $hom_breakdown_ref, $expected_hom_total) = @_; 
    print STDOUT "There were $expected_ref_total homozygous ref sites\n";
    $self->print_breakdown($expected_ref_total,$ref_breakdown_ref);
    print STDOUT "There were $expected_het_total heterozygous calls (could include bi-allelic calls)\n";
    $self->print_breakdown($expected_het_total,$het_breakdown_ref);
    print STDOUT "There were $expected_hom_total homozygous calls\n";
    $self->print_breakdown($expected_hom_total, $hom_breakdown_ref);
}

sub print_breakdown {
    my ($self, $total, $hash) = @_;
    #first print predominant class (match)
    if(exists($hash->{'match'})) {
        print STDOUT "\tMatching expected Genotype\n";
        foreach my $type (keys %{$hash->{'match'}}) {
            printf STDOUT (
                           "\t\t%s\t%d\t%0.2f\t%0.2f\n",
                           $type,
                           $hash->{'match'}{$type}{'n'},
                           $hash->{'match'}{$type}{'n'}/$total*100,
                           $hash->{'match'}{$type}{'depth'}/$hash->{'match'}{$type}{'n'}
                          );
        }
    }
    if(exists($hash->{'partial match'})) {
        print STDOUT "\tPartially Matching expected Genotype\n";
        foreach my $type (keys %{$hash->{'partial match'}}) {
            printf STDOUT (
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
        print STDOUT "\tMismatching expected Genotype\n";
        foreach my $type (keys %{$hash->{'mismatch'}}) {
            printf STDOUT (
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
    my ($self, $ref_breakdown_ref, $expected_ref_total, $het_breakdown_ref, $expected_het_total, $hom_breakdown_ref, $expected_hom_total) = @_; 
    my $report_txt;
    # store txt report
    $report_txt .= "There were $expected_ref_total homozygous ref sites\n";
    $report_txt .= $self->get_breakdown($expected_ref_total,$ref_breakdown_ref);
    $report_txt .= "There were $expected_het_total heterozygous calls (could include bi-allelic calls)\n";
    $report_txt .= $self->get_breakdown($expected_het_total,$het_breakdown_ref);
    $report_txt .=  "There were $expected_hom_total homozygous calls\n";
    $report_txt .= $self->get_breakdown($expected_hom_total, $hom_breakdown_ref);

    $self->_report_txt($report_txt);

}

sub get_breakdown {
    my ($self, $total, $hash) = @_;
    my $breakdown;
    #first print predominant class (match)
    if (exists($hash->{'match'})) {
        $breakdown .= "\tMatching expected Genotype\n";
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
        $breakdown .= "\tPartially Matching expected Genotype\n";
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
        $breakdown .= "\tMismatching expected Genotype\n";
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

#Create hashes of expected SNVs

sub create_expected_snv_hashes {
    my ($self,$expected_fh) = @_;

    #create db instance to check the ref seq
    my $refdb = Bio::DB::Fasta->new($self->REFDIR);

    #create temporary counter variables
    my $total_expected_homozygous_ref_positions;
    my $total_expected_homozygous_snv_positions;
    my $total_expected_heterozygote_snvs;
    
    my %heterozygous_snv_at;
    my %homozygous_snv_at;
    my %reference_snv_at;

    #there is no header on this file
    #Format is tab separated
    #Chr\tPos\tPos\tAllele1\tAllele2\tPlatform1_Allele1\tPlatform1_Allele2\tPlatform2_Allele1\tPlatform2_Allele2
    while(my $line = $expected_fh->getline) {
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
                    $self->error_message("Inconsistent types within a platform on a homozygous SNV at " . $expected_fh->input_line_number);
                    next;
                }
                if($allele1_type1 eq 'ref') {
                    if($allele1 eq $ref_a) {
                        #it was in our reference as ref
                        $reference_snv_at{$chr}{$pos} = $allele1.$allele2;
                        $total_expected_homozygous_ref_positions++;
                    }
                    else {
                        #expected SNV reference base is altered
                        #So it is not actually a reference base
                        $self->error_message("expected SNV reference doesn't match B36 reference sequence");
                    }
                }
                else {
                    #assuming homozygous SNV at this position
                    if($allele1 ne $ref_a) {
                        $homozygous_snv_at{$chr}{$pos} = $allele1.$allele2;
                        $total_expected_homozygous_snv_positions++;
                    }
                    else {
                        $self->error_message("expected SNV is listed as reference in B36 sequence");
                    }
                }

            }
            else {
                #platforms disagree
                if($allele1_type1 ne $allele2_type1 || $allele1_type2 ne $allele2_type2) {
                    $self->error_message("Inconsistent types within a platform on a homozygous SNV at ".$expected_fh->input_line_number);
                    next;
                }
                #Check if the allele matches the reference base in B36
                if($allele1 ne $ref_a) {
                    #It's a SNV!
                    $homozygous_snv_at{$chr}{$pos} = $allele1.$allele2;
                    $total_expected_homozygous_snv_positions++;
                }
                else {
                    #it was in our reference as ref
                    $reference_snv_at{$chr}{$pos} = $allele1.$allele2;
                    $total_expected_homozygous_ref_positions++;
                }

            }
        }
        else {
            #heterozygous site
            #
            #Check that the platforms agree
            if($allele1_type1 eq $allele1_type2 && $allele2_type1 eq $allele2_type2) {
                #het site
                #check that the allele is actually a snv
                if($allele1_type1 eq $allele2_type1) {
                    #non-ref bi-allelic SNV unlikely and unhandled. Let the user know
                    $self->error_message("Heterozygous snv where both alleles are non-reference detected at $chr\t$pos. Not added to hash");
                    #ignore
                }
                elsif($allele1_type1 eq 'SNV' ) {
                    #$total_expected++;
                    $total_expected_heterozygote_snvs++;
                    $heterozygous_snv_at{$chr}{$pos} = $allele1.$allele2;
                }
                elsif($allele2_type1 eq 'SNV' ) {
                    #$total_expected++;
                    $total_expected_heterozygote_snvs++;
                    $heterozygous_snv_at{$chr}{$pos} = $allele2.$allele1;
                }
                else {
                    #something is up, neither alleles is labeled as SNV
                    $self->error_message("Supposedly heterozygous SNV not labeled as such at " . $expected_fh->input_line_number);
                }
                    
            }
            else {
                $self->error_message("Platforms disagree on reference at line ".$expected_fh->input_line_number);
            }
        }
    }
    #set the class variables
    $self->total_expected_heterozygote_snvs($total_expected_heterozygote_snvs);
    $self->total_expected_homozygous_ref_positions($total_expected_homozygous_ref_positions);
    $self->total_expected_homozygous_snvs($total_expected_homozygous_snv_positions);
    

    return (\%heterozygous_snv_at,\%homozygous_snv_at,\%reference_snv_at);
}

sub define_snv_call {
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

sub compare_expected_to_call {
    my ($self, $expected_alleles, $call) = @_;
    my $alleles = Genome::Info::IUB->iub_to_string($call);
    
    if($expected_alleles eq $alleles || scalar(reverse($expected_alleles)) eq $alleles) {
        return 'match';
    }
    elsif( $self->overlap($expected_alleles, $call)) {
        return 'partial match';
    }
    else {
        return 'mismatch';
    }
}
    
#idea borrowed from ssmith's gmt snv intersect
sub overlap {

    my ($self, $expected_alleles, $call) = @_;
    my $iub_string = Genome::Info::IUB->iub_to_string($call);
    
    my $num_bases_overlapping = 0;
    my %expected_bases = map {$_ => 1} (split //, $expected_alleles); #make a list of unique bases in the expected call
    
    for my $base (keys %expected_bases) {
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
    my ($self, $ref_breakdown_ref, $expected_ref_total, $het_breakdown_ref, $expected_het_total, $hom_breakdown_ref, $expected_hom_total) = @_; 

    my $report_xml;
    $self->{_xml} = XML::LibXML->createDocument;
    unless ( $self->{_xml} ) {
        $self->error_message("Can't create XML object");
        return;
    }

    # main node
    $self->{_main_node} = $self->_xml->createElement('expected-snv-intersection-report')
      or Carp::confess('Create main node to XML');
    $self->_xml->addChild( $self->_main_node )
      or Carp::confess('Add main node to XML');
    $self->_xml->setDocumentElement( $self->_main_node );
    
    # assemble expected-homozygous-ref node
    $self->{_expected_hom_ref_node} =  $self->_xml->createElement("expected-homozygous-ref");
    $self->_expected_hom_ref_node->addChild( $self->_xml->createAttribute("total", $expected_ref_total) );
    $self->_add_variants_node($expected_ref_total, $ref_breakdown_ref, $self->_expected_hom_ref_node);

    # assemble expected-heterozygous-snv node
    $self->{_expected_het_snv_node} =  $self->_xml->createElement("expected-heterozygous-snv");
    $self->_expected_het_snv_node->addChild( $self->_xml->createAttribute("total", $expected_het_total) );
    $self->_add_variants_node($expected_het_total, $het_breakdown_ref, $self->_expected_het_snv_node);

    # assemble expected-homozygous-snv node   
    $self->{_expected_hom_snv_node} =  $self->_xml->createElement("expected-homozygous-snv");
    $self->_expected_hom_snv_node->addChild( $self->_xml->createAttribute("total", $expected_hom_total) );
    $self->_add_variants_node($expected_hom_total, $hom_breakdown_ref, $self->_expected_hom_snv_node);

    # add nodes
    $self->_main_node->addChild( $self->_expected_hom_ref_node );
    $self->_main_node->addChild( $self->_expected_het_snv_node );
    $self->_main_node->addChild( $self->_expected_hom_snv_node );
    
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

sub _expected_hom_ref_node {
    return $_[0]->{_expected_hom_ref_node};
}

sub _expected_het_snv_node {
    return $_[0]->{_expected_het_snv_node};    
}

sub _expected_hom_snv_node {
    return $_[0]->{_expected_hom_snv_node};
}

1;
