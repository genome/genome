#####################################################################################
#   perl module used to parse in a polyscan.out file for polyscan3.0
#*** MG::IO::Polyscan.pm ***#
# Copyright (C) 
#    Ken Chen (kchen@watson.wustl.edu), 2007
# All rights reserved!
#####################################################################################

package MG::IO::Polyscan;

use MG::IO::Polyscan::Contig;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
use warnings;
use Log::Log4perl;
use IO::File;

=head1 NAME

 MG::IO::Polyscan -- read and write information from/to a Polyscan file
 MG::IO::Polyscan::SNP
 MG::IO::Polyscan::INDEL
 MG::IO::Polyscan::Contig

=head1 SYNOPSIS

 my $ps=new MG::IO::Polyscan();
 $ps->readIn($inPolyscanfile);  #read in original polyphred.out

or

 my $pp=new MG::IO::Polyscan(polyscan => $inPolyscanfile);

=head2 Synopsis subheading

 MG::IO::Polyscan::SNP
 MG::IO::Polyscan::INDEL
 MG::IO::Polyscan::Contig

=head1 DESCRIPTION

This module reads and writes from a Polyscan file.

=head1 Constructor and Initialization

=item (object) new (polyscan => Polyscan_filename)

=head2 arguments (optional)

=item yamlfile

 (Default: $0.yaml)     The YAML file takes the same arguments as given here
                        (except the 'yamlfile' argument).  See the synopsis
                        for examples.

=item polyscan

                        Polyscan (3.0) filename

=head2 Methods

=item examethod(argument)

Description of example method.

=cut

#####################################################################################
# MG::IO::Polyscan
#####################################################################################

sub new {
 	my ($class, %arg) = @_;
=cut
	my ($config_object) = new MG::Config(yamlfile => $arg{yamlfile});
	my ($config) = $config_object->Config();
=cut
	my $self =
	  {
	   _logger          => $arg{logger},
	   _polyscan    	=> $arg{polyscan},
       fasta            => $arg{fasta},
       begin_position   => $arg{begin_position},
       end_position     => $arg{end_position},
       chromosome       => $arg{chromosome},
       #|| $config->{polyscan}
       #|| undef,
	   _param=>undef,
	   _ncontigs=>undef,
	   _contigs=>undef
	  };
 	bless($self, $class || ref($class));
    $self->{_assembly_project_name} = $arg{assembly_project_name} if exists $arg{assembly_project_name};

	if (defined($self->{_polyscan})) {
	  $self->readIn($self->{_polyscan}); # read in original polyscan.out
	}
 	return $self;
}

#-------------------------------------------------
sub setParam{
  my ($self, $r_param)=@_;
  $self->{_param}=$r_param if defined($r_param);
}

sub getParam {
  my ($self)=@_;
  return $self->{_param};
}

#-------------------------------------------------
sub setNContigs{
  my ($self, $nc)=@_;
  $self->{_ncontigs}=$nc if defined($nc);
}

#-------------------------------------------------
sub getNContigs{
  my ($self)=@_;
  return $self->{_ncontigs};
}

#-------------------------------------------------
sub setContigs{
  my ($self, $r_contigs)=@_;
  $self->{_contigs}=$r_contigs if defined($r_contigs);
}

#-------------------------------------------------
sub getContig{
  my ($self, $id_contig)=@_;
  my @contigs=@{$self->{_contigs}};
  return $contigs[$id_contig];
}

#-------------------------------------------------
sub readIn{
  my ($self, $infile)=@_;
  my($log) = $self->{_logger};
  unless (open(IN,"<$infile")) {
      #$log->error("unable to open $infile");
    return;
  }

  do {
    $_=<IN>; chomp;
  }until($_=~/BEGIN\_PARAMETERS/ || eof(IN));
  my @para;

  do {
    $_=<IN>; chomp;
    push @para, $_;
  }until($_=~/END\_PARAMETERS/ || eof(IN));

  pop @para;
  $self->{_param}=\@para;

  return if(eof(IN));  #return if polyscan file is empty

  do {
    $_=<IN>; chomp;
  }until($_=~/BEGIN\_CONTIG/ || eof(IN));

  my @contigs;
  do {
    my $contig=new MG::IO::Polyscan::Contig();
    $contig->ReadIn(\*IN);
    push @contigs,$contig;

    do {
      $_=<IN>;
    }until( eof(IN) || $_=~/BEGIN\_CONTIG/ );

  } until (eof(IN));

  $self->setContigs(\@contigs);
  $self->setNContigs($#contigs+1);
}


#-------------------------------------------------
sub getSNPTable{
  my ($self)=@_;
  my %tab_SNP;
  my $ncontigs=$self->{_ncontigs};
  if($ncontigs>0){
    my @contigs=@{$self->{_contigs}};
    foreach my $r_c (@contigs) {
      my $ctab_SNP=$r_c->getSNP_READSITES();
      foreach my $cpos(keys %{$ctab_SNP}){
	foreach my $rd(keys %{$$ctab_SNP{$cpos}}){
	  my ($refpos,$rdpos, $read, $allel1, $allel2, $score, $var_tag)=$$ctab_SNP{$cpos}{$rd}->get();
	  my @gt=($allel1,$allel2,$score);
	  $tab_SNP{$cpos}{$rd}=\@gt;
	}
      }
    }
  }
  return \%tab_SNP;
}

#-------------------------------------------------
sub printOut {
  my ($self, $outfile)=@_;
  my $fout;
  my($log) = $self->{_logger};
  if(defined $outfile){
    unless (open(OUT,">$outfile")) {
        #$log->error("unable to open $outfile");
      return;
    }
    $fout=\*OUT;
  }
  else{
    $fout=\*STDOUT;
  }
  print  $fout '<BEGIN_PARAMETERS>'."\n";

  my @hs=@{$self->{_param}};
  foreach my $line (@hs) {
    print $fout "$line\n";
  }
  print  $fout "<END_PARAMETERS>\n\n\n";


  my $ncontigs=$self->{_ncontigs};
  if($ncontigs>0){
    my @contigs=@{$self->{_contigs}};
    foreach my $r_c (@contigs) {
      $r_c->printOut($fout);
    }
  }

  if(defined $outfile){
    close($fout);
  }
}

#TODO these methods are duplicated in Polyphred.pm, fix
sub get_pcr_name_from_read_name{
    my ($self, $read_name) = @_;
    my ($sequence_read) = GSC::Sequence::Read->get(trace_name => $read_name);
    unless ($sequence_read){
        warn "no sequence_read for $read_name";
        return undef;
    }
    my $pcr_name = $sequence_read->template_id;
    unless ($pcr_name){
        warn "no template_id for sequence_read of $read_name";
        return undef;
    }
    return $pcr_name;
}

sub get_sample_name_from_pcr_name{
    my ($self, $pcr_name) = @_;
    my $pcr = GSC::DNA->get(dna_name => $pcr_name);
    return unless $pcr;
    unless ($pcr){
        warn "no GSC::DNA for dna_name => $pcr";
        return undef;
    }
    my $genomic_dna = $pcr->get_first_ancestor_with_type('genomic dna');
    unless ($genomic_dna){
        warn "no genomic dna for GSC::DNA for $pcr";
        return undef;
    }
    my $sample_name = $genomic_dna->dna_name;
    unless ($sample_name){
        warn "no sample name from genomic dna from $pcr";
        return undef;
    }
}

sub get_gene_name_from_amplicon_tag{
    my ($self, $amplicon_tag) = @_;
    unless ($amplicon_tag){
        warn "no amplicon tag passed in!";
        return undef;
    }
    unless ($amplicon_tag->isa('GSC::Sequence::Tag::Amplicon')){
        warn "passed in tag is not a GSC::Sequence::Tag::Amplicon!";
        return undef;
    }
    my $gene_tag = $amplicon_tag->get_gene_tag;
    unless ($gene_tag){
        warn "no gene tag returned associated with amplicon tag";
        return undef;
    }
    my $gene = GSC::Gene->get(ref_seq_id => $gene_tag->seq_id);
    unless ($gene){
        warn "no gene returned associated with gene tag";
        return undef;
    }
    my $gene_name = $gene->gene_name;
    return $gene_name;
}        

sub get_reference_base_string_from_fasta {
    my $self = shift;
    my $fasta = shift;

    my $fasta_fh = IO::File->new($fasta);
    unless ($fasta_fh) {
        die "Could not get file handle from $fasta\n";
    }

    # reference base string is the second line of the fasta file, throw away the first line and return the second
    $fasta_fh->getline;
    my @lines = $fasta_fh->getlines;
    chomp @lines;
    
    my $reference_base_string = join('', @lines);
    return $reference_base_string;
}

#THIS METHOD RETURNS A SNP AND INDEL HASH WITH SUMMARIZED INFORMATION FOR EACH SNP, THIS IS CURRENTLY USED ONLY FOR THE COMBINE AND ANNOTATE VARIANTS PIPELINE.  IF THIS CHANGES TO BECOME MORE GENERIC, REPLACE AMPLICON REFERENCES WITH SOMETHING MORE APPROPRIATE
sub collate_sample_group_mutations {
    my $self = shift;
    my (@snps, @indels, $snp_hash, $indel_hash);
    
    my $fasta = $self->{fasta};
    my $reference_sequence = $self->get_reference_base_string_from_fasta($fasta);
    my $genomic_offset = $self->{begin_position};
    my $chromosome = $self->{chromosome};
    my $assembly_length = length($reference_sequence);
    

    for my $contig (@{$self->{_contigs}}){
        #skip contigs which aren't the "main contig"
        #if($contig->{_NAME} !~ /$assembly_name/ and $contig->{_NAME} ne "Contig1"){
        #    next;
        #}

#INITIAL FILTER
#we have to loop through each read here
#and the next step will group reads based on pcr_product

#SNPS
        my $ref = $contig->{_SNP_READSITES};
#composition of ref:
#HASH
#con_pos=> #   HASH
#   read_name => Polyscan::snp
#       _allel1
#       _allel2
#       _logger
#       _rdpos
#       _read
#       _refpos #TODO this may not work with padded assemblies
#       _score
#       _var_tag ( heterozygous )
        for my $con_pos (keys %$ref){
            my $ref_sequence = substr($reference_sequence, $con_pos - 1,1);
            for my $genotype (values %{$ref->{$con_pos}}){ 

                
                next unless defined $genotype->{_score}; #next if both alleles match reference
                next unless defined $genotype->{_var_tag}; #next if both alleles match reference

                #next unless $genotype->{_score} > 0;

               
                my $pcr_name = $self->get_pcr_name_from_read_name($genotype->{_read});
                warn "no pcr_name for read name: " . $genotype->{_read} and next unless $pcr_name;
                my $sample_name = $self->get_sample_name_from_pcr_name($pcr_name);
                warn "no sample for read name: " . $genotype->{_read} and next unless $sample_name;
                
                my $genomic_coord = ($genomic_offset + $con_pos -1);
                
                my $filtered_genotype;
                #positionally determined
                $filtered_genotype->{reference}         = uc($ref_sequence);
                $filtered_genotype->{chromosome}        = $chromosome;
                $filtered_genotype->{variation_type}    = 'SNP';
                $filtered_genotype->{read_type}         = 'sanger'; #TODO 
                
                #read name determined
                $filtered_genotype->{pcr_product_name}  = $pcr_name;
                $filtered_genotype->{sample_name}       = $sample_name;
                $filtered_genotype->{begin_position}             = $genomic_coord;
                $filtered_genotype->{end_position}              = $genomic_coord;
               
                #copied from genotype
                $filtered_genotype->{allele1}           = uc($genotype->{_allel1});
                $filtered_genotype->{allele2}           = uc($genotype->{_allel2});
                $filtered_genotype->{allele1_type}      = $genotype->{_allel1} eq $ref_sequence? 'REF' : 'SNP';
                $filtered_genotype->{allele2_type}      = $genotype->{_allel2} eq $ref_sequence? 'REF' : 'SNP';
                $filtered_genotype->{score}             = $genotype->{_score};
                $filtered_genotype->{read}              = $genotype->{_read};
                $filtered_genotype->{variation_tag}     = $genotype->{_var_tag};
                
                $filtered_genotype->{filename} = $self->{_polyscan};
                $filtered_genotype->{con_pos} = $con_pos;

                push @{ $snp_hash->{$genomic_coord}->{$pcr_name} }, $filtered_genotype;
            }
        }

#BUILD READ COUNTS
        for my $position ( keys %$snp_hash ){
            my %read_support_hash;
            my %product_support_hash;
            my $total_reads;  #these totals aren't currently used, maybe in the future.  If so, we'd combine these with the indel totals before adding a key to the hash-adukes
            my $total_pcr_products;
            my @pcr_product_genotypes;
            
            for my $pcr_product ( keys %{ $snp_hash->{$position} } ){
                #choose a genotype for each pcr_product
                my @genotypes = @{$snp_hash->{$position}->{$pcr_product}};

                my $highest_score= -5;  #setting this lower than 0 because some snps have -1 score, not sure if they should even be returned, but for now they will be-adukes
                my $pcr_product_genotype;
                
                for my $genotype (@genotypes){

                    $total_reads++;
                    
                    #track supported variant alleles by read
                    my @supported_variants;
                    if ($genotype->{allele1_type} eq 'SNP'){
                        push @supported_variants, $genotype->{allele1};
                    }elsif($genotype->{allele2_type} eq 'SNP'){
                        push @supported_variants, $genotype->{allele2};
                    }
                    for my $supported_variant (@supported_variants){
                        $read_support_hash{$supported_variant}++;
                    }

                    #find highest scoring snp for this pcr_product
                    my $score = $genotype->{score};
                    if ($score > $highest_score){
                        $pcr_product_genotype = $genotype;
                        $highest_score = $score;
                    }
                }
                
                $total_pcr_products++;

                #track supported variant alleles by pcr product
                my @pcr_product_supported_variants;
                if ($pcr_product_genotype->{allele1_type} eq 'SNP'){
                    push @pcr_product_supported_variants, $pcr_product_genotype->{allele1};
                }elsif($pcr_product_genotype->{allele2_type} eq 'SNP'){
                    push @pcr_product_supported_variants, $pcr_product_genotype->{allele2};
                }
                for my $pcr_product_supported_variant (@pcr_product_supported_variants){
                    $product_support_hash{$pcr_product_supported_variant}++;
                }
                
                #store chosen pcr_product genotype
                push @pcr_product_genotypes, $pcr_product_genotype;
            }

            for my $pcr_product_genotype ( @pcr_product_genotypes ){
                if ($pcr_product_genotype->{allele1_type} eq 'SNP'){
                    my $variant = $pcr_product_genotype->{allele1};
                    my $read_count = $read_support_hash{$variant};
                    my $pcr_count = $product_support_hash{$variant};
                    $pcr_product_genotype->{allele1_read_support} = $read_count;
                    $pcr_product_genotype->{allele1_pcr_product_support} = $pcr_count;
                }
                if ($pcr_product_genotype->{allele2_type} eq 'SNP'){
                    my $variant = $pcr_product_genotype->{allele2};
                    my $read_count = $read_support_hash{$variant};
                    my $pcr_count = $product_support_hash{$variant};
                    $pcr_product_genotype->{allele2_read_support} = $read_count;
                    $pcr_product_genotype->{allele2_pcr_product_support} = $pcr_count;
                }
                $pcr_product_genotype->{allele1_read_support} ||= 0;
                $pcr_product_genotype->{allele2_read_support} ||= 0;
                $pcr_product_genotype->{allele1_pcr_product_support} ||= 0;
                $pcr_product_genotype->{allele2_pcr_product_support} ||= 0;
            }

            push @snps, @pcr_product_genotypes
        }

#INDELS
#INITIAL FILTER, actually _INDEL_READSITES is already filtered to just
#reads w/ indels. we are just appending some info
        my $indel_ref = $contig->{_INDEL_READSITES};
#structure = HASH
#position ->read Polyscan::indel
#       '_comment' => undef
#       '_logger' => undef
#       '_rdpos' => 247
#       '_read' => 'H_GP-0113tPCR0003202_02ta.g1'
#       '_refpos' => 184  #TODO this may not work with padded assemblies
#       '_refpos2' => undef
#       '_score' => 39
#       '_seqindel' => 'G'
#       '_size' => 1
#       '_spratio' => 0.53
#       '_type' => 'hetDel'

        for my $con_pos (keys %$indel_ref){

            my $ref_sequence = substr($reference_sequence,$con_pos - 1,1);
            for my $genotype(values %{$indel_ref->{$con_pos}}){
                
                next unless defined $genotype->{_score};
                #next unless $genotype->{_score} > 0;
                
                my $pcr_name = $self->get_pcr_name_from_read_name( $genotype->{_read} );
                warn "no pcr_name for read name: " . $genotype->{_read} and next unless $pcr_name;
                my $sample_name = $self->get_sample_name_from_pcr_name($pcr_name);
                warn "no sample for read name: " . $genotype->{_read} and next unless $sample_name;

                my $genomic_coord = ($genomic_offset + $con_pos -1);
                
                my $filtered_genotype;

                if ($genotype->{_type} =~ /del/i) {
                    $filtered_genotype->{variation_type}='DEL';
                }
                elsif ($genotype->{_type} =~ /ins/i) {
                    $filtered_genotype->{variation_type}='INS';
                }
                else {
                    warn "Indel variation type: ".$genotype->{_type}. "cannot be matched...skipping\n";
                    print Dumper($genotype);
                    next;
                }
                $filtered_genotype->{allele1_type} = 'REF';  #TODO this is a change, but i think it's more accurate.-adukes
                $filtered_genotype->{allele2_type} = $filtered_genotype->{variation_type};
                $filtered_genotype->{read_type}         = 'sanger'; #TODO 
                
                $filtered_genotype->{pcr_product_name}= $pcr_name;
                $filtered_genotype->{sample_name}= $sample_name;
                $filtered_genotype->{chromosome} = $chromosome;
                
                $filtered_genotype->{length} = $genotype->{_size};
                $filtered_genotype->{begin_position} = $genomic_coord;
                $filtered_genotype->{end_position} = $genomic_coord + $filtered_genotype->{length};
            
                $filtered_genotype->{reference} = uc($ref_sequence);
                $filtered_genotype->{variation_tag} = $genotype->{_type};
                $filtered_genotype->{allele1} = uc(substr($reference_sequence,($con_pos - 1), $genotype->{_size}));
                $filtered_genotype->{allele2} = uc($genotype->{_seqindel});
                $filtered_genotype->{score} = $genotype->{_score};
                $filtered_genotype->{read} = $genotype->{_read};

                $filtered_genotype->{filename} = $self->{_polyscan};
                $filtered_genotype->{con_pos} = $con_pos;

                push @{$indel_hash->{$genomic_coord}->{$pcr_name}}, $filtered_genotype;
            }
        }

#INDELS
#BUILD READ COUNTS
        for my $position ( keys %$indel_hash ){
            my %read_support_deletion_hash;
            my %read_support_insertion_hash;
            my %product_support_insertion_hash;
            my %product_support_deletion_hash;
            my $total_reads;
            my $total_pcr_products;
            my @pcr_product_genotypes;
            
            for my $pcr_product ( keys %{ $indel_hash->{$position} } ){
                #choose a genotype for each pcr_product
                my @genotypes = @{$indel_hash->{$position}->{$pcr_product}};

                my $highest_score=-5;
                my $pcr_product_genotype;

                for my $genotype (@genotypes){
                    $total_reads++;

                    #track supported variant alleles by read
                    my @supported_insertions;
                    my @supported_deletions;
                    if($genotype->{allele2_type} eq 'INS'){
                        push @supported_insertions, $genotype->{allele2};
                    }elsif($genotype->{allele2_type} eq 'DEL'){
                        push @supported_deletions, $genotype->{allele2};
                    }

                    for my $supported_deletion (@supported_deletions){
                        $read_support_deletion_hash{$supported_deletion}++;
                    }
                    for my $supported_insertion (@supported_insertions){
                        $read_support_insertion_hash{$supported_insertion}++;
                    }

                    #find highest scoring indel for this pcr_product
                    my $score = $genotype->{score};
                    if ($score > $highest_score){
                        $pcr_product_genotype = $genotype;
                        $highest_score = $score;
                    }
                }

                $total_pcr_products++;

                #track supported variant alleles by pcr product
                my @pcr_product_supported_insertions;
                my @pcr_product_supported_deletions;
                if ($pcr_product_genotype->{allele2_type} eq 'INS'){
                    push @pcr_product_supported_insertions, $pcr_product_genotype->{allele2};
                }elsif($pcr_product_genotype->{allele2_type} eq 'DEL'){
                    push @pcr_product_supported_deletions, $pcr_product_genotype->{allele2};
                }
                for my $pcr_product_supported_insertion (@pcr_product_supported_insertions){
                    $product_support_insertion_hash{$pcr_product_supported_insertion}++;
                }
                for my $pcr_product_supported_deletion (@pcr_product_supported_deletions){
                    $product_support_deletion_hash{$pcr_product_supported_deletion}++;
                }

                #store chosen pcr_product genotype
                push @pcr_product_genotypes, $pcr_product_genotype;
            }

            for my $pcr_product_genotype ( @pcr_product_genotypes ){
                if ($pcr_product_genotype->{allele2_type} eq 'INS'){
                    my $insertion = $pcr_product_genotype->{allele2};
                    my $read_count = $read_support_insertion_hash{$insertion};
                    my $pcr_count = $product_support_insertion_hash{$insertion};
                    $pcr_product_genotype->{allele2_read_support} = $read_count;
                    $pcr_product_genotype->{allele2_pcr_product_support} = $pcr_count;
                }
                if ($pcr_product_genotype->{allele2_type} eq 'DEL'){
                    my $deletion = $pcr_product_genotype->{allele2};
                    my $read_count = $read_support_deletion_hash{$deletion};
                    my $pcr_count = $product_support_deletion_hash{$deletion};
                    $pcr_product_genotype->{allele2_read_support} = $read_count;
                    $pcr_product_genotype->{allele2_pcr_product_support} = $pcr_count;
                }
                $pcr_product_genotype->{allele1_read_support} = 0;
                $pcr_product_genotype->{allele1_pcr_product_support} = 0;

                $pcr_product_genotype->{allele2_read_support} ||= 0;
                $pcr_product_genotype->{allele2_pcr_product_support} ||= 0;
            }
            push @indels, @pcr_product_genotypes
        }
    }

    $self->{_snps_array}=\@snps;
    $self->{_indels_array}=\@indels;

    return (\@snps,\@indels); 
}

1;

#-------------------------------------------------

=head2 EXPORT

None by default.

=head1 SEE ALSO

YAML
Log::Log4perl
Log::MedSeq
MG::IO::Polyphred

=head1 FILES

Example

=head1 BUGS

Maybe you'll find some. Let me know.

=head1 AUTHOR

Ken Chen, E<lt>kchen@watson.wustl.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007 by Ken Chen.  All rights reserved.

=cut

__END__

1;
# $Header$
