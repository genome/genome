package Genome::Model::Tools::Somatic::PindelStats;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use IO::File;

my %positions;
my %insertions;
my %deletions;
my %size_type_hist;
my %tumor_support_hist;
my %ref_support_hist;
my %event_support_ratio;
my $dbsnp_hits;

class Genome::Model::Tools::Somatic::PindelStats {
    is => 'Command',
    has => [
        indels_all_sequences_bed_file =>{
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Indel sites to assemble in annotator input format',
        },
        pindel_output_directory => {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => "location of the pindel output_directory.",
        },
        refseq =>{
            type => 'String',
            is_optional => 1,
            default => Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fasta',
            doc => "reference sequence to use for reference assembly",
        },
        use_old_pindel => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Run on pindel 0.2 or 0.1',
        },
        dbsnp_concordance => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Set this to cause dbsnp concordance to be calculated for each event.',
        },
        _dbsnp_insertions => {
            type => 'String',
            is_optional => 1,
            default => '/gscmnt/ams1102/info/info/dbsnp130_indels/insertions_start_stop_adjusted_dbsnp130',
            doc => 'dbsnp insertion file',
        },
        _dbsnp_deletions => {
            type => 'String',
            is_optional => 1,
            default => '/gscmnt/ams1102/info/info/dbsnp130_indels/deletions_adjusted_dbsnp130',
            doc => 'dbsnp deletion file',
        },
    ]
};

sub execute {
    my $self = shift;
    $| = 1;
    #my $file = $self->indels_all_sequences_bed_file;
    my $dir = $self->pindel_output_directory;
    #my $reference_fasta = $self->refseq;

    #my $fh = IO::File->new($file);

    #my %indels;
    #my %answers;
    if($self->dbsnp_concordance){
        my $ifh = IO::File->new($self->_dbsnp_insertions);
        while (my $line = $ifh->getline) {
            chomp $line;
            my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
            next unless ($allele =~ m/-/);
            $allele = substr($allele, 2);
            $insertions{$chr}{$start}{$stop}{'allele'}=$allele;
            $insertions{$chr}{$start}{$stop}{'id'}=$id;
        }
        $ifh->close;
        my $dfh = IO::File->new($self->_dbsnp_deletions);
        while (my $line = $dfh->getline) {
            chomp $line;
            my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
            next unless ($allele =~ m/-/);
            $allele = substr($allele, 2);
            $deletions{$chr}{$start}{$stop}{'allele'}=$allele;
            $deletions{$chr}{$start}{$stop}{'id'}=$id;
        }
        $dfh->close;
    }
    #my $chr = 1;
    #print "CHR\tSTART\tSTOP\tREF\tVAR\tINDEL_SUPPORT\tREFERENCE_SUPPORT\t%+STRAND\tDBSNP_ID\n";
    my @chromosomes = qw| 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y|;
    for my $chr (@chromosomes){
        $self->process_file($chr);
    }

}
    #for my $size (sort {$a <=> $b} (keys(%{$tumor_support_hist{'D'}}))){
    #    
    #}

=cut
    print "Histogram of event size for Deletions\n";
    $self->display_histogram($size_type_hist{"D"});
    print "\n\n";
    print "Histogram of event size for Insertions\n\n";
    $self->display_histogram($size_type_hist{"I"});
    my %histogram;
    for my $type (sort(keys(%size_type_hist))){
        my $total = 0;
        for my $size (sort {$a <=> $b} (keys(%{$size_type_hist{$type}}))){
            $histogram{$size}+= $size_type_hist{$type}{$size};
        }
    }
    print "\n\n";
    print "Histogram of event size for both INS/DEL\n";
    my $total_indels = $self->display_histogram(\%histogram);
    
    print "\n\n";
    print "Histogram of read support in tumor for Deletions\n";
    $self->display_histogram($tumor_support_hist{'D'});
    print "\n\n";
    print "Histogram of read support in tumor for Insertions\n";
    $self->display_histogram($tumor_support_hist{'I'});
    my %hgram;
    for my $type (sort(keys(%tumor_support_hist))){
        my $total = 0;
        for my $size (sort {$a <=> $b} (keys(%{$tumor_support_hist{$type}}))){
            $hgram{$size}+= $tumor_support_hist{$type}{$size};
        }
    }
    print "\n\n";
    print "Histogram of read support in tumor for both INS/DEL\n";
    $self->display_histogram(\%hgram);
    $self->bin_hash(\%hgram,(4,9));
    return 1;

    print "\n\n";
    print "Histogram of reads supporting reference for Deletions\n";
    $self->display_histogram($ref_support_hist{'D'});
    print "\n\n";
    print "Histogram of reads supporting reference for Insertions\n";
    $self->display_histogram($ref_support_hist{'I'});
    my %hgrammer;
    for my $type (sort(keys(%ref_support_hist))){
        my $total = 0;
        for my $size (sort {$a <=> $b} (keys(%{$ref_support_hist{$type}}))){
            $hgrammer{$size}+= $ref_support_hist{$type}{$size};
        }
    }
    print "\n\n";
    print "Histogram of reads supporting reference for both INS/DEL\n";
    $self->display_histogram(\%hgrammer);
    if($self->dbsnp_concordance){
        print "dbSNP concordance =  ( dbSNP hits:".$dbsnp_hits." / total indels:".$total_indels.") = ".($dbsnp_hits/$total_indels)."\n";
    }
    print "\n\n";
    print "Histogram of ratios of variant support vs ref support reads\n";
    $self->display_histogram(\%event_support_ratio);
    
}

    
    print "Histogram of event sizes for all events\n";
    print "size\toccurence\n";
    print "=============================\n";
    my $total;
    for my $size (sort {$a <=> $b} (keys(%histogram))){
        print $size."\t".$histogram{$size}."\n";
        $total+=$histogram{$size};
    }
    print "=============================\n";
    print "Total events: ".$total."\n\n";

    $total = 0;
    print "Histogram of # of reads supporting tumor event\n";
    print "# of reads\toccurence\n";
    print "=============================\n";
    for my $number (sort {$a <=> $b} (keys(%tumor_support_hist))){
        print $number."\t".$tumor_support_hist{$number}."\n";
        $total+= $tumor_support_hist{$number};
    }
    print "=============================\n";
    print "Total events: ".$total."\n\n";

}
=cut

sub display_histogram {
    my $self = shift;
    my $histogram = shift;
    my $total = 0;
    for my $index(sort {$a <=> $b} (keys(%{$histogram}))){
        print $index."\t".$histogram->{$index}."\n";
        $total+=$histogram->{$index};
        #$histogram{$size}+= $size_type_hist{$type}{$size};
    }
    print "============================\n";
    print "Total events: ".$total."\n";
    return $total;
}

sub bin_hash {
    my $self = shift;
    my $hash = shift;
    my %hash = %{$hash};
    my @bins = @_;
    for my $b (@bins){
        print "bin size = ".$b."\n";
    }
    my %bins;
    for my $val (sort {$a <=> $b} (keys(%hash))){
        if($val>$bins[-1]){
            $bins{$bins[-1]+1}+=$hash{$val};
        }else{
            for my $index(0..scalar(@bins)){
                if($val < $bins[$index]){
                    #if($index==1){
                    #    $bins{0}+=$hash{$val};
                    #    last;
                    #}else{
                        $bins{$bins[$index]}+=$hash{$val};
                        last;
                   # }
                }
            }
        }
    }
    print "Bins\n";
    $self->display_histogram(\%bins); 
    return \%bins;
}

sub process_file {
    my $self = shift;
    my $chr = shift;
    my $dir = $self->pindel_output_directory;
    my $reference_fasta = $self->refseq;
    my $filename = $dir."/".$chr."/indels_all_sequences";
    unless(-s $filename){
        print $filename." had zero size, skipping.\n";
        return;
    }
    my $pindel_output = IO::File->new($filename);
    my $pindel_config = $dir."/".$chr."/pindel.config";
    my $pconf = IO::File->new($pindel_config);
    $pconf->getline;
    my $tumor_bam = $pconf->getline;
    ($tumor_bam) = split /\s/, $tumor_bam;
    unless(-s $tumor_bam){
        die "couldnt find tumor bam reference in pindel.config at ".$tumor_bam."\n";
    }
    my %events;
    my ($chrom,$pos,$size,$type);

    while(my $line = $pindel_output->getline){
        my $normal_support=0;
        my $read = 0;
        if($line =~ m/^#+$/){
            my $call = $pindel_output->getline;
            if($call =~ m/^Chr/){
                while(1) {
                    $line = $pindel_output->getline;
                    if($line =~ m/#####/) {
                        $call = $pindel_output->getline;
                    } 
                    if($call !~ m/^Chr/) {
                        last;
                    }          
                }
            }
            my $reference = $pindel_output->getline;
            my @call_fields = split /\s/, $call;
            my $type = $call_fields[1];
            my $size = $call_fields[2];   #12
            my $pos_strand = 0;
            my $neg_strand = 0;
            my $mod = ($call =~ m/BP_range/) ? 2: -1;
            my $support;
            if($self->use_old_pindel){
                $support = ($type eq "I") ? $call_fields[10+$mod] : $call_fields[12+$mod];
            } else {
                $support = $call_fields[12+$mod];
            }
            unless(defined($support)){
                print "No support. Call was:   ".$call."\n";
                die;
            }
            for (1..$support){
                $line = $pindel_output->getline;
                if($line =~ m/normal/) {
                    $normal_support=1;
                }
                if($line =~ m/\+/){
                    $pos_strand++;
                }else{
                    $neg_strand++;
                }
                $read=$line;
            }

#charris speed hack
            my ($chr,$start,$stop);
            if($self->use_old_pindel){
                $chr = ($type eq "I") ? $call_fields[4] : $call_fields[6];
                $start= ($type eq "I") ? $call_fields[6] : $call_fields[8];
                $stop = ($type eq "I") ? $call_fields[7] : $call_fields[9];
            } else {
                $chr = $call_fields[6];
                $start= $call_fields[8];
                $stop = $call_fields[9];
            }
            if($type eq 'I') {
                $start = $start - 1;
            }            
#charris speed hack

            my $type_and_size = $type."/".$size;
            $events{$chr}{$start}{$type_and_size}{'neg'}=$neg_strand;
            $events{$chr}{$start}{$type_and_size}{'pos'}=$pos_strand;
            if($normal_support){
                $events{$chr}{$start}{$type_and_size}{'normal'}=$normal_support;
            }else{
                if($self->dbsnp_concordance){
                    $events{$chr}{$start}{$type_and_size}{'var'}=$self->get_variant_allele($call, $reference, $read);
                } else {
                    $events{$chr}{$start}{$type_and_size}{'var'}=1;
                }
            }
        }
    }
    for $chrom (sort {$a cmp $b} (keys(%events))){
        for $pos (sort{$a <=> $b} (keys( %{$events{$chrom}}))){
            for my $type_and_size (sort(keys( %{$events{$chrom}{$pos}}))){
                unless(exists($events{$chrom}{$pos}{$type_and_size}{'normal'})){
                    my $pos_strand = $events{$chrom}{$pos}{$type_and_size}{'pos'};
                    my $neg_strand = $events{$chrom}{$pos}{$type_and_size}{'neg'};
                    if($self->dbsnp_concordance){
                        my $var = $events{$chrom}{$pos}{$type_and_size}{'var'};
                        
                        my $dbsnp = '-';
                        if(defined($var)){
                            $dbsnp = $self->dbsnp_lookup($var);
                        }
                        if($dbsnp ne '-'){
                            $dbsnp_hits++;
                        }
                    }
                    my $pos_percent=0;
                    if($neg_strand==0){
                        $pos_percent = 1.0;
                    } else {
                        $pos_percent = sprintf("%.2f", $pos_strand / ($pos_strand + $neg_strand));
                    }
                    my $answer = "neg = ".$neg_strand."\tpos = ".$pos_strand." which gives % pos str = ".$pos_percent."\n";
                    my $reads = $pos_strand + $neg_strand;
                    my ($type,$size) = split /\//, $type_and_size;
                    #print "checking on ".$chrom."\t".$pos."\t".$type_and_size."\n";
                    #$tumor_support_hist{$type}{'reads'}{$reads}++;
                    $tumor_support_hist{$type}{$size}{$reads}++;
                    $size_type_hist{$type}{$size}++;
                    
                    #my @stop = keys(%{$positions{$chrom}{$pos}});
                    #unless(scalar(@stop)==1){
                    #    
                    #    die "too many stop positions at ".$chrom." ".$pos."\n";
                    #}
                    my $stop = ($type eq 'I') ? $pos+2 : $pos + $size;
=cut
                    if($size > 100){
                        next;
                    }

                    my @results = `samtools view $tumor_bam $chrom:$pos-$stop | grep -v "XT:A:M"`;
                    my $read_support=0;
                    for my $result (@results){
                        #print $result;
                        chomp $result;
                        my @details = split /\t/, $result;
                        if($result =~ /NM:i:(\d+)/){
                            if($1 > 2){
                                next;
                            }
                        }
                        unless($details[5] =~ m/[IDS]/){
                            if(($details[3] > ($pos - 40))&&($details[3] < ($pos -10))){
                                $read_support++;
                                #print "cigar = ".$details[5]."\n";
                            }
                        }
                    }
                    $ref_support_hist{$type}{$read_support}++;
                    my $ratio = $reads/($reads+$read_support);
                    $ratio = int($ratio * 100);
                    $event_support_ratio{$ratio/100}++;
                    #my $dbsnp_id = $self->dbsnp_lookup($events{$chrom}{$pos}{$type_and_size}{'bed'});
                    #my $output = $events{$chrom}{$pos}{$type_and_size}{'bed'}."\t".$reads."\t".$read_support."\t".$pos_percent."\t$dbsnp_id\n";
                    #print $output;
=cut
                }
            }
        }
    }
}

sub dbsnp_lookup {
    my $self=shift;
    my $bed_line =shift;
    my $dbsnp_id="-";
    chomp $bed_line;
    my ($chr, $start, $stop, $type,$var) = split "\t", $bed_line;
    if($type eq "I") {
        if(exists($insertions{$chr}{$start}{$stop}{'allele'})) {
            if ($var eq $insertions{$chr}{$start}{$stop}{'allele'}) {
                $dbsnp_id=$insertions{$chr}{$start}{$stop}{'id'};
            }
        }
    }
    else {        
        if(exists($deletions{$chr}{$start}{$stop}{'allele'})) {
            #if ($ref eq $deletions{$chr}{$start}{$stop}{'allele'}) {
            $dbsnp_id=$deletions{$chr}{$start}{$stop}{'id'};
           # }
        } 
    }
    return $dbsnp_id;
}


sub get_variant_allele {
    my $self = shift;
    my ($call, $reference, $first_read) = @_;
    #parse out call bullshit
    chomp $call;
    my @call_fields = split /\s+/, $call;
    my $type = $call_fields[1];
    my $size = $call_fields[2];
    $DB::single=1 if $size == 1445;
    my ($chr,$start,$stop);
    if($self->use_old_pindel){
        $chr = ($type eq "I") ? $call_fields[4] : $call_fields[6];
        $start= ($type eq "I") ? $call_fields[6] : $call_fields[8];
        $stop = ($type eq "I") ? $call_fields[7] : $call_fields[9];
    } else {
        $chr = $call_fields[6];
        $start= $call_fields[8];
        $stop = $call_fields[9];
    }
    my $support = $call_fields[-1];
    my ($ref, $var);
    if($type =~ m/D/) {
        $var =0;
    }
    elsif($type =~ m/I/) {
        $start = $start - 1;
        $ref=0;
        my ($letters_until_space) =   ($reference =~ m/^([ACGTN]+) /);
        my $offset_into_first_read = length($letters_until_space);
        $var = substr($first_read, $offset_into_first_read, $size);
    }
    if($size >= 100) {
        #my $big_fh = $self->_big_output_fh;
        #$big_fh->print("$chr\t$start\t$stop\t$size\t$support\n");
        return undef;
    }
    return join("\t",($chr,$start,$stop,$type,$var));
}

sub parse {
    my $self = shift;
    my $reference_fasta = $self->refseq;
    my ($call, $reference, $first_read) = @_;
    #parse out call bullshit
    chomp $call;
    my @call_fields = split /\s+/, $call;
    my $type = $call_fields[1];
    my $size = $call_fields[2];
    $DB::single=1 if $size == 1445;
    my ($chr,$start,$stop);
    if($self->use_old_pindel){
        $chr = ($type eq "I") ? $call_fields[4] : $call_fields[6];
        $start= ($type eq "I") ? $call_fields[6] : $call_fields[8];
        $stop = ($type eq "I") ? $call_fields[7] : $call_fields[9];
    } else {
        $chr = $call_fields[6];
        $start= $call_fields[8];
        $stop = $call_fields[9];
    }
    my $support = $call_fields[-1];
    my ($ref, $var);
    if($type =~ m/D/) {
        $var =0;
        ###Make pindels coordinates(which seem to be last undeleted base and first undeleted base) 
        ###conform to our annotators requirements

        ###also deletions which don't contain their full sequence should be dumped to separate file
        $stop = $stop - 1;
        my $allele_string;
        my $start_for_faidx = $start+1;
        my $sam_default = Genome::Model::Tools::Sam->path_for_samtools_version;
        my $faidx_cmd = "$sam_default faidx " . $reference_fasta . " $chr:$start_for_faidx-$stop";
        #my $faidx_cmd = "$sam_default faidx " . $reference_fasta . " $chr:$start-$stop";
        my @faidx_return= `$faidx_cmd`;
        shift(@faidx_return);
        chomp @faidx_return;
        $allele_string = join("",@faidx_return);

        $ref = $allele_string;
    }
    elsif($type =~ m/I/) {
        $start = $start - 1;
        $ref=0;
        my ($letters_until_space) =   ($reference =~ m/^([ACGTN]+) /);
        my $offset_into_first_read = length($letters_until_space);
        $var = substr($first_read, $offset_into_first_read, $size);
    }
    if($size >= 100) {
        #my $big_fh = $self->_big_output_fh;
        #$big_fh->print("$chr\t$start\t$stop\t$size\t$support\n");
        return undef;
    }
    return ($chr,$start,$stop,$ref,$var);
}
