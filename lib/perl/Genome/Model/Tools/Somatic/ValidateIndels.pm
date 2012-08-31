package Genome::Model::Tools::Somatic::ValidateIndels;

use strict;
use warnings;

use Genome;
use IO::File;

#my %positions;
#my %insertions;
#my %deletions;
my %ranges;

my $use_ranges=0;

class Genome::Model::Tools::Somatic::ValidateIndels {
    is => 'Command',
    has => [
        valid_indel_bed_file =>{
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Bed format file containing validated true events',
        },
        valid_indel_annotation_file =>{
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Annotation format file containing validated true events',
        },
        called_indel_bed_file => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => "Bed format file containing indels as called by the indel caller being tested",
        },
        called_indel_annotation_file => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => "Annotation format file containing indels as called by the indel caller being tested",
        },
        fuzzy_output_file => {
            type => 'String',
            is_optional => 1,
            is_output => 1,
            doc => "The potential fuzzy matches will be dumped here",
        },
        missed_output_file => {
            type => 'String',
            is_optional => 1,
            is_output => 1,
            doc => "The calls which missed entirely will be dumped here",
        },
        allele_mismatch_output_file => {
            type => 'String',
            is_optional => 1,
            is_output => 1,
            doc => "The calls which hit positionally but had different alleles will be dumped here",
        },
        intersection_range => {
            type => 'Integer',
            is_input => 1,
            is_optional => 0,
            default => -1,
            doc => "Use this to specify the range within to look for fuzzy concordance. Default is -1, which performs no fuzzy concordance.",
        },
        refseq =>{
            type => 'String',
            is_optional => 1,
            default => Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fasta',
            doc => "reference sequence to use for reference assembly",
        },
    ],
    has_param => [
         lsf_queue => {
             default_value => 'tcga',
         }, 
         lsf_resource => {
             default_value => "-M 6000000 -R 'select[type==LINUX64 && mem>16000] rusage[mem=16000]'",
         },
     ],
};

#Internally, this tool uses bed format to represent chr start stop of events.
sub load_bed_file {
    my $self = shift;
    my $file = shift;
    my %data;

    my $fh = Genome::Sys->open_file_for_reading($file);
    while(my $line = $fh->getline){
        chomp $line;
        my ($chr,$start,$stop,$ref,$var,$lower_range,$upper_range) = split "\t",$line;

        if($ref =~ m/\//){
            $upper_range = $lower_range;
            $lower_range = $var;
            ($ref,$var) = split "/", $ref;
        }
        my $refvar = "$ref/$var";
        $data{$chr}{$start}{$stop}=$refvar;
        if(defined($lower_range)){
            $use_ranges=1;
            my $ranger = $lower_range.",".$upper_range;
            $data{$chr}{$start}{ranges}=$ranger;
        }
    }
    $fh->close;
    return \%data;
}

sub load_annotation_file {
    my $self = shift;
    my $file = shift;
    my %data;

    my $fh = Genome::Sys->open_file_for_reading($file);
    while(my $line = $fh->getline){
        chomp $line;
        my ($chr,$start,$stop,$ref,$var,$lower_range,$upper_range) = split "\t",$line;
        if($ref =~ m/\//){
            $upper_range = $lower_range;
            $lower_range = $var;
            ($ref,$var) = split "/", $ref;
        }
        if(defined($lower_range)){
            $use_ranges = 1;
        }
        if($var ne '0'){
            $start++;
            $lower_range++ if $lower_range;
        } else {
            $start--;
            $lower_range-- if $lower_range;
        }
        my $refvar = "$ref/$var";
        $data{$chr}{$start}{$stop}=$refvar;
        if(defined($lower_range)){
            my $ranger = $lower_range.",".$upper_range if $lower_range;
            $data{$chr}{$start}{ranges}=$ranger if $lower_range;
        }
    }
    $fh->close;
    return \%data;
}
sub execute {
    my $self = shift;

    unless(defined($self->valid_indel_bed_file) xor defined($self->valid_indel_annotation_file)){
        $self->error_message("You must either specify valid_indel_bed_file or valid_indel_annotation_file, but not both.");
        die $self->error_message;
    }
    unless(defined($self->called_indel_bed_file) xor defined($self->called_indel_annotation_file)){
        $self->error_message("You must either specify called_indel_bed_file or called_indel_annotation_file, but not both.");
        die $self->error_message;
    }
    
    my %valid_data = %{(defined($self->valid_indel_bed_file))? $self->load_bed_file($self->valid_indel_bed_file): $self->load_annotation_file($self->valid_indel_annotation_file)};
    my %called_data = %{(defined($self->called_indel_bed_file))? $self->load_bed_file($self->called_indel_bed_file): $self->load_annotation_file($self->called_indel_annotation_file)};

    my $valid_total = 0;
    my $called_total = 0;
    my $total_fail=0;
    my $fuzzy_hits=0;
    my $fuzzy_hits_exact_allele=0;
    my $exact_hits=0;
    my $within_range=0;
    my $allele_mismatch_hits=0;
    my %allele_mismatch;
    my %perfect_match;
    my %fuzzy_candidates;
    my %total_fail;
    my %touched_valids;
    for my $chr (sort(keys(%called_data))){
        for my $start (sort(keys(%{$called_data{$chr}}))){
            for my $stop (sort(keys(%{$called_data{$chr}{$start}}))){
                if($stop eq 'ranges'){ next;}
                $called_total++;
                if(exists($valid_data{$chr}{$start})&&defined($valid_data{$chr}{$start}{$stop})){
                    if($valid_data{$chr}{$start}{$stop} eq $called_data{$chr}{$start}{$stop}){
                        $perfect_match{$chr}{$start}{$stop} = $called_data{$chr}{$start}{$stop};
                        $touched_valids{$chr}{$start}{$stop} = 1;
                        $exact_hits++;
                    }else{
                        $allele_mismatch{$chr}{$start}{$stop} = 1;
                        $touched_valids{$chr}{$start}{$stop} = 1;
                    
                        $allele_mismatch_hits++;
                    }
                } else {
                    unless($self->intersection_range == -1){
                        my $fuzzy = $self->find_intersection(join(",",($chr,$start,$stop)),\%valid_data);
                        unless(defined($fuzzy)){
                            $total_fail{$chr}{$start}{$stop} = $called_data{$chr}{$start}{$stop};
                            $total_fail++;
                            next;
                        }
                        my ($ch,$st,$sp) = split ",",$fuzzy;
                        $fuzzy_candidates{$chr}{$start}{$stop} = join("\t",($called_data{$chr}{$start}{$stop},$ch,$st,$sp,$valid_data{$ch}{$st}{$sp}));
                        if($called_data{$chr}{$start}{$stop} eq $valid_data{$ch}{$st}{$sp}){
                            $fuzzy_hits_exact_allele++;
                        } else {
                            $fuzzy_hits++;
                        }
                        $touched_valids{$ch}{$st}{$st} = 1;
                        if($use_ranges){
                            my ($lower_range, $upper_range) = split ",", $called_data{$chr}{$start}{'ranges'};
                            if($self->inside($chr,$lower_range,$upper_range,$ch,$st,$sp)){
                                $perfect_match{$chr}{$start}{$stop} = $called_data{$chr}{$start}{$stop};
                                $within_range++;
                            }
                        }
                    }
                }
            }
        }
    }
    for my $chr (keys(%valid_data)){
        for my $start (keys(%{$valid_data{$chr}})){
            $valid_total++;
        }
    }
    my $touched=0;
    for my $chr (keys(%touched_valids)){
        for my $start (keys(%{$touched_valids{$chr}})){
            $touched++;
        }
    }

    print "The total hits in results =                                 ".($exact_hits+$within_range+$allele_mismatch_hits)."\n";
    print "    Of which, the total EXACTLY matching hits was =         ".$exact_hits."\n";
    print "    Of which, the total within bp range =                   ".$within_range."\n" if $use_ranges;
    print "    Of which, Alleles different, but positional exact hit = ".$allele_mismatch_hits."\n";
    print "Total Calls Read = ".$called_total."\n";
    print "True Positive Hit Ratio = ".(int(($exact_hits/$called_total)*100)/100)."\n";
    #print "For a hits ratio of ".(int(($exact_hits/$valid_total)*100)/100)."\n";
    print "Total validation events = ".$valid_total."\n";
    print "Fuzzy Hits, within a range of ".$self->intersection_range." = ".$fuzzy_hits."\n" unless $self->intersection_range==-1;
    print "Fuzzy Hits, within a range of ".$self->intersection_range.", but having exactly matching allels = ".$fuzzy_hits_exact_allele."\n" unless $self->intersection_range==-1;
    print "Calls which had a BP Range that included the valid indel = ".$within_range."\n" if $use_ranges;
    print "Fuzzy Hits - BP Range hits = ".(($fuzzy_hits_exact_allele+$fuzzy_hits)-$within_range)."\n" unless $self->intersection_range==-1;
    print "Complete misses = ".$total_fail."\n" unless $self->intersection_range ==-1;

    print "Total exclusive number of valid indels touched = ".$touched."\n";

    #print Data::Dumper::Dumper(%fuzzy_candidates);
    $self->dump_to_file(\%fuzzy_candidates,$self->fuzzy_output_file) if $self->fuzzy_output_file;
    $self->dump_to_file(\%total_fail,$self->missed_output_file) if $self->missed_output_file;
    $self->dump_to_file(\%allele_mismatch,$self->allele_mismatch_output_file) if $self->allele_mismatch_output_file;

    return 1;
}

sub dump_to_file {
    my $self = shift;
    my $hash = shift;
    my $file = shift;
    my $fh = Genome::Sys->open_file_for_writing($file);
    for my $chr (sort(keys(%{$hash}))){
        for my $start (sort {$a <=> $b} keys(%{$hash->{$chr}})){
            for my $stop (keys(%{$hash->{$chr}{$start}})){
                my $stuff = join("\t",($chr,$start,$stop,$hash->{$chr}{$start}{$stop}))."\n";
                print $fh $stuff;
            }
        }
    }
    $fh->close;
}


sub find_intersection {
    my $self = shift;
    my ($chr,$start,$stop,$ref,$var) = split ",",shift;
    my $hash = shift;
    my $range = $self->intersection_range;
    $start = $start - $range;
    $stop = $stop + $range;
    unless(defined($hash->{$chr})){
        return undef;
    }
    my @answers;
    for my $st (sort {$a <=> $b} keys(%{$hash->{$chr}})){
        my @sp = keys(%{$hash->{$chr}{$st}});
        unless(scalar(@sp)==1){
            #return undef;
        }
        my $sp = $sp[0];
        my $ans=0;
        if(($start == $st)||($stop==$st)){
            $ans = 1;
        } elsif (($start == $sp)||($stop==$sp)){
            $ans = 1;
        } elsif (($st>$start)&&($st<$stop)){
            $ans = 1;
        } elsif (($sp<$stop)&&($sp>$start)){
            $ans = 1;
        }
        if($ans==1){
            push @answers, join(",",($chr,$st,$sp));
        }
    }
    if(scalar(@answers)>1){
        $self->error_message("Found multiple intersections for $chr:$start-$stop ".Data::Dumper::Dumper(@answers));
        warn $self->error_message;
        #die $self->error_message;
        #return undef;
    }elsif(scalar(@answers)==0){
        #return undef;
    }
    return $answers[0];
}

sub in_bp_range {
    my $self = shift;
    my ($chrom,$start,$stop) = split ",", shift;
    my ($lower_range,$upper_range) = split ",", shift;
    unless(defined($lower_range)){
        return;
    }
}

sub inside {
    my $self = shift;
    my ($chr,$lower_range,$upper_range,$ch,$st,$sp) = @_;
    unless($chr eq $ch){
        return undef;
    }
    if(($st>=$lower_range)&&($sp<=$upper_range)){
        return 1;
    }
    return undef;
}
