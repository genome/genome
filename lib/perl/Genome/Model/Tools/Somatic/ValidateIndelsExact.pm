package Genome::Model::Tools::Somatic::ValidateIndelsExact;

use strict;
use warnings;

use Genome;
use IO::File;

#my %positions;
#my %insertions;
#my %deletions;
my %ranges;

my $use_ranges=0;

class Genome::Model::Tools::Somatic::ValidateIndelsExact {
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
        output_file => {
            type => 'String',
            is_optional => 0,
            is_output => 1,
            doc => "The match results will be dumped here",
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
        intersection_range => {
            type => 'Integer',
            is_input => 1,
            is_optional => 0,
            default => -1,
            doc => "Use this to specify the range within to look for fuzzy concordance. Default is -1, which performs no fuzzy concordance. Use -999 to find fuzzy hits within a distance set by each indel length.",
        },
        reference_build	=> {
            is => 'Text',
            doc => "36 or 37",
            is_optional => 1,
            is_input => 1,
            example_values => [36],
        },
    ],
    has_param => [
         lsf_queue => {
             default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
         }, 
         lsf_resource => {
             default_value => "-M 6000000 -R 'select[type==LINUX64 && mem>16000] rusage[mem=16000]'",
         },
     ],
};

sub help_brief {
    "Compares the exact location of indels, but first left-shifts to normalize for indels in repetitive regions"
}

sub help_detail {
    <<'HELP';
    Compares the exact location of indels, but first left-shifts to normalize for indels in repetitive regions
HELP
}

#Internally, this tool uses bed format to represent chr start stop of events.
sub load_bed_file {
    my $self = shift;
    my $file = shift;
    my %data;

    my $fh = Genome::Sys->open_file_for_reading($file);
    $DB::single=1;
    while(my $line = $fh->getline){
        chomp $line;
        my ($chr,$start,$stop,$ref,$var,$lower_range,$upper_range) = split "\t",$line;

        if($ref =~ m/\//){
            $upper_range = $lower_range;
            $lower_range = $var;
            ($ref,$var) = split "/", $ref;
        }
        if($ref eq '*'){
                $ref = '0';
        }
        if($var eq '*'){
                $var = '0';
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

if ($start == 56215246) { print "3-$chr,$start,$stop,$ref,$var\n";}
if ($start == 56215247) { print "3b-$chr,$start,$stop,$ref,$var\n";}
        if($ref =~ m/\//){
            $upper_range = $lower_range;
            $lower_range = $var;
            ($ref,$var) = split "/", $ref;
        }

        if($ref eq '*'){
                $ref = '0';
        }
        if($var eq '*'){
                $var = '0';
        }
        
        if($ref eq 'O'){
                $ref = '0';
                print "wtf? O??I'm changing this ref to 0.\n";
        }
        if($var eq 'O'){
                $var = '0';
                print "wtf? O??I'm changing this var to 0.\n";
        }

        if(defined($lower_range)){
            $use_ranges = 1;
        }
        if($var =~ /[ACTGN]/i && $ref !~ /[ACTGN]/i) {
            $start++;
            $lower_range++ if $lower_range;
        }
        elsif($ref =~ /[ACTGN]/i && $var !~ /[ACTGN]/i) {
            $start--;
            $lower_range-- if $lower_range;
        }
        else {
            $self->error_message("This line was not an indel: $line");
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

sub leftshift_indels {
    my $self = shift;
    my $file = shift;
    my $filetype = shift;

    my $reference_build = $self->reference_build;
    my $reference_build_fasta;
    if ($reference_build =~ m/36/) {
            my $reference_build_fasta_object= Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
            $reference_build_fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif ($reference_build =~ m/37/) {
            my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
            $reference_build_fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    else {
           die "Please specify either build 36 or 37";
    }

    ## Build temp file of leftshifted indels##
    my ($tfh2,$temp_path2) = Genome::Sys->create_temp_file;
    unless($tfh2) {
            $self->error_message("Unable to create temporary file $!");
            die;
    }
    $temp_path2 =~ s/\:/\\\:/g;

    my $normalize_cmd;
    if ($filetype eq 'bed') {
        ## Build temp file of annotation format##
        my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
        unless($tfh) {
                $self->error_message("Unable to create temporary file $!");
                die;
        }
        $temp_path =~ s/\:/\\\:/g;

        my $fh = Genome::Sys->open_file_for_reading($file);
        while(my $line = $fh->getline){
                chomp $line;
                my ($chr,$start,$stop,$ref,$var,$lower_range,$upper_range) = split "\t",$line;
                if($ref =~ m/\//){
                    $upper_range = $lower_range;
                    $lower_range = $var;
                    ($ref,$var) = split "/", $ref;
                }
if ($start == 56215246) { print "1-$chr\t$start\t$stop\t$ref\t$var\n";}
if ($start == 56215247) { print "1b-$chr\t$start\t$stop\t$ref\t$var\n";}

                if($start == $stop) { #insertion
                        $stop = $stop+1;
                        print $tfh "$chr\t$start\t$stop\t$ref\t$var\n";
                }
                else {
                        $start = $start+1;
                        print $tfh "$chr\t$start\t$stop\t$ref\t$var\n";
                } 
if ($start == 56215246) { print "2-$chr\t$start\t$stop\t$ref\t$var\n";}
if ($start == 56215247) { print "2b-$chr\t$start\t$stop\t$ref\t$var\n";}

        }
        $normalize_cmd = "gmt analysis normalize-indel-location --annotation-input-file $temp_path --reference-fasta $reference_build_fasta --output-file $temp_path2";
        close($tfh);
    }
    else {
        $normalize_cmd = "gmt analysis normalize-indel-location --annotation-input-file $file --reference-fasta $reference_build_fasta --output-file $temp_path2";
    }
print "$normalize_cmd";
    system($normalize_cmd);
    close($tfh2);
    return $temp_path2;

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

    my $fh;
    if ($self->output_file) {
        $fh = Genome::Sys->open_file_for_writing($self->output_file);
    }

    print "Left-shifting valid indel data\n";
    my $valid_data_leftshift_file = (defined($self->valid_indel_bed_file))? $self->leftshift_indels($self->valid_indel_bed_file,"bed"): $self->leftshift_indels($self->valid_indel_annotation_file,"annotation");
    print "Loading valid indel data\n";
    my %valid_data = %{$self->load_annotation_file($valid_data_leftshift_file)};

    print "Left-shifting called indel data\n";
    my $called_data_leftshift_file = (defined($self->called_indel_bed_file))? $self->leftshift_indels($self->called_indel_bed_file,"bed"): $self->leftshift_indels($self->called_indel_annotation_file,"annotation");
    print "Loading called indel data\n";
    my %called_data = %{$self->load_annotation_file($called_data_leftshift_file)};
    print "Loaded all data, running comparisons\n";

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
                    }
                    elsif(uc($valid_data{$chr}{$start}{$stop}) eq uc($called_data{$chr}{$start}{$stop})){
                        $perfect_match{$chr}{$start}{$stop} = $called_data{$chr}{$start}{$stop};
                        $touched_valids{$chr}{$start}{$stop} = 1;
                        $exact_hits++;
                    }
                    else{
                        $allele_mismatch{$chr}{$start}{$stop} = "Called:$called_data{$chr}{$start}{$stop}\tValid:$valid_data{$chr}{$start}{$stop}";
print "Called:$called_data{$chr}{$start}{$stop}\tValid:$valid_data{$chr}{$start}{$stop}\n";
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

    if ($self->output_file) {
        print $fh "Total Calls\tTotal Validation Events\tTotal Hits\tExact Hits\tExact Position But Different Alleles\tFuzzy Hits\tComplete Misses\tFuzzy Indels Touched\n";
        print $fh $called_total."\t";
        print $fh $valid_total."\t";
        my $total_hits = ($exact_hits+$within_range+$allele_mismatch_hits);
        print $fh $total_hits."\t";
        print $fh $exact_hits."\t";
        print $fh $allele_mismatch_hits."\t";
        print $fh (($fuzzy_hits_exact_allele+$fuzzy_hits)) unless $self->intersection_range < 0;
        print $fh "\t";
        print $fh $total_fail unless $self->intersection_range < 0;
        print $fh "\t";
        print $fh $touched."\n";
        print $fh "Exact Hits\n";

        for my $chr (sort keys(%perfect_match)){
            for my $start (sort {$a <=> $b} keys(%{$perfect_match{$chr}})){
                for my $stop (sort {$a <=> $b} keys(%{$perfect_match{$chr}{$start}})){
                    my $stuff = join("\t",($chr,$start,$stop,$perfect_match{$chr}{$start}{$stop}))."\n";
                    print $fh $stuff;
                }
            }
        }
        print $fh "Position Hits (Allele Mismatch)\n";
        for my $chr (sort keys(%allele_mismatch)){
            for my $start (sort {$a <=> $b} keys(%{$allele_mismatch{$chr}})){
                for my $stop (sort {$a <=> $b} keys(%{$allele_mismatch{$chr}{$start}})){
                    my $stuff = join("\t",($chr,$start,$stop,$allele_mismatch{$chr}{$start}{$stop}))."\n";
                    print $fh $stuff;
                }
            }
        }
        $fh->close;
    }


    #print Data::Dumper::Dumper(%fuzzy_candidates);
    $self->dump_to_file(\%fuzzy_candidates,$self->fuzzy_output_file) if $self->fuzzy_output_file;
    $self->dump_to_file(\%total_fail,$self->missed_output_file) if $self->missed_output_file;



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
    if ($range == -999) {
        if (($stop - $start) == 0) {
            $range = length($var);
            unless (defined($var)) {
                print "$chr,$start,$stop,$ref,$var\n";
            }
        }
        else {
            $range = ($stop - $start);
        }
    }
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
