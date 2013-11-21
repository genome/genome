package Genome::Model::Tools::Somatic::CalculatePindelReadSupport;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use IO::File;

my %positions;
my %insertions;
my %deletions;

class Genome::Model::Tools::Somatic::CalculatePindelReadSupport {
    is => 'Command',
    has => [
        indels_all_sequences_bed_file =>{
            type => 'String',
            is_optional => 0,
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
            example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fasta'],
            doc => "reference sequence to use for reference assembly",
        },
        use_old_pindel => {
            type => 'Boolean',
            is_input => 1,
            is_optional => 1,
            default => 1,
            doc => 'Run on pindel 0.2 or 0.1',
        },
        germline_events => {
            type => 'Boolean',
            is_input => 1,
            is_optional => 1,
            default => 0,
            doc => 'Use this to calculate support for germline events, not tumor events',
        },
        _dbsnp_insertions => {
            type => 'String',
            is_optional => 1,
            example_values => ['/gscmnt/ams1102/info/info/dbsnp130_indels/insertions_start_stop_adjusted_dbsnp130'],
            doc => 'dbsnp insertion file',
        },
        _dbsnp_deletions => {
            type => 'String',
            is_optional => 1,
            example_values => ['/gscmnt/ams1102/info/info/dbsnp130_indels/deletions_adjusted_dbsnp130'],
            doc => 'dbsnp deletion file',
        },
        _output_filename => {
            calculate_from => [ 'indels_all_sequences_bed_file'], 
            calculate => q| $indels_all_sequences_bed_file.".read_support"|,
            is_output => 1,
        },
    ],
    has_param => [
         lsf_queue => {
             default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
         }, 
         lsf_resource => {
             default_value => "-M 16000000 -R 'select[type==LINUX64 && mem>16000] rusage[mem=16000]'",
         },
     ],
};

sub execute {
    my $self = shift;
    if(-s $self->_output_filename){
        $self->error_message("Output for calculate pindel read support at ".$self->_output_filename." already exists. Skipping.");
        return 1;
    }
    my $output = Genome::Sys->open_file_for_writing($self->_output_filename);
    my %indels;
    my %answers;

    my $ifh = Genome::Sys->open_file_for_reading($self->_dbsnp_insertions); #IO::File->new($self->_dbsnp_insertions);
    while (my $line = $ifh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
        next unless ($allele =~ m/-/);
        $allele = substr($allele, 2);
        $insertions{$chr}{$start}{$stop}{'allele'}=$allele;
        $insertions{$chr}{$start}{$stop}{'id'}=$id;
    }
    $ifh->close;
    my $dfh = Genome::Sys->open_file_for_reading($self->_dbsnp_deletions);#IO::File->new($self->_dbsnp_deletions);


    while (my $line = $dfh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
        next unless ($allele =~ m/-/);
        $allele = substr($allele, 2);
        $deletions{$chr}{$start}{$stop}{'allele'}=$allele;
        $deletions{$chr}{$start}{$stop}{'id'}=$id;
    }
    $dfh->close;

    my $fh = Genome::Sys->open_file_for_reading($self->indels_all_sequences_bed_file);
    while (<$fh>){
        my $line = $_;
        my ($chr,$start,$stop,$refvar) = split /\t/, $line;
        $positions{$chr}{$start}{$stop}=$refvar;
        my ($ref,$var) = split "/", $refvar;
        $indels{$chr}{$start} = $refvar;
    }

    #print $output "CHR\tSTART\tSTOP\tREF/VAR\tINDEL_SUPPORT\tREFERENCE_SUPPORT\t%+STRAND\tDBSNP_ID\n";
    for my $chr (sort(keys(%indels))){
        my %indels_by_chr = %{$indels{$chr}};
        $self->process_file($chr, \%indels_by_chr, $output);

    }
    $output->close;

}

sub process_file {
    my $self = shift;
    my $chr = shift;
    my $indels_by_chrom = shift;
    my $output = shift;
    my $dir = $self->pindel_output_directory;
    my $filename = $dir."/".$chr."/indels_all_sequences";
    my $pindel_output = Genome::Sys->open_file_for_reading($filename); #IO::File->new($filename);
    my $pindel_config = $dir."/".$chr."/pindel.config";
    my $pconf = Genome::Sys->open_file_for_reading($pindel_config);  #IO::File->new($pindel_config);
    my $tumor_bam = $pconf->getline;
    my $normal_bam;
    if($tumor_bam =~ m/normal/){
        ($normal_bam)= split /\s/, $tumor_bam;
        $tumor_bam = $pconf->getline;
        unless($tumor_bam =~ m/tumor/){
            die $self->error_message("Could not locate a tumor bam in the pindel config file at ".$pindel_config);
        }
    }
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
#charris speed hack
            unless(exists $indels_by_chrom->{$start}){
                next;
            }
            my @bed_line = $self->parse($call, $reference, $read);
            next unless scalar(@bed_line)>1;
            unless((@bed_line)&& scalar(@bed_line)==4){
                next;
            }
            my $type_and_size = $type."/".$size;
            $events{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'neg'}=$neg_strand;
            $events{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'pos'}=$pos_strand;
            $events{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'bed'}=join("\t",@bed_line);
            if($normal_support){
                $events{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'normal'}=$normal_support;
            }
        }
    }
    for $chrom (sort {$a cmp $b} (keys(%events))){
        for $pos (sort{$a <=> $b} (keys( %{$events{$chrom}}))){
            for my $type_and_size (sort(keys( %{$events{$chrom}{$pos}}))){
                unless(exists($events{$chrom}{$pos}{$type_and_size}{'normal'})&&(not $self->germline_events)){
                    my $pos_strand = $events{$chrom}{$pos}{$type_and_size}{'pos'};
                    my $neg_strand = $events{$chrom}{$pos}{$type_and_size}{'neg'};
                    my $pos_percent=0;
                    if($neg_strand==0){
                        $pos_percent = 1.0;
                    } else {
                        $pos_percent = sprintf("%.2f", $pos_strand / ($pos_strand + $neg_strand));
                    }
                    my $answer = "neg = ".$neg_strand."\tpos = ".$pos_strand." which gives % pos str = ".$pos_percent."\n";
                    my $reads = $pos_strand + $neg_strand;
                    my @stop = keys(%{$positions{$chrom}{$pos}});
                    #unless(scalar(@stop)==1){
                    #    
                    #    die "too many stop positions at ".$chrom." ".$pos."\n";
                    #}
                    my ($type,$size) = split /\//, $type_and_size;
                    #my $stop = ($type eq 'I') ? $pos+2 : $pos + $size;
                    my $stop = $pos;
                    #my @results = `samtools view $tumor_bam $chrom:$pos-$stop | grep -v "XT:A:M"`;
                    my @results = `samtools view $tumor_bam $chrom:$pos-$stop`;
                    my $tumor_read_support=0;
                    if($self->germline_events){
                        push @results, `samtools view $tumor_bam $chrom:$pos-$stop`;
                    }
                    my $read_support=0;
                    for my $result (@results){
                        #print $result;
                        chomp $result;
                        my @details = split /\t/, $result;
                        if($details[5] =~ m/[ID]/){
                            $tumor_read_support++;
                        }
                    }
                    @results = `samtools view $normal_bam $chrom:$pos-$stop`;
                    my $normal_read_support=0;
                    for my $result (@results){
                        #print $result;
                        chomp $result;
                        my @details = split /\t/, $result;
                        if($details[5] =~ m/[ID]/){
                            $normal_read_support++;
                        }
                    }
                    my $dbsnp_id = $self->dbsnp_lookup($events{$chrom}{$pos}{$type_and_size}{'bed'});
                    my $bed_output = $events{$chrom}{$pos}{$type_and_size}{'bed'}."\t".$reads."\t".$tumor_read_support."\t".$normal_read_support."\t".$pos_percent."\t$dbsnp_id\n";
                    print $output $bed_output;
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
    my ($chr, $start, $stop, $refvar) = split "\t", $bed_line;
    my ($ref,$var) = split "/",$refvar;
    if($ref eq "0") {
        if(exists($insertions{$chr}{$start}{$stop}{'allele'})) {
            if ($var eq $insertions{$chr}{$start}{$stop}{'allele'}) {
                $dbsnp_id=$insertions{$chr}{$start}{$stop}{'id'};
            }
        }
    }
    else {        
        if(exists($deletions{$chr}{$start}{$stop}{'allele'})) {
            if ($ref eq $deletions{$chr}{$start}{$stop}{'allele'}) {
                $dbsnp_id=$deletions{$chr}{$start}{$stop}{'id'};
            }
        } 
    }
    return $dbsnp_id;
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
        $stop = $stop - 1;
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
    my $refvar = "$ref/$var";
    return ($chr,$start,$stop,$refvar);
}

