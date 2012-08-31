package Genome::Model::Tools::Somatic::PindelToBed;

use warnings;
use strict;

use Genome;
use Workflow;

class Genome::Model::Tools::Somatic::PindelToBed {
    is => ['Genome::Model::Tools::Somatic'],
    has => [
        reference_fasta => {
            is => 'String',
            default=> Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa',
            doc => "The reference fasta file used to look up the reference sequence with samtools faidx. This is necessary because pindel will truncate long reference sequences.",
        },
        use_old_pindel => {
            type => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => 'Run on pindel 0.2 or 0.1',
        },
        include_normal => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Include events which have some or all normal support alongside events with only tumor support',
        },
        germline_only => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Include only events which have some normal support.',
        },
        include_bp_ranges => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Include pindels calculated bp_range for the location of the indel.',
        },
        output => {
            type => 'String',
            is_optional => 0,
            doc => 'The location of the output bed file.',
        },
        source => {
            type => 'String',
            is_optional => 0,
            doc => 'The pindel indels_all_sequences file or directory containing the run-pindel pipeline output, to be converted to bed.'
        },
    ],
    has_transient_optional => [
        _big_output_fh => {
            is => 'IO::File',
            doc => 'Filehandle for the output of large events',
        },
        _input_fh => {
            type =>'IO::File',
            doc => 'Filehandle for input',
        },
        _output_fh => {
            type =>'IO::File',
            doc =>'Filehandle for output',
        },

    ],
    has_param => [
        # Make workflow choose 64 bit blades, this is needed for samtools faidx
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1] -M 4000000',
        },
        lsf_queue => {
            default_value => 'long'
        }, 
    ],
};

sub help_brief {
    "Transforms a file from pindel output format to bed format";
}

sub help_synopsis {
    return <<"EOS"
gmt bed convert indel pindel-to-bam --input-file pindel.outfile --output-file pindel.adapted 
EOS
}

sub help_detail {                           
    return <<EOS 
Transforms a file from pindel output format to bed format
EOS
}

sub execute {
    my $self = shift;
    if($self->include_normal && $self->germline_only){
        $self->error_message("Cannot include normal AND germline only.");
        die $self->error_message;
    }
    # test architecture to make sure we can run (needed for samtools faidx)
    unless (`uname -a` =~ /x86_64/) { 
       $self->error_message("Must run on a 64 bit machine");
       die;
    }
    if(-d $self->source){
        $self->status_message("Using ".$self->source." as the pindel directory, proceeding to cat all chroms and run on that.");
        my @grob = glob($self->source."/*");
        my @dirs;
        my $temp_file = Genome::Sys->create_temp_file_path;
        for my $t (@grob){
            if(-d $t){
                if(-s $t."/indels_all_sequences"){
                    my $cmd = "cat ".$t."/indels_all_sequences >> ".$temp_file;
                    unless(Genome::Sys->shellcmd( cmd => $cmd)){
                        $self->error_message("Failed to cat files. command looked like: ".$cmd);
                        die $self->error_message;
                    }
                }
            }
        }
        $self->_input_fh(Genome::Sys->open_file_for_reading($temp_file));
    } else {
        $self->_input_fh(Genome::Sys->open_file_for_reading($self->source));
    }
    
    if($self->_big_output_fh) {
        return 1; #Already initialized
    }
    $self->_output_fh(Genome::Sys->open_file_for_writing($self->output));
    
    my $big_output = $self->output . ".big_deletions";
    my $big_output_fh;
    eval {
        $big_output_fh = Genome::Sys->open_file_for_writing($big_output);
        $self->_big_output_fh($big_output_fh);
    };
    
    if($@) {
        $self->error_message('Failed to open file. ' . $@);
        $self->close_filehandles;
        return;
    }
    my $input_fh = $self->_input_fh;
    my %events;
    my %ranges;
    my ($chrom,$pos,$size,$type);
    my $type_and_size;
    while(my $line = $input_fh->getline){
        my $normal_support=0;
        my $read = 0;
        if($line =~ m/^#+$/){
            my $call = $input_fh->getline;
            my $reference = $input_fh->getline;
            my @call_fields = split /\s/, $call;
            $type = $call_fields[1];
            $size = $call_fields[2];   #12
            my $mod = ($call =~ m/BP_range/) ? 2: -1;
            #my $support = ($type eq "I") ? $call_fields[10+$mod] : $call_fields[12+$mod];
            ###charris patch for use old pindel  
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
            my $lower_range = ($type eq "I") ? $call_fields[7+$mod] : $call_fields[9+$mod];
            my $upper_range = ($type eq "I") ? $call_fields[8+$mod] : $call_fields[10+$mod];
            ### end charris patch
            for (1..$support){
                $line = $input_fh->getline;
                if($line =~ m/normal/) {
                    $normal_support++;
                }
                $read=$line;
            }
            
            my @bed_line = $self->parse($call, $reference, $read);
            unless((@bed_line)&& scalar(@bed_line)==5){
                next;
            }
            $type_and_size = $type."/".$size;
            $events{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'bed'}=join(",",@bed_line);
            if($self->include_bp_ranges){
                $ranges{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'lower'}=$lower_range;
                $ranges{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'upper'}=$upper_range;
            }
            if($normal_support){
                $events{$bed_line[0]}{$bed_line[1]}{$type_and_size}{'normal'}=$normal_support;
            }
        }
    }

    for $chrom (sort {$a cmp $b} (keys(%events))){
        for $pos (sort{$a <=> $b} (keys( %{$events{$chrom}}))){
            for $type_and_size (sort(keys( %{$events{$chrom}{$pos}}))){
                if((not (exists($events{$chrom}{$pos}{$type_and_size}{'normal'}))||$self->include_normal)&& not $self->germline_only){
                    my @bed = split ",", $events{$chrom}{$pos}{$type_and_size}{'bed'};
                    if($self->include_bp_ranges){
                        push @bed , $ranges{$chrom}{$pos}{$type_and_size}{'lower'};
                        push @bed , $ranges{$chrom}{$pos}{$type_and_size}{'upper'};
                    }
                    my $bed = join(",",@bed);
                    $self->write_bed_line(split ",", $bed);
                }
                if($self->germline_only){
                    if(exists($events{$chrom}{$pos}{$type_and_size}{'normal'})){
                        my @bed = split ",", $events{$chrom}{$pos}{$type_and_size}{'bed'};
                        if($self->include_bp_ranges){
                            push @bed , $ranges{$chrom}{$pos}{$type_and_size}{'lower'};
                            push @bed , $ranges{$chrom}{$pos}{$type_and_size}{'upper'};
                        }
                        my $bed = join(",",@bed);
                        $self->write_bed_line(split ",", $bed);
                    }
                }
            }
        }
    }

    $big_output_fh = $self->_big_output_fh;
    close($big_output_fh) if $big_output_fh;

    return 1;
}

sub parse {
    my $self=shift;
    #my $reference_fasta = $self->refseq;
    my ($call, $reference, $first_read) = @_;
    #parse out call bullshit
    chomp $call;
    my @call_fields = split /\s+/, $call;
    my $type = $call_fields[1];
    my $size = $call_fields[2];
    ####use old pindel patch######
    my ($chr, $start, $stop);
    if($self->use_old_pindel){
        $chr = ($type eq "I") ? $call_fields[4] : $call_fields[6];
        $start= ($type eq "I") ? $call_fields[6] : $call_fields[8];
        $stop = ($type eq "I") ? $call_fields[7] : $call_fields[9];
    } else {
        $chr = $call_fields[6];
        $start= $call_fields[8];
        $stop = $call_fields[9];
    }
    ####end charris use old pindel patch
    my $support = $call_fields[-1];
    my ($ref, $var);
    if($type =~ m/D/) {
        $var =0;
        ###Make pindels coordinates(which seem to be last undeleted base and first undeleted base) 
        ###conform to our annotators requirements
        $stop = $stop -1;
        ###also deletions which don't contain their full sequence should be dumped to separate file
        my $allele_string;
        my $start_for_faidx = $start+1; 
        my $sam_default = Genome::Model::Tools::Sam->path_for_samtools_version;
        my $faidx_cmd = "$sam_default faidx " . $self->reference_fasta . " $chr:$start_for_faidx-$stop"; 
        my @faidx_return= `$faidx_cmd`;
        shift(@faidx_return);
        chomp @faidx_return;
        $allele_string = join("",@faidx_return);

        $ref = $allele_string;
    }
    elsif($type =~ m/I/) {
        #misunderstanding of bed format
        #0 based numbers teh gaps so an insertion of any number of bases between base 10 and 11 in 1base
        #is 10 10 in bed format
        #$start = $start - 1;
        $ref=0;
        my ($letters_until_space) =   ($reference =~ m/^([ACGTN]+) /);
        my $offset_into_first_read = length($letters_until_space);
        $var = substr($first_read, $offset_into_first_read, $size);
        $stop = $stop - 1;
    }
    if($size >= 100) {
        my $big_fh = $self->_big_output_fh;
        $big_fh->print("$chr\t$start\t$stop\t$size\t$support\n");
        return undef;
    }
    return ($chr,$start,$stop,$ref,$var);
}
sub write_bed_line {
    my $self = shift;
    #start is zero-based index of first base in the event.
    #stop is zero-based index of first base *after* the event (e.g. for a SNV these will differ by one)
   
     #my ($chromosome, $start, $stop, $reference, $variant) = @_;
    my @values = @_;
    
    my $output_fh = $self->_output_fh;
    
    #my $name = join('/', $reference, $variant);
    my $name = join('/', $values[3], $values[4]);
    my @columns;
    for my $index (0..scalar(@values)){
        if($index==3){
            push @columns, $name;
        } elsif ($index==4){
        } else {
            if(defined($values[$index])){
                push @columns, $values[$index];
            }
        }
    }
    #print $output_fh join("\t", $chromosome, $start, $stop, $name), "\n";
    print $output_fh join("\t", @columns), "\n";
    
    return 1;
}
