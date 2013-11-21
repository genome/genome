package Genome::Model::Tools::Somatic::AssembleIndelBed;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use IO::File;
use Cwd qw( abs_path );
my $SAM_DEFAULT = Genome::Model::Tools::Sam->default_samtools_version;



class Genome::Model::Tools::Somatic::AssembleIndelBed {
    is => 'Command',
    has => [
        indel_file =>{
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Indel sites to assemble in annotator input format',
        },
        bam_file =>{
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'File from which to retrieve reads',
        },
        buffer_size =>{
            type => 'Integer',
            is_optional => 1,
            default => 100,
            doc => 'Size in bp around start and end of the indel to include for reads for assembly',
        },
        data_directory =>   {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => "Location to dump individual chr results etc",
        },
        refseq =>{
            type => 'String',
            is_optional => 1,
            example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fasta'],
            doc => "reference sequence to use for reference assembly",
        },
        assembly_indel_list =>{
            type => 'String',
            is_input=>1,
            is_output=>1,
            doc => "List of assembly results",
        },
        sam_version => {
            is  => 'String',
            doc => "samtools version to be used, default is $SAM_DEFAULT",
            default_value => $SAM_DEFAULT,
            is_optional => 1,
        },
        lsf_resource => {
            is_param => 1,
            default_value => 'rusage[mem=2000] select[type==LINUX64 & mem > 2000] span[hosts=1]',
        },
        lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        } 
    ]
};




sub dir_for_chrom {
    my $self=shift;
    my $chr = shift;
    if($chr) {
        return $self->data_directory . "/$chr";
    }
    else {
        return $self->data_directory . "/";
    }
}



sub execute {
    my $self=shift;
    my %DONE;
    my $DONE = \%DONE;

    #test architecture to make sure we can run samtools
    #copied from G::M::T::Maq""Align.t 
    
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        return;
    }
    my $bam_file = $self->bam_file;
    unless(-s $bam_file) {
        $self->error_message("$bam_file does not exist or has zero size");
        return;
    }
    my $dir = $self->data_directory;
    unless(-d $dir) {
        $self->error_message("$dir is not a directory");
        return;
    }

    my $refseq = $self->refseq;
    unless(-s $refseq) {
        $self->error_message("$refseq does not exist or has zero size");
        return;
    }


    my $output_file = $self->assembly_indel_list;
    my $output_fh = IO::File->new($output_file, ">");
    unless($output_fh) {
        $self->error_message("Couldn't open $output_file: $!"); 
        return;
    }


    my $indel_file = $self->indel_file;
    my $fh = IO::File->new($indel_file, "r");
    unless($fh) {
        $self->error_message("Couldn't open $indel_file: $!"); 
        return;
    }

    my $sam_pathname = Genome::Model::Tools::Sam->path_for_samtools_version($self->sam_version);

    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $stupid_ref_var_string_gabe_makes_us_use) = split /\t/, $line;
        my ($ref,$var) = split /\//, $stupid_ref_var_string_gabe_makes_us_use;
        my $dir = $self->dir_for_chrom($chr);
        unless(-e $dir) {
            `mkdir -p $dir`;
        }
        next if $chr eq 'chromosome_name';
        if($ref && $ref ne '-' && $var && $var ne '-') {
            $self->error_message("No indel found at $fh->input_line_number. Skipping...");
            next;
        }

        my $region_start = $start - $self->buffer_size;
        my $region_stop = $stop + $self->buffer_size;
        my @reads = `$sam_pathname view $bam_file $chr:$region_start-$region_stop | cut -f1,10`;
        if(@reads) {
            my $prefix = "$dir/${chr}_${start}_${stop}";
            my $read_file = "$prefix.reads.fa";
            my $fa_fh = IO::File->new( $read_file,"w");
            unless($fa_fh) {
                $self->error_message("Unable to open $read_file for writing");
                return;
            }
            foreach my $read (@reads) {
                chomp $read;
                next unless $read;
                my ($readname, $sequence) = split /\t/, $read;
                print $fa_fh ">$readname\n$sequence\n";
            }
            $fa_fh->close;

            #make reference fasta
            my $ref_file = "$prefix.ref.fa";
            my $ref_fh = IO::File->new($ref_file,"w");
            unless($ref_fh) {
                $self->error_message("Unable to open $ref_file for writing");
                return;
            }
            my @contig = `$sam_pathname faidx $refseq $chr:$region_start-$region_stop`;
            if(@contig) {
                print $ref_fh @contig;
            }
            $ref_fh->close;

            print `/gsc/scripts/pkg/bio/tigra/installed/local_var_asm_wrapper.sh $read_file`; #assemble the reads
            `cross_match $read_file.contigs.fa $ref_file -bandwidth 20 -minmatch 20 -minscore 25 -penalty -4 -discrep_lists -tags -gap_init -4 -gap_ext -1 > $prefix.stat`;
            my $cmd = "gmt parse crossmatch --chr-pos ${chr}_${region_start} --crossmatch=$prefix.stat --min-indel-size=1";
            print "$cmd\n";
            my ($stupid_header1, $stupid_header2, $result) = `$cmd`;
            if(defined $result && $result =~ /\S+/) {
                my $output = $self->thing($result, $DONE);
                $output_fh->print($output);
                print $result . "\n";
            }
            print STDERR "########################\n";
        }
    }

    return 1;
}


1;

sub thing {
    my $self = shift;
    my $somatic_event = shift;
    my $DONE = shift;
    my @fields = split /\t/, $somatic_event;
    my $chr = $fields[0];
    if($chr =~m /:/) {
        my ($new_chr, ) = split ":", $chr;
        $chr = $new_chr;
        $fields[0]=$chr;
    }
    my $pos = $fields[1];
    my $size = $fields[7];
    my $type = $fields[8];
    my $start;
    my $stop;
    if ($type =~ m/DEL/) {
        $start=$pos; 
        $stop= $pos + $size;
    }
    elsif($type =~ m/INS/ ) {
        $start=$pos;
        $stop=$pos;#+1;
    }
    my ($ref, $var) = $self->generate_alleles($DONE, @fields);
    my $output_line = join("\t",$chr,$start,$stop,$ref,$var);

    #if line has already been printed, go to next event
    if ($DONE->{$output_line}) {
        next;
    }
    #else, record printed event
    else {
        $DONE->{$output_line} = 1;
        return "$output_line\n";
    }
}

#returns a list of the reference allele and the variant allele in that order 
sub generate_alleles {
    my ($self, $DONE, @assembled_event_fields) = @_;
    my ($chr,$pos,$contig_pos,$size,$type,$contig_name,$reference_name) = @assembled_event_fields[0,1,5,7,8,9,11];

    if($type =~ /INS/) {
        my ($original_position) = $reference_name =~ m/_(\d+)/;
        $original_position += 100; #regenerate the original position
            my $glob_pattern = $self->data_directory . "/$chr/$chr" . "_" . $original_position . "_*.reads.fa.contigs.fa";
        my ($contig_filename,@others) = glob($glob_pattern);

        #if there is more than one glob_pattern found
        if(@others) {

            #if any of the filenames have already been used, find one which hasn't been used
            if (exists($DONE->{'glob'}{$chr}{$original_position})) {
                #roll through the filenames
                for my $filename ($contig_filename, @others) {
                    #if it's been used, go to the next one
                    if (scalar grep { m/^$filename$/ } @{$DONE->{'glob'}{$chr}{$original_position}}) {
                        next;
                    }
                    #if it hasn't been used, push it on the array of used filenames and use it
                    else {
                        push @{$DONE->{'glob'}{$chr}{$original_position}},$filename;
                        $contig_filename = $filename;
                        last;
                    }
                }
            }
            #if the hash key doesn't exist, none of the filenames have been used, so record this first one
            else {
                $DONE->{'glob'}{$chr}{$original_position} = [$contig_filename];
            }               

            $self->error_message("Hey, there are multiple contig files this thingy could belong too. Choosing the first filename and logging its use so that it is avoided if glob'ed again.");
            $self->error_message("Glob Pattern: $glob_pattern; Used file $contig_filename");
            #$self->error_message("Choosing first file in the hope that the calls ended up identical.");
        }

        my $contig_fh = IO::File->new($contig_filename, "r");
        unless($contig_fh) {
            $self->error_message("Unable to open $contig_filename");
            $self->error_message("Glob Pattern: $glob_pattern");
            $self->error_message("contig data for event not found: @assembled_event_fields");
            die;
        }
        while(my $line = $contig_fh->getline) {
            chomp $line;
            next if($line !~ /^>$contig_name/);
            my $sequence = '';
            while($line = $contig_fh->getline) {
                chomp $line;
                $sequence .= $line if $line !~ /^>Contig/;;
                if($line =~ /^>Contig/ || $contig_fh->eof) {
                    #inserted sequence first nucleotide is the base reported in the assembly results
                    return (0, substr $sequence,$contig_pos-1,$size);
                }
            }
        }
        $self->error_message("Unable to find _${contig_name}_ in $contig_filename");
        return;
    }
    else {
        #fetch reference sequence
        #for this the first deleted base should be the first base listed in the file
        my $ref_seq = $self->refseq;
        my $end = $pos + $size - 1;
        my ($header,$sequence) = `samtools faidx $ref_seq ${chr}:$pos-$end`;
        chomp $sequence;
        unless($sequence) {
            $self->error_message("Unable to retrieve sequence from reference");
            return;
        }
        return ($sequence, 0);
    }

}



sub help_brief {
    "Scans a snp file and finds adjacent sites. Then identifies if these are DNPs and annotates appropriately."
}


#copied directly from Ken
sub ComputeTigraN50{
    my ($self,$contigfile)=@_;
    my @sizes;
    my $totalsize=0;
    open(CF,"<$contigfile") || die "unable to open $contigfile\n";
    while(<CF>){
        chomp;
        next unless(/^\>/);
        next if(/\*/);
        my ($id,$size,$depth,$ioc,@extra)=split /\s+/;
        next if($size<=50 && $depth<3 && $ioc=~/0/);
        push @sizes,$size;
        $totalsize+=$size;
    }
    close(CF);
    my $cum_size=0;
    my $halfnuc_size=$totalsize/2;
    my $N50size;
    @sizes=sort {$a<=>$b} @sizes;
    while($cum_size<$halfnuc_size){
        $N50size=pop @sizes;
        $cum_size+=$N50size;
    }
    return $N50size;
}
