package Genome::Model::Tools::Somatic::IntersectAssembledIndels;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use IO::File;
use Cwd qw( abs_path );
my $SAM_DEFAULT = Genome::Model::Tools::Sam->default_samtools_version;



class Genome::Model::Tools::Somatic::IntersectAssembledIndels {
    is => 'Command',
    has => [
    tumor_indel_file =>
    {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'Result of Tumor Assembly',
    },
    normal_indel_file =>
    {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'Result of Normal Assembly',
    },
#    input_to_assembly =>
#    {
#        type => 'String',
#        is_optional => 1,
#        is_input => 1,
#        doc => 'Original (Samtools|Pindel) Calls from Tumor',
#    },
    somatic_output_list =>
    {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        is_output => 1,
        doc => 'Output file',
    },
    germline_output_list => 
    {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'Output File',
    },
    tumor_assembly_data_directory =>
    {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'data directory of tumor reads and contigs',
    },
    reference =>
    {
        type => 'String',
        is_optional => 1,
        is_input => 1,
        example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa'],
        doc => 'the reference sequence the reads were aligned to. Used to generate deletion alleles',
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

sub execute {
    my $self = shift;
    my @somatic_events;
    my @germline_events;
    unless( -s $self->tumor_indel_file && -s $self->normal_indel_file) {
        $self->error_message("Tumor or normal input not found. Aborting\n");
        return;
    }
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        return;
    }

    
    my $tumorfh = IO::File->new($self->tumor_indel_file);
    my $normalfh = IO::File->new($self->normal_indel_file);
    my %normal_hash;
    while (my $line = $normalfh->getline) {
        chomp $line;
        my ($chr,$pos) = split "\t", $line;
        if ($chr =~ m/:/) {
            
            $self->warning_message("Janky tumor output: $line\n");
            my ($new_chr, )  = split ":", $chr;
            $chr = $new_chr;
        }
        $normal_hash{$chr}{$pos}=1;
    }
    while (my $line = $tumorfh->getline) {
        chomp $line;
        my ($chr,$pos) = split "\t", $line;
        if ($chr =~ m/:/) {
            $self->warning_message("Janky normal output: $line\n");
            my ($new_chr, )  = split ":", $chr;
            $chr=$new_chr;
        }
        if(exists($normal_hash{$chr}{$pos})) {
            push @germline_events, $line;
        }
        else {
            push @somatic_events, $line;
        }
    }

    #a hash to record what has been collated and ouput
    my %DONE;
    my $DONE = \%DONE;
    
    $self->collate_and_output_daterz($self->somatic_output_list, $DONE, @somatic_events);
    $self->collate_and_output_daterz($self->germline_output_list, $DONE, @germline_events);

}

sub collate_and_output_daterz {
    my $self=shift;
    my $output_filename= shift;
    my $DONE = shift;
    my @somatic_events = @_;


    my $output_fh = IO::File->new($output_filename, ">");

    # unless(-s $self->input_to_assembly) {
    #    $self->warning_message("No original calls supplied or file not found. Outputting assembly output"); 
    #    $output_fh->print(join("\n", @somatic_events));
    # return;
#    }
    #load original calls
    #my $allelefh = IO::File->new($self->input_to_assembly);
    #my @original_calls = $allelefh->getlines;


    my $original_call;
    my $original_pos;
    for my $somatic_event (@somatic_events) {
        chomp $somatic_event;
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
            $stop=$pos+1;
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
            $output_fh->print("$output_line\n");
        }
    }
}

#returns a list of the reference allele and the variant allele in that order 
sub generate_alleles {
    my ($self, $DONE, @assembled_event_fields) = @_;
    my ($chr,$pos,$contig_pos,$size,$type,$contig_name,$reference_name) = @assembled_event_fields[0,1,5,7,8,9,11];

    if($type =~ /INS/) {
        my ($original_position) = $reference_name =~ m/_(\d+)/;
        $original_position += 100; #regenerate the original position
        my $glob_pattern = $self->tumor_assembly_data_directory . "/$chr/$chr" . "_" . $original_position . "_*.reads.fa.contigs.fa";
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
        my $ref_seq = $self->reference;
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

