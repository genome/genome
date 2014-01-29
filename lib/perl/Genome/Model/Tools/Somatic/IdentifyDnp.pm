package Genome::Model::Tools::Somatic::IdentifyDnp;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;
use IO::File;
use POSIX;

class Genome::Model::Tools::Somatic::IdentifyDnp {
    is => 'Command',
    has => [
        annotation_input_file => {
            type => 'String',
            is_optional => 0,
            doc => 'List of sites in annotation input format to look for DNPs. This must be sorted by chromosome and coordinate.',
            default => '',
        },
        proportion => {
            type => 'Float',
            is_optional => 1,
            default => 0.1,
            doc => 'Proportion of reads supporting the DNP required for the site to be considered a DNP',
        },
        bam_file => {
            type => 'String',
            is_optional => 0,
            doc => 'File from which to retrieve reads. Must be indexed.',
        },
    ],
    has_optional => [
        output_file => {
            is => 'Text',
            doc => 'Output file; STDOUT if not specified.',
        },
    ],
};
#Code should operate as follows
#Scan the snp file
#Upon finding adjacent snps
#retrieve reads for the two snps
#Check and see if they are in cis ie in the same read together

sub execute {
    my $self = shift;

    #TODO Add checks on the files and architecture
    unless (POSIX::uname =~ /64/) {
        $self->error_message("This script requires a 64-bit system to run samtools");
        return;
    }

    unless(-e $self->bam_file && !-z $self->bam_file) {
        $self->error_message($self->bam_file . " does not exist or is of zero size");
        return;
    }

    my $fh = Genome::Sys->open_file_for_reading($self->annotation_input_file);

    my $output_fh;
    if ($self->output_file) {
        $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    }
    else {
        $output_fh = *STDOUT;
    }

    my $last_pos = undef;
    my $last_chr = undef;
    my $last_var = undef;
    my $last_ref = undef;
    my $last_cns = undef;
    my @lines = (); #line buffer to store last few lines

    #the following logic assumes that you have a single position per line
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $pos, $stop, $ref, $cns, $type, @rest) = split /\t/, $line;

        my @variants = Genome::Info::IUB::variant_alleles_for_iub($ref, $cns);
        unless(@variants == 1) {
            warn "Multiple variant alleles not yet supported. Only taking first";
        }

        if($last_chr && $last_pos) {
            if($chr eq $last_chr) {
                #TODO when the annotator supports DNPs that are non-adjacent. This can be upped to 2 and the output line changed appropriately
                #FIXME this doesn't support tri-nucleotide polymorphisms and up.
                if(($pos - $last_pos) == 1) {
                    #potential DNP
                    $self->debug_message("Potential DNP found at $chr:$pos-$last_pos\n");

                    if($self->is_dnp($chr, $last_pos, $last_var, $pos, $variants[0])) {
                        my $ref_col = $last_ref . $ref;
                        my $var_col = $last_var . $variants[0];
                        print $output_fh join("\t", $chr, $last_pos, $pos, $ref_col, $var_col, 'DNP', @rest) . "\n";
                        @lines = ();    #erase the DNP line
                        $line = '';
                    }
                }
            }
        }
        print $output_fh shift @lines,"\n" if @lines;
        push @lines, $line if $line;

        $last_chr = $chr;
        $last_pos = $pos;
        $last_var = $variants[0];
        $last_ref = $ref;
        $last_cns = $cns;
    }
    print $output_fh shift @lines,"\n" if @lines;

    return 1;
}


1;

sub help_brief {
    "Scans an annotation file, finds adjacent sites and then identifies if these are DNPs."
}

sub help_detail {
    <<'HELP';
This is a simple script which operated by identifying adjacent sites in a SORTED annotation file (really, just chr, start, stop, ref, cns, and type are need in that order) and then tries to determine if the alleles are linked in the same reads. If they are then the original lines are not printed and a new line, wrapping them into a single DNP event, is printed. This line is NOT reannotated, and will need to be annotated separately after running. Non-DNP sites are simply printed as is. Currently, only DNPs are examined and tri-nucleotide polymorphisms on up will not be properly identified. The input file was intended to be the pre-annotation "adapted file" of SNVs from the somatic pipeline file and, therefore, should NOT contain indels.
HELP
}


#I can't say why this is just a wrapper for the other program
sub is_dnp {
    my ($self, $chr, $pos1, $base1, $pos2, $base2) = @_;
    return $self->_determine_dnp_from_bam_reads($self->bam_file,$chr,$pos1,$base1,$pos2,$base2);
}

#This grabs the reads overlapping the positions
#and checks to see if they contain both potential DNP bases
sub _determine_dnp_from_bam_reads {
    my ($self, $alignment_file, $chr, $pos1, $base1, $pos2, $base2) = @_;

    unless(open(SAMTOOLS, "samtools view $alignment_file $chr:$pos1-$pos2 |")) {
        $self->error_message("Unable to open pipe to samtools view");
        return;
    }

    my ($reads, $reads_supporting_dnp) = (0,0);
    while( <SAMTOOLS> ) {
        chomp;
        my ($qname, $flag, $rname, $pos_read, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $RG, $MF, @rest_of_fields) = split /\t/;
        next if($mapq == 0); #only count q1 and above

        my $offset1 = $self->_calculate_offset($pos1, $pos_read, $cigar);
        next unless defined $offset1; #skip deletions
        my $offset2 = $self->_calculate_offset($pos2, $pos_read, $cigar);
        next unless defined $offset2; #skip deletions
        $reads++;

        if(uc(substr($seq,$offset1,1)) eq uc($base1) && uc(substr($seq,$offset2,1)) eq uc($base2)) {
            $reads_supporting_dnp++;
        }


    }

    unless(close(SAMTOOLS)) {
        $self->error_message("Error running samtools");
        return -1;
    }

    if($reads > 0 and $reads_supporting_dnp/$reads > $self->proportion) {
        return 1;
    }
    else {
        return 0;
    }
}

#this calculates the offset of a position into a seqeunce string based on the CIGAR string specifying the alignment
#these are some tests used to test if I got this right.

#use Test::Simple tests => 2;
#my $fake_read_seq = "ACTATCG";
#my $fake_read_pos = 228;
#my $fake_cigar = "3M1D4M";
##
#my $offset = calculate_offset(231,$fake_read_pos, $fake_cigar);
#ok(!defined($offset));
#$offset = calculate_offset(232,$fake_read_pos, $fake_cigar);
#ok(substr($fake_read_seq,$offset,1) eq "A");
#exit;
#
sub _calculate_offset {
    my $self = shift;
    my $pos = shift;
    my $read_pos = shift;
    my $cigar = shift;

    my $current_offset=0;
    my $current_pos=$read_pos;
    my @ops = $cigar =~ m/([0-9]+)([MIDNSHP])/g;
    OP:
    while(my ($cigar_len, $cigar_op) =  splice @ops, 0, 2 ) {
        my $new_offset;
        my $last_pos=$current_pos;

        if($cigar_op eq 'M') {
            $current_pos+=$cigar_len;
            $current_offset+=$cigar_len;
        }
        elsif($cigar_op eq 'I') {
            $current_offset+=$cigar_len;
        }
        elsif($cigar_op eq 'D') {
            $current_pos+=$cigar_len;

        }
        elsif($cigar_op eq 'N') {
            #this is the same as a deletion for returning a base from the read
            $current_pos += $cigar_len;
        }
        elsif($cigar_op eq 'S') {
            #soft clipping means the bases are in the read, but the position (I think) of the read starts at the first unclipped base
            #Functionally this is like an insertion at the beginning of the read
            $current_offset+=$cigar_len;
        }
        elsif($cigar_op eq 'H') {
            #hard clipping means the bases are not in the read and the position of the read starts at the first unclipped base
            #Shouldn't do anything in this case, but ignore it
        }
        else {
            die("CIGAR operation $cigar_op currently unsupported by this module");
        }

        if($pos < $current_pos && $pos >= $last_pos) {
            if($cigar_op eq 'M') {
                my $final_adjustment = $current_pos - $pos;
                return $current_offset - $final_adjustment;
            }
            else {
                return;
            }
        }
    }
    #position didn't cross the read
    return;
}

