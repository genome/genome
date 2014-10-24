package Genome::InstrumentData::AlignmentResult::MultiAligner;

use strict;
use warnings;
use List::MoreUtils 'any';

use Genome;

class Genome::InstrumentData::AlignmentResult::MultiAligner {
    is => 'Genome::InstrumentData::AlignmentResult',
    
    has_constant => [
        aligner_name => { value => 'multi aligner', is_param=>1 },
    ],
    has_transient_optional => [
        any_fillmd => { is => 'Integer' },
    ],
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[tmp>90000 && mem>10000] span[hosts=1] rusage[tmp=90000, mem=10000]' -M 10000000 -n 4";
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return "multi-aligner " . $self->aligner_params;
}

sub _run_aligner {
    my $self = shift;
    my @aligners = split /\s*}\s*/, $self->aligner_params or die "Didn't get any aligners to align with";
    print "\nStarting multi-alignment\n\n";

    for my $aln (@aligners) {
        my $aligner = $self->get_alignment_module($aln,@_) or die "Failed to make '$aln' alignment module: $!";

        *Genome::InstrumentData::AlignmentResult::Bwa::supports_streaming_to_bam = sub { 0};
    
        print join("\n", ('-=' x 20, "Starting $aln alignment",'-=' x 20)),"\n";
        $aligner->_run_aligner(@_);
        $self->any_fillmd(1) if $aligner->fillmd_for_sam;
        $aligner->delete;
        # unaligned reads => new input fastqs
        @_ = $self->extract_unaligned_reads(scalar(@_)) unless $aln eq $aligners[-1];
    }
    
    print "\nFinishing multi-alignment\n\n";
    return 1;
}

sub extract_unaligned_reads {
    my $self = shift;
    my $num_fastqs = shift;
    print "Extracting FastQ files from all_sequences.sam ...\n";
    my $tmp_dir = $self->temp_scratch_directory;
    my %readnames = ();
    open FORWARD, ">$tmp_dir/unaligned_fwd.fq" or die "Couldn't open file for writing: $!";
    open REVERSE, ">$tmp_dir/unaligned_rev.fq" or die "Couldn't open file for writing: $!";
    if (-s "$tmp_dir/all_sequences_unaligned.sam"){
        open SAMFILE, "$tmp_dir/all_sequences_unaligned.sam" or die "Couldn't open sam file. Strange: $!";
    } else {
        open SAMFILE, "$tmp_dir/all_sequences.sam" or die "Couldn't open sam file. Strange: $!";
    }
    while (<SAMFILE>){
        (my $name,my $flags) = split;
        next unless ($flags & 0x0004) || ($flags & 0x0008);
        # book-keeping to make sure we don't add duplicates
        if ($readnames{$name}){
            if ($readnames{$name}{$flags}){
                next;
            } else {
                $readnames{$name}{$flags} = 1;
            }
        } else {
            $readnames{$name} = {$flags=>1};
        }
        # direct it to file
        if ($flags & 0x40){ # first of the pair
            print FORWARD sam2fastq(1,$_);
        } elsif ($flags & 0x80){ # second of the pair
            print REVERSE sam2fastq(2,$_);
        } else { #if ($num_fastqs == 1) { # only aligning one
            print FORWARD sam2fastq(1,$_);
        }# else {
        #    die "Bad flag format: $flags";
        #}
    }
    close FORWARD;
    close REVERSE;
    close SAMFILE;
    unlink "$tmp_dir/all_sequences_unaligned.sam"; # just in case, so we don't accidentally read the same unaligned next time
    return ("$tmp_dir/unaligned_fwd.fq",) if $num_fastqs == 1;
    return ("$tmp_dir/unaligned_fwd.fq","$tmp_dir/unaligned_rev.fq");
}

sub sam2fastq {
    my $readnum = shift; # 1 for fwd, 2 for rev
    (my $name,my $read,my $qual) = (split)[0,9,10];
    $name .= "/$readnum"; # note: doesn't mess with the #0 ending (might need to be first 5 chars of fwd read)
    return join "\n", ('@'.$name,$read,'+'.$name,$qual,'')
}

sub get_alignment_module {
    my $self = shift;
    my $aln = shift;
    my ($aln_name, @opts) = $aln =~ /\S*['"].*?['"]|\S+/g;
    my %params = ();
    foreach (@opts) {
        my ($key,$val) = split "=";
        $val =~ s/^['"]|['"]$//g; # strip any quotes
        if ($key =~ /version/){
            $params{aligner_version} = $val;
        } elsif ($key =~ /params/){
            $params{aligner_params} = $val;
        } else {
            print "Warning: Ignoring unknown param: $_\n";
        }
    }
    my $aln_class = Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aln_name) or die "Not a valid AlignmentResult: $aln_name";
    # for whatever reason, it dies unless you have aligner_version set, even though a default exists
    unless (defined($params{aligner_version})){
        $params{aligner_version} = "Genome::Model::Tools::$aln_class"->default_version;
    }
    my $alignment_result_class_name = "Genome::InstrumentData::AlignmentResult::$aln_class";
    eval "use $alignment_result_class_name";
    # sneaky way to instantiate an AlignmentResult without starting all over
    my $alignment_result = $alignment_result_class_name->__define__(
                                          instrument_data_id => $self->instrument_data->id,
                                          aligner_name => $aln_name,
                                          reference_build => $self->reference_build, 
                                          temp_scratch_directory => $self->temp_scratch_directory,
                                          temp_staging_directory => $self->temp_staging_directory,
                                          %params );

    $alignment_result->lookup_hash($alignment_result->calculate_lookup_hash);
    return $alignment_result;
}

# true if any of the component aligners' fillmd is true
sub fillmd_for_sam { 
    my $self = shift;
    return $self->any_fillmd;
}

sub _check_read_count {
    return 1;
}

sub supports_streaming_to_bam {
    0;
}
