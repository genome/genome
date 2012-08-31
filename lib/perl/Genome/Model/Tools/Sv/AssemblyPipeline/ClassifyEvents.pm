
#
#
# 
#
# Parse output from svCaptureValidation.pl
# put in one file doing 
#    'cat  BRC*.csv |  grep -v OUTER  > allPatientsAllSvCalls.extend10.end1.sub1.indel1.score50.110411.csv'
#
# Partition events into somatic, germline, noEvent and notClear
#
# non-event:  $normalSv == 0 &&  $tumorSv == 0
# germline:   $tumorSv >= $MinGermlineSvReads && $normalSv >= $MinGermlineSvReads 
# somatic:    $tumorSv >= $MinTumorSvReads && $tumorSv > $normalSv && $normalSv <= $MaxNormalSvReads && $pvalue <= $Max_pValue 
#             and not germline in any other sample
# ambiguous:  none of above
#
#
#
#   If there is evidence for SV in normal, it is presumed germline so don't need to correct the number of SV-supporting reads in tumor.  
#   If there was tumor contamination in normal, then need to deal with that
#   
#   Need input file that has fraction tumor contamination in normal
#   format: $samplePrefix#  $fractionContamination
#

package Genome::Model::Tools::Sv::AssemblyPipeline::ClassifyEvents;    

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Genome;    

# can also include optional parameters is 'has => [ ]'
#     like this:   foo => { is => 'Text', doc => "", is_optional => 1, default_value => 5 },             
class Genome::Model::Tools::Sv::AssemblyPipeline::ClassifyEvents {
    is => 'Command',                       
    has => [           
    readcount_file => { is => 'Text',  doc => "Output file from 'gmt sv assembly-pipeline remap-reads'.  Output from all patients can be catenated together into one file." },
    ],

    has_optional => [
    output_file_prefix => { 
        is => 'Text',    
        doc => "'.somatic', '.germline', '.ambiguous', and '.noevent' will be appended to this prefix to produce output files. Default will be the input to '--readcount-file'.",
        #calculate_from => 'readcount_file',
        #calculate => q|$readcount_file|,
        is_optional => 1,
    },
    tumor_purity_file => { 
        is => 'Text',
        #doc => "File with two columns; patientId and fraction of tumor cells in tumor sample (1 is pure tumor).",
        doc => "Used only when there is tumor contamination in normal sample. File with two columns; patientId and tumor purity expressed as fraction of tumor cells in tumor sample (1 is pure tumor).",
        is_optional => 1 
    },
    tumor_in_normal_file => { 
        is => 'Text',
        #doc => "File with two columns; patientId and fraction of normal contamination in tumor (0 is no contamination)."
        doc => "Used only when there is tumor contamination in normal sample. File with two columns; patientId and amount of tumor contamination in normal (0 is no contamination)."
    },
    min_tumor_sv_reads => { 
        is => 'Integer',    
        doc => "Minimum number of SV-supporting reads in tumor sample to call somatic.",
        default => 5 
    },
    max_normal_sv_reads => { 
        is => 'Integer',    
        doc => "Maximum number of SV-supporting reads in normal sample to call somatic.",
        default => 5 
    },
    min_germline_sv_reads => { 
        is => 'Integer',    
        doc => "Minimum number of SV-supporting reads in tumor and normal to call germline.",
        default => 20 
    },
    max_pvalue => { 
        is => 'Number',    
        doc => "Maximum p-value from Fisher's exact comparing tumor/normal SV-supporting reads required before calling event somatic.",
        default => 0.05  
    },
    ], 
};

sub _is_hidden_in_docs { $ENV{GENOME_EXAMPLES} ? () : 1 }

sub sub_command_sort_position { -2 }

sub help_brief {                           
    "Given number of SV-supporting reads in tumor and normal, classifies events as germline, somatic or ambiguous."                 
}


sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
Requires output file from 'gmt sv assembly-pipeline remap-reads', listing number of SV-supporting reads in tumor and normal samples.  Creates four files:
 ..output_file_prefix.somatic -- somatic events
 ..output_file_prefix.germline --germline events
 ..output_file_prefix.ambiguous --events that are neither somatic or germline (or which has breakpoints identical to those of a germline event in a different patient)
 ..output_file_prefix.noevent --Ask John Wallis
EOS
}

sub execute {                              
    my $self = shift;

    #parse input parameters
    my $ReadCountFile = $self->readcount_file;
    my $MinTumorSvReads = $self->min_tumor_sv_reads;
    my $MaxNormalSvReads = $self->max_normal_sv_reads;
    my $MinGermlineSvReads = $self->min_germline_sv_reads;
    my $Max_pValue = $self->max_pvalue;

    #set up output files and temp dir
    #my $Dir = $self->output_dir;
    my $prefix = $self->output_file_prefix;
    unless (defined $prefix) { $prefix = $ReadCountFile; }
    #my $date = `date +%F`; chomp $date; $date =~ s/\-//g; if ( $date =~ /20(\d{6})/ ) { $date = $1; }
    my $tempAmbiguousFile = Genome::Sys->create_temp_file_path();
    my $tempSomaticFile = Genome::Sys->create_temp_file_path();
    my $noEventFile = $prefix . ".noevents";
    my $ambiguousFile = $prefix . ".ambiguous";
    my $finalSomaticFile = $prefix . ".somatic";
    my $germlineFile = $prefix . ".germline";

    #my $noEventFile = "$Dir/allPatientsNoEvent.$date.csv";
    #my $ambiguousFile = "$Dir/allPatientsAmbiguousAndPresumedGermline.$date.csv";
    #my $tempAmbiguousFile = "$Dir/tempAmbiguousFile.$date.csv";
    #my $tempSomaticFile = "$Dir/tempSomaticFile.$date.csv";
    #my $finalSomaticFile = "$Dir/allPatientsSomatic.$date.csv";
    #my $germlineFile = "$Dir/allPatientsGermline.$date.csv";

    open(NON, "> $noEventFile");
    open(AMB, "> $tempAmbiguousFile");
    open(SOM, "> $tempSomaticFile");
    open(GERM, "> $germlineFile");

    my ( @entireFile, $line, $normalTotal, $tumorTotal, $normalSv, $tumorSv, $chrA, $bpA, $chrB, $bpB, $event, 
        $pvalue, %germline, $patientId );

    # Get a ref to hash that has patient ID and fraction of tumor contamination in normal
    my ( $FractionTumorInNormal, $TumorPurity );
    if ( defined $self->tumor_purity_file && defined $self->tumor_in_normal_file ) {
        $FractionTumorInNormal = getContamination($self->tumor_in_normal_file);
        $TumorPurity = getContamination($self->tumor_purity_file);
    }

    my $header;
    open(IN, "< $ReadCountFile") || die "Could not open '$ReadCountFile': $!";
    @entireFile = <IN>;
    close IN;
    foreach $line (@entireFile) {
        if ( $line =~ /^#/ ) { $header = $line; print NON $header; print GERM $header; next; }
        chomp $line;
        if ( $line =~ /no\s+fasta\s+sequence/ ) { print NON "$line\n"; next; }
        $normalTotal = $patientId = $tumorTotal = $normalSv = $tumorSv  = undef;
        if ( $line =~ /(\S+).normal\.totalReads\:(\d+)/i ) { $patientId = $1; $normalTotal = $2; }
        if ( $line =~ /$patientId.tumor\.totalReads\:(\d+)/i ) { $tumorTotal = $1; }
        if ( $line =~ /$patientId.normal.svReadCount\:(\d+)/i ) { $normalSv = $1; }
        if ( $line =~ /$patientId.tumor.svReadCount\:(\d+)/i ) { $tumorSv = $1; }
        ( $normalTotal =~ /\d+/ && $tumorTotal =~ /\d+/ && $normalSv =~ /\d+/ && $tumorSv  =~ /\d+/ ) ||
        die "Did not get tumor and/or normal total reads and SV read count: '$line'";

        # if there are no SV-supporting reads, it is a non-event
        if ( $normalSv == 0 &&  $tumorSv == 0 ) { 
            print NON "$line\t$patientId\n";
            next;
        }

        # Correct normal SV read count for the number that is due to contamination of normal sample with tumor cells
        my ( $correctedNormalSv, );
        if ( defined $FractionTumorInNormal && defined $TumorPurity ) {
            (defined $$FractionTumorInNormal{$patientId}) || die "Did not get amount of tumor contamination in normal for '$patientId'";
            $correctedNormalSv = adjustReadCount($tumorTotal, $tumorSv, $normalTotal, $normalSv, $$FractionTumorInNormal{$patientId}, $$TumorPurity{$patientId});
        } else {
            $correctedNormalSv = $normalSv;
        }

        # if neither tissue has minimum number of SV reads, it is "can't tell"
        if ( $correctedNormalSv < $MinTumorSvReads &&  $tumorSv < $MinTumorSvReads  ) { 
            print AMB "$line\t$patientId\tNormalCorrectedSvReadCount:$correctedNormalSv\t-\n";
            next;
        }

        $pvalue = pValue($normalTotal, $correctedNormalSv, $tumorTotal, $tumorSv);

        # Somatic: 
        if ( $tumorSv >= $MinTumorSvReads && $tumorSv > $correctedNormalSv && $correctedNormalSv <= $MaxNormalSvReads && $pvalue <= $Max_pValue  ) {
            print SOM "$line\t$patientId\tNormalCorrectedSvReadCount:$correctedNormalSv\tTumor_normal:$pvalue\n";

            # Germline
        } elsif ( $tumorSv >= $MinGermlineSvReads && $correctedNormalSv >= $MinGermlineSvReads ) {
            print GERM "$line\t$patientId\tNormalCorrectedSvReadCount:$correctedNormalSv\tTumor_normal:$pvalue\n";

            # Ambiguous.  Not clearly germline and not clearly somatic
        } else {
            print AMB "$line\t$patientId\tNormalCorrectedSvReadCount:$correctedNormalSv\tTumor_normal:$pvalue\n";
        }
    }

    close NON; close GERM; close SOM;


    # Now see if any of the somatic calls are identical to any of the germline ones
    # If so, remove them from somatic file and put in ambiguous file
    open(GERM, "< $germlineFile");
    @entireFile = <GERM>;
    close GERM;
    foreach $line (@entireFile) {
        chomp $line;
        # 2.28	2	99234794	99234794	2	99235019	99235019	DEL
        (undef, $chrA, $bpA, undef, $chrB, $bpB, undef, $event) = split /\s+/, $line;
        $germline{$chrA}{$bpA}{$chrB}{$bpB}{$event} = 1;
    }
    open(SOM, "< $tempSomaticFile");
    open(NEW_SOM, "> $finalSomaticFile");
    print NEW_SOM $header;
    @entireFile = <SOM>;
    close SOM;
    foreach $line (@entireFile) {
        chomp $line;
        # 2.28	2	99234794	99234794	2	99235019	99235019	DEL
        (undef, $chrA, $bpA, undef, $chrB, $bpB, undef, $event) = split /\s+/, $line;
        if ( defined $germline{$chrA}{$bpA}{$chrB}{$bpB}{$event} ) {
            print AMB "$line\tsomatic_germline\n";
        } else {
            print NEW_SOM "$line\n";
        }
    }
    close NEW_SOM; close AMB; 

    # Now see which of the ambiguous could be called germline and annotate them as such
    open(AMB, "< $tempAmbiguousFile") || die "Could not open '$tempAmbiguousFile' for input: $!";
    @entireFile = <AMB>;
    close AMB;
    open(NEW_AMB, "> $ambiguousFile") || die "Could not open '$ambiguousFile' for output: $!";
    print NEW_AMB $header;
    foreach $line (@entireFile) {
        chomp $line;
        # 2.28	2	99234794	99234794	2	99235019	99235019	DEL
        (undef, $chrA, $bpA, undef, $chrB, $bpB, undef, $event) = split /\s+/, $line;
        print NEW_AMB "$line";
        if ( defined $germline{$chrA}{$bpA}{$chrB}{$bpB}{$event} && $line !~ /somatic_germline/ ) {
            print NEW_AMB "\tgermline";
        } else {
            print NEW_AMB "\t-";
        }
        print NEW_AMB "\n";
    }
    close NEW_AMB;

    return 1;                              
}



sub pValue {
    # The world's most inefficient way to get a pvalue
    # Input: total normal reads, normal SV reads, total tumor reads, tumor SV reads
    # Return: pvalue from Fisher's Exact test to see if the SV reads from
    #         tumor and normal are significantly different given the total number of
    #         reads from each

    my ($normalTotal, $normalSv, $tumorTotal, $tumorSv) = @_;

    my $pvalue;
    my $rFile = "/tmp/rFile".rand();
    open(OUT, "> $rFile") || die "Could not open '$rFile': $!";
    print OUT "x=matrix(c($tumorTotal,$tumorSv,$normalTotal,$normalSv),ncol=2)\nfisher.test(x)\n";
    close OUT;
    open(ROUT, "R --silent --slave <  $rFile |");
    my @rOutput = <ROUT>;
    close ROUT;
    foreach (@rOutput) { 
        if ( $_ =~ /p-value\s+[=<]\s+(\S+)/ ) { $pvalue = $1; }
    }
    unlink $rFile;

    return $pvalue;
}



sub adjustReadCount {
    # Subtract the number of SV-supporting reads in normal that are thought to be due to tumor contamination
    my ( $tumorTotal, $tumorSv, $normalTotal, $normalSv,  $fractionTumorInNormal, $tumorPurity ) = @_;

    my ( $svReadsInNormalDueToTumor, $adjustedNormalSvReadCount  );
    if ( $tumorPurity == 0 || $tumorTotal == 0 ) {
        $svReadsInNormalDueToTumor = 0;
    } else {
        $svReadsInNormalDueToTumor = $fractionTumorInNormal/$tumorPurity * $tumorSv/$tumorTotal * $normalTotal;
    }
    $adjustedNormalSvReadCount = $normalSv - $svReadsInNormalDueToTumor; 
    if ( $adjustedNormalSvReadCount < 0 ) { $adjustedNormalSvReadCount = 0; }

    # Round off to nearest integer
    $adjustedNormalSvReadCount += 0.5;
    $adjustedNormalSvReadCount = int($adjustedNormalSvReadCount );

    return $adjustedNormalSvReadCount;
=head
I would do:

expected number of sv reads in normal due to tumor = [(Number SV in Tumor / total reads in tumor ) * (TumorContaminationOfNormal / 
TumorPurity) * NumberTotalNormalReads]

ie if there are:
100 reads in the normal and tumor
10 SV reads in the tumor
Tumor purity of 0.5
Tumor contamination of the normal of 0.1

I would calculate the number of expected SV reads in the normal as:

(10/100) * (0.1/0.5) * 100 = 2
=cut
}

sub getContamination {
    # Input: file with two columns; sample name and fraction contamination
    # Return: ref to hash with key = sample name, value = fraction contamination
    my $contaminationFile = $_[0];
    my ( @entireFile, $line, $fraction, $sampleName, %sampleToFraction );
    open(IN, "< $contaminationFile") || confess "Could not open '$contaminationFile': $!";
    @entireFile = <IN>;
    close IN;
    foreach $line ( @entireFile ) {
        chomp $line;
        if ( $line =~ /^\s*$/ ) { next; }
        ($sampleName, $fraction) = split /\s+/, $line;
        #if ( $sampleName =~ /AML0(\d)/ ) { $sampleName = "AML$1"; }
        ($fraction =~ /0?\.\d+/ || $fraction == 0 || $fraction == 1 ) 
        || confess "Unexpected format for fraction in '$contaminationFile': \$fraction = '$fraction' in '$line'";
        if ( defined $sampleToFraction{$sampleName} && $sampleToFraction{$sampleName} != $fraction ) {
            confess "Two values for '$sampleName' in '$contaminationFile' '$fraction' and '$sampleToFraction{$sampleName}'";
        }
        $sampleToFraction{$sampleName} = $fraction;
    }
    return \%sampleToFraction;
}

1;
