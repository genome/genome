package Genome::Model::Tools::Analysis::Maf::AddUuids;

use strict;
use warnings;
use Genome;
use IO::File;
use JSON;


class Genome::Model::Tools::Analysis::Maf::AddUuids {
    is => 'Command',
    has => [
        maf_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'maf file containing proper TCGA barcodes',
        },

        output_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'output file',
        },
        
        ]
};


sub help_brief {
    "Queries TCGA servers for UUIDs that match your barcodes, appends the appropriate UUID columns"
}

sub help_detail {
    "Queries TCGA servers for UUIDs that match your barcodes, appends the appropriate UUID columns. It always places the UUIDs in column 33 and 34, as per the spec. If you have extra columns on your maf, those will be shifted right to accomodate the UUID columns."
}


################################################################################################

sub execute {
    my $self = shift;
    my $maf_file = $self->maf_file;
    my $output_file = $self->output_file;

    

    my %barcodes;

    my $inFh = IO::File->new( $maf_file ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        
        next if $line =~ /^#/;
        next if $line =~ /^Hugo_Symbol/;

        my @F = split("\t",$line);
        unless (($F[15] =~ /TCGA-/) && ($F[16] =~ /TCGA-/)){
            die("tcga barcodes need to be present in columns 16 and 17");
        }
        $barcodes{$F[15]} = 0;
        $barcodes{$F[16]} = 0;
    }
    close($inFh);


    #grab the results from TCGA as JSON
    my $barcodeString = join(",",keys(%barcodes));    
    my $content = `curl -k -X POST --data \"$barcodeString\" --header \"Content-Type: text/plain\" https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/barcode/batch`;

    unless($content){
        die("no content returned from API query");
    }
    if($content =~ /errorMessage/){
        die("API query failed with the following response:\n$content\n");
    }

    my $json = new JSON;
    my $json_text = $json->allow_nonref->utf8->relaxed->decode($content);
    
    #store each barcode/uuid pair in a hash
    my %uuids;
    my $i;
    foreach my $pair (@{$json_text->{uuidMapping}}){
        $uuids{$pair->{barcode}} = $pair->{uuid};
    }

    open(OUTFILE,">$output_file");
    #go back through the maf file, append UUIDS
    $inFh = IO::File->new( $maf_file ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        
        #skip def line
        if ($line =~ /^#/){
            print OUTFILE $line . "\n";
            next
        }

        my @F = split("\t",$line);

        my @extraCols;
        if(defined($F[32])){
            @extraCols=@F[32..$#F];
        }

        if ($line =~ /^Hugo_Symbol/){
            $F[32] = "Tumor_Sample_UUID";
            $F[33] = "Matched_Norm_Sample_UUID";
            print OUTFILE join("\t",(@F[0..33],@extraCols)) . "\n";
            next;
        }
 
        $F[32] = "";
        $F[33] = "";

        if(defined($uuids{$F[15]})){
            $F[32] = $uuids{$F[15]};
        }
        if(defined($uuids{$F[16]})){
            $F[33] = $uuids{$F[16]};
        }
        print OUTFILE join("\t",(@F[0..33],@extraCols)) . "\n";
    }
    close($inFh);

}




1;

