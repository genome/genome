package Genome::Model::Tools::Array::CreateGenotypesFromAffyCalls;

#use strict;
use warnings;

use Genome;
use Command;
use Text::CSV_XS;
use Sort::Naturally qw( nsort );

class Genome::Model::Tools::Array::CreateGenotypesFromAffyCalls {
    is => 'Command',
    has => [
    annotation_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "Annotation file for the array information being converted",
    },
    call_file => 
    {
        type => 'String',
        is_optional => 0,
        doc => "File of calls describing which allele of the probe was detected",
    },
    output_filename => {
        type => 'String',
        is_optional => 0,
        doc => "Filename of the output file",
    },

    ]
};


sub execute {
    my $self=shift;

    #TODO Some basic file checks
    print "about to call create_call_file_hash\n";
    my $call_href = $self->create_call_file_hash;

    $call_href = $self->convert_to_genotype($call_href);

    #Convert each call into a genotype and write to a new file by chromosome and position
    #create file handles to write out each sample
    my $file = $self->output_filename;
    my $filehandle = new IO::File "$file", "w";
    unless(defined($filehandle)) {
        $self->error_message("Couldn't open filehandle for file: $file");
        return;
    }   
    for my $chromosome (nsort keys %$call_href) {
        for my $position (sort {$a <=> $b} keys %{$call_href->{$chromosome}}) {
            print $filehandle "$chromosome\t$position\t",$call_href->{$chromosome}{$position},"\n";
        }
    }

   $filehandle->close; 
    
    return 1;
}

sub create_call_file_hash {
    my $self = shift;
    my $file = $self->call_file;
    my %call_hash;

    my $fh = new IO::File "$file", "r";
    unless(defined($fh)) {
        return 0;
    }
    #This is somewhat different as the convention seems to have switched to having a single file per sample
    while(my $line = $fh->getline) {
        next if $line =~ /^#%/; #these are comments in affy apt-chp-to-text output

        my ($PS_ID, $call, $confidence) = split /\s+/, $line; #confidence may not be in file depending on format
        $call_hash{$PS_ID} = $call; 
    }
    return \%call_hash;
}

sub convert_to_genotype {
    my ($self, $calls) = @_;

    my $csv = new Text::CSV_XS;
    my $file = $self->annotation_file;
    my $afh = new IO::File "$file","r";

    my %new_calls;
    
    LINE:
    while(my $line = <$afh>) {
        chomp ($line);    

        #Skip comment lines
        next if($line =~ /^ \# (.*) $/xi);

        # Header
        #
        # "Probe Set ID","Affy SNP ID","dbSNP RS ID","Chromosome","Physical
        # Position","Strand","ChrX pseudo-autosomal region
        # 1","Cytoband","Flank","Allele A","Allele B","Associated Gene","Genetic
        # Map","Microsatellite","Fragment Enzyme Length Start Stop","Allele
        # Frequencies","Heterozygous Allele Frequencies","Number of
        # individuals/Number of chromosomes","In Hapmap","Strand Versus dbSNP","Copy
        # Number Variation","Probe Count","ChrX pseudo-autosomal region 2","In Final
        # List","Minor Allele","Minor Allele Frequency"

        $csv->parse($line);

        my
        ($snp_id,$as_id,$dbsnp,$chr,$phys,$strand,$pseudo,$cyto,$flank,$alleleA,$alleleB,$gene,$field)
        = $csv->fields();

        #next if($snp_id !~ /^$snprange$/xi);
        next if $phys =~ /\-/; #exclude ambiguous sites

        if(exists($calls->{$snp_id})) {
            my $call = $calls->{$snp_id};
            if(($call eq 'AA')||($call eq  '0')) {
                $call = $alleleA x 2;
            }
            elsif(($call eq 'AB')||($call eq '1')) {
                $call = $alleleA . $alleleB;
            }
            elsif(($call eq 'BB')||($call eq '2')) {
                $call = $alleleB x 2;
            }
            elsif(($call eq 'NC')||($call eq '-1')) {
                $call = '--';
            }
            else {
                warn "Unrecognized genotype call $call. Skipping...\n";
                next LINE;
            }
            $call =~ tr/ACTGactg/TGACtgac/ if $strand eq '-';
            $new_calls{$chr}{$phys} = $call;
            delete $calls->{$snp_id};
        }

    }
    return \%new_calls;
}



1;

sub help_brief {
    "Converts Affy genotype call file into actual base calls"
}
