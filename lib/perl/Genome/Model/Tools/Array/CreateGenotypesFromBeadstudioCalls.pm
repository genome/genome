package Genome::Model::Tools::Array::CreateGenotypesFromBeadstudioCalls;

use strict;
use warnings;

use Genome;
use Command;
use Text::CSV_XS;
use Sort::Naturally qw( nsort );

class Genome::Model::Tools::Array::CreateGenotypesFromBeadstudioCalls {
    is => 'Command',
    has => [
    genotype_file => 
    {
        type => 'String',
        is_optional => 0,
        doc => "File of forward strand describing which allele of the probe was detected",
    },
    output_directory =>
    {
        type => 'String',
        is_optional => 0,
        doc => "Directory to place individual files for each header in the file",
    },
    ucsc_array_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => 'File from UCSC indicating which strand and alleles the Infinium platform probes'
    },
    #illumina_manifest_file =>
    #{
    #    type => 'String',
    #    is_optional => 0,
    #    doc => 'File from Illumina indicate genotypes and stranding information',
    #},
    ]
};


sub execute {
    my $self=shift;

    #TODO Some basic file checks
    $| = 1;
    my $call_href = $self->create_call_file_hash;
    unless($call_href) {
        $self->error_message("Unable to parse genotype file and create call file hash");
        return;
    }
    #my $probe_href = $self->create_probe_hash;
    $self->status_message("Finished creating call file hash");
    #Store the header information for creating filenames later
    my @filehandles = split /\t/, $call_href->{ID};
    for my $name (@filehandles) {
        $name =~s/\r//;
    }

    $call_href = $self->convert_to_genotype($call_href);#, $probe_href);

    #Convert each call into a genotype and write to a new file by chromosome and position
    #create file handles to write out each sample
    my $out_dir = $self->output_directory;
    for my $file (@filehandles) {
        my $filehandle = new IO::File "$out_dir/$file.genotype", "w";
        unless(defined($filehandle)) {
            $self->error_message("Couldn't open filehandle for file: $out_dir/$file.genotype");
            return;
        }   
        $file = $filehandle;
    }

    for my $chromosome (nsort keys %$call_href) {
        for my $position (sort {$a <=> $b} keys %{$call_href->{$chromosome}}) {
            my $i;
            my @calls = split /\t/, $call_href->{$chromosome}{$position};
            for($i = 0; $i < scalar(@calls); $i++) {
                print {$filehandles[$i]} "$chromosome\t$position\t",$calls[$i],"\n";
            }
        }
    }

    map { $_->close; } @filehandles;


    return 1;
}

sub create_call_file_hash {
    my $self = shift;
    my $file = $self->genotype_file;
    my %call_hash;

    my $fh = new IO::File "$file", "r";
    unless(defined($fh)) {
        $self->error_message("Unable to open $file for reading");
        return;
    }
    #Skip through the header lines
    my $observed_content_header = 0;
    while(my $line = $fh->getline) {
        last if $line =~ /\[Data\]/;
    }
    if($fh->eof) {
        $self->error_message("Unexpected file format. Never encountered [Data] flag");
        return;
    }
    my $expected_calls = undef;

    while(my $line = $fh->getline) {
        chomp $line;
        my ($ID, @calls) = split /[\t,]/, $line; #supporting comma or tab separated. APparently this varies
        map {$_ =~ s/(\D+)\|.*$/$1/g;} @calls;

        if(defined($expected_calls) && $expected_calls != scalar(@calls)) {
            $self->error_message("Unexpected number of calls");
            return;
        }
        else {
            $expected_calls = scalar(@calls);
        }
        $ID = $ID eq q{} ? 'ID' : $ID;  #make sure that on first line there is an actual label
        my $call_string = join "\t", @calls; #doing this because arrays in perl apparently utilize an obscene amount of mem. ~38 bytes per element in this case. Which means for this program which has about ~1million entries. We will be using almost 3.5 G of ram just to store eveything in an array when addining into the hash, the overhead is much greater. In contrast to store the whole thing as a string, though computationally a bit expensive, results in much less overhead.
        $call_hash{$ID} = $call_string; 
    }
    return \%call_hash;
}

sub create_probe_hash {
    my ($self) = @_;

    my $file = $self->illumina_manifest_file;
    my $fh = new IO::File "$file", "r";
    unless(defined($fh)) {
        return 0;
    }

    my %probe;

    while(my $line = $fh->getline) {
        last if $line =~ /^\[Assay\]/;    
    }

    #skip final header line
    $fh->getline;
    
    while(my $line = $fh->getline) {
        my ($illumina_id, $name, $illumina_strand, $alleles) = split /,/, $line;
        $probe{$name} = { id => $illumina_id,
            strand => substr($illumina_strand, 0, 1),
            alleles => $alleles,
        };
    }

    return \%probe;
}

sub convert_to_genotype {
    my ($self, $calls,) = @_;

    my $csv = new Text::CSV_XS({sep_char => "\t"}); #tab separated
    my $file = $self->ucsc_array_file;
    my $afh = new IO::File "$file","r";

    my %new_calls;
    my %prev_scanned_probe;

    while(my $line = <$afh>) {
        chomp ($line);    

        #File is of the format Name Chr Position

        $csv->parse($line);

        my ($bin,$chr, $start0, $pos, $snp_id, $score, $strand, $observed_alleles) = $csv->fields();

        next if($bin =~ /^\#bin/xi); #skip header
        $chr =~ s/^chr//;   #adjust from UCSC chromosome notation
        $chr =~ s/^M$/MT/; #further adjust chr
        if(exists($calls->{$snp_id}) && !exists($prev_scanned_probe{$snp_id})) {
            #Expecting forward strand calls
            #check that this is the case
            #
            #THIS DID NOT HELP, BUT THE LOGIC SEEMS LIKE IT MAY PROVE USEFUL AT A LATER DATA
            #
            #if(exists($probe->{$snp_id})) {
            #    #determine forward strand
            #    my ($name_il_strand, $name_dbSNP_strand) = $probe->{$snp_id}{id} =~ /\_(\D)\_(\D)\_/;
            #    if($name_dbSNP_strand eq 'R') {
            #        $alleles_altered = $name_il_strand eq $probe->{$snp_id}{strand} ? 1 : 0;
            #    }
            #    elsif($name_dbSNP_strand eq 'F') {
            #        $alleles_altered = $name_il_strand eq $probe->{$snp_id}{strand} ? 0 : 1;
            #    }
            #    else {
            #        $self->error_message("Parsing error. Unknown dbSNP strand");
            #        die;
            #    }

            #}
            #else {
            #    $self->error_message("Probe id $snp_id not found in Illumina probe file");
            #    die;
            #}
            ##UCSC doesn't adjust
            #if($alleles_altered && !$self->is_dbsnp($snp_id)) {
            #    $observed_alleles =~ tr/ACTGactg/TGACtgac/;
            #}
                


            ##Check alleles to make sure that what we're getting from Illumina matches what we get from UCSC
            #if($alleles_altered) {
            #   $observed_alleles =~ tr/ACTGactg/TGACtgac/;
            #}
            if($self->is_dbsnp($snp_id)) {
                #ucsc only retrieves things from dbSNP. THey don't alter anything else
                unless($self->contains_expected_alleles($observed_alleles, $calls->{$snp_id})) {
                    my $alleles = $calls->{$snp_id};                                
                    $self->error_message("Unexpected alleles for probe $snp_id. Expected $observed_alleles. Got $alleles");
                    next;
                }
                if($strand eq '-' ) {
                    #adjust stranding to +/- instead of forward/reverse
                    #it's on the - strand in our file
                    $calls->{$snp_id} =~ tr/ACTGactg/TGACtgac/;
                }
            }
            $new_calls{$chr}{$pos} = $calls->{$snp_id};
            delete $calls->{$snp_id};
        }
        if(exists($prev_scanned_probe{$snp_id})) {
            my @first_pos = @{$prev_scanned_probe{$snp_id}};
            delete $new_calls{$first_pos[0]}{$first_pos[1]};
        }
        $prev_scanned_probe{$snp_id} = [$chr,$pos];

    }
    return \%new_calls;
}

sub contains_expected_alleles {
    my ($self, $expected_alleles, $reported_calls_string) = @_;

    my %expected_alleles = map { uc($_) => 1 } split /\//, $expected_alleles;
    my @reported_alleles = grep {$_ !~ /\s+/ } map { split // } $reported_calls_string;
    foreach my $allele (@reported_alleles) {
        unless(exists($expected_alleles{uc($allele)}) || $allele eq '-') {
            return;
        }
    }
    return 1;
}

sub is_dbsnp {
    my ($self, $probe_name) = @_;
    if($probe_name =~ /^rs/) {
        return 1;
    }
    else {
        return 0;
    }
}


1;

sub help_brief {
    "Converts an Illumina Beadstudio Forward Strand Report file to genotype file";
}

