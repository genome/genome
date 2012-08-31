package Genome::Model::Tools::Vcf::JoinVcf;

use strict;
use warnings;
use Genome;
use File::stat;
use File::Basename;
use DateTime;
use POSIX;
use Sort::Naturally;

use constant {
    CHROM => 0,
    POS => 1,
    ID => 2,
    REF => 3,
    ALT => 4,
    QUAL => 5,
    FILTER => 6,
    INFO => 7,
    FORMAT => 8,
    SAMPLE => 9,
};

class Genome::Model::Tools::Vcf::JoinVcf {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            is_optional => 0,
            doc => "Output merged VCF",
        },

        vcf_file_a => {
            is => 'Text',
            is_optional => 0,
            is_input => 1,
            doc => 'First file to merge',
        },

        vcf_file_a_source => {
            is => 'Text',
            is_optional => 0,
            is_input => 1,
            doc => 'Source of the first file to merge',
        },

        vcf_file_b => {
            is => 'Text',
            is_optional => 0,
            is_input => 1,
            doc => 'Second file to merge',
        },

        vcf_file_b_source => {
            is => 'Text',
            is_optional => 0,
            is_input => 1,
            doc => 'Source of the second file to merge',
        },

        intersection => {
            is => 'Boolean',
            doc => 'Set this to cause non-passing filter records to be propagated above passing',
            is_input => 1,
            default => 0,
        },
        use_gzip_files => {
            is => 'Boolean',
            doc => "Set this to use gzip input files and output gzip data",
            is_input => 1,
            default => 0,
        },
    ],
    has_transient_optional => [
        _vcf_a_fh => {
            is => 'IO::File',
            doc => 'File handle to vcf_a',
        },
        _vcf_b_fh => {
            is => 'IO::File',
            doc => 'File handle to vcf_b',
        },
        _output_fh => {
            is => 'IO::File',
            doc => 'File handle to output_file',
        },
        _header => {
            doc => 'This is a hash containing the combined header',
        },
        _format_fields => {
            doc => 'The fields used to determine the format fields in the output'
        },
        _info_fields => {
            doc => 'The fields used to determine the info fields in output',
        },
        _a_source => {
            doc => 'The detector that made the call for vcf_a',
        },
        _b_source => {
            doc => 'The detector that made the call for vcf_a',
        },
    ],
};


sub help_synopsis {
    <<'HELP';
Merge multiple VCFs - keep the FORMAT lines from files in desc order.
HELP
}


sub execute {
    my $self = shift;

    my $vcf_a = $self->vcf_file_a;
    unless(-s $vcf_a){
        die $self->error_message("Could not locate vcf_file_a at: ".$vcf_a);
    }
    my $vcf_b = $self->vcf_file_b;
    unless(-s $vcf_b){
        die $self->error_message("Could not locate vcf_file_b at: ".$vcf_b);
    }

    my ($vcf_a_fh, $vcf_b_fh, $output_fh);
    my $output = $self->output_file;
    if ($self->use_gzip_files) {
        $vcf_a_fh = Genome::Sys->open_gzip_file_for_reading($vcf_a);
        $vcf_b_fh = Genome::Sys->open_gzip_file_for_reading($vcf_b);
        $output_fh = Genome::Sys->open_gzip_file_for_writing($output);
    } else {
        $vcf_a_fh = Genome::Sys->open_file_for_reading($vcf_a);
        $vcf_b_fh = Genome::Sys->open_file_for_reading($vcf_b);
        $output_fh = Genome::Sys->open_file_for_writing($output);
    }

    $self->_vcf_a_fh($vcf_a_fh);
    $self->_vcf_b_fh($vcf_b_fh);
    $self->_output_fh($output_fh);

    $self->merge_headers; 
    $self->process_format;
    $self->process_info;
    $self->print_header;
    
    $self->process_records;

    $self->_vcf_a_fh->close;
    $self->_vcf_b_fh->close;
    $self->_output_fh->close;

    return 1;
}

sub process_info {
    my $self = shift;
    my $h = $self->_header;
    my %info_fields;
    for my $f (sort(keys(%{$h->{INFO}}))){
        my $line = $h->{INFO}->{$f};
        my (undef,$id) = split /\=/,$line;
        ($id,undef) = split /\,/, $id;
        $info_fields{$id}=1;
        print "ID=".$id."\n";
    }
    $self->_info_fields(\%info_fields);

    return 1;
}

sub process_format {
    my $self = shift;
    my $h = $self->_header;
    my %format_fields;
    for my $f (sort(keys(%{$h->{FORMAT}}))){
        my $line = $h->{FORMAT}->{$f};
        my (undef,$id) = split /\=/,$line;
        ($id,undef) = split /\,/, $id;
        $format_fields{$id}=1;
    }
    $self->_format_fields(\%format_fields);

    return 1;
}

#smartly combine two headers and stuff a hash of it into $self->_header
sub merge_headers {
    my $self = shift;
    my %header;
    my %filter;
    my %info;
    my %format;
    my @file_handles = ($self->_vcf_a_fh,$self->_vcf_b_fh);
    FILE: for my $fh (@file_handles){
        while(my $line = $fh->getline){
            chomp $line;
            if($line =~ m/^##/){
                $line =~ s/^##//;
                my ($tag,@data) = split /\=/, $line;
                my $data = join("=",@data);
                if(exists($header{$tag})){
                    if($header{$tag} eq $data){
                        next;
                    } else {
                        if($tag =~ m/FILTER/){
                            my $key = $data[1];
                            if(exists($header{FILTER}{$key})){
                                unless($header{FILTER}{$key} eq $data){
                                    #TODO turn this on when it's "safe" to do so...
                                    #die $self->error_message("Cannot merge FILTER tags.");
                                }
                            } else {
                                $header{FILTER}{$key} = $data;
                            }
                        } elsif ( $tag =~ m/FORMAT/){
                            my $key = $data[1];
                            if(exists($header{FORMAT}{$key})){
                                unless($header{FORMAT}{$key} eq $data){
                                    die $self->error_message("Cannot merge FORMAT tags.");
                                }
                            } else {
                                $header{FORMAT}{$key} = $data;
                            }

                        } elsif ( $tag =~ m/INFO/){
                            my $key = $data[1];
                            if(exists($header{INFO}{$key})){
                                unless($header{INFO}{$key} eq $data){
                                    die $self->error_message("Cannot merge INFO tags.");
                                }
                            } else {
                                $header{INFO}{$key} = $data;
                            }

                        } elsif($tag =~ m/source/){
                            $header{$tag} .= ",".$data;
                        } elsif ($tag =~ m/fileformat/){
                            die $self->error_message("Cannot continue, trying to merge files with different VCF versions! .. see header tags \"fileformat\"");
                        }
                    }
                } else {
                    my $key;
                    if($tag =~ m/FILTER/){ 
                        $header{$tag} = \%filter;
                        $key = $data[1];
                        $header{FILTER}{$key} = $data;
                    }elsif ($tag =~ m/FORMAT/){ 
                        $header{$tag} = \%format;
                        $key = $data[1];
                        $header{FORMAT}{$key} = $data;
                    }elsif ($tag =~ m/INFO/){ 
                        $header{$tag} = \%info;
                        $key = $data[1];
                        $header{INFO}{$key} = $data;
                    } elsif ($tag =~ m/filedate/i){
                        my $dt = DateTime->now;
                        my $month = $dt->month;
                        if($month < 10){
                            $month = "0".$month;
                        }
                        my $day = $dt->day;
                        if($day <10){
                            $day = "0".$day;
                        }
                        $header{$tag} = $dt->year.$month.$day;
                    } else {
                        $header{$tag} = $data;
                    }
                }
            } else {
                if( $line =~ m/^#CHROM/){
                    $line =~ s/^#//;
                    my $header = $line;
                    if(exists($header{CHROM})){
                        unless($header eq $header{CHROM}){
                            die $self->error_message("Column headers are not equal!");
                        }
                    } else {
                        $header{CHROM} = $header;
                    }
                    next FILE;
                } else {
                    die $self->error_message("Parsed header but did not find column headers!");
                }
            }
        }
    }
    my $info = "<ID=VC,Number=.,Type=String,Description=\"Variant caller\">";
    my $key = "VC,Number";
    $header{INFO}{$key} = $info;
    
    $self->_header(\%header);
    return 1;
}

# Grab a hash of the header from $self->_header and print it to the _output_fh
sub print_header {
    my $self = shift;
    my %header = %{$self->_header};
    my $h = \%header;
    my $fh = $self->_output_fh;

    $self->print_and_delete_from_hash($h,"fileformat",$fh);
    $self->print_and_delete_from_hash($h,"fileDate",$fh);
    $self->print_and_delete_from_hash($h,"source",$fh);
    $self->print_and_delete_from_hash($h,"reference",$fh);
    $self->print_and_delete_from_hash($h,"phasing",$fh);
    $self->print_and_delete_from_hash($h,"FILTER",$fh);
    $self->print_and_delete_from_hash($h,"FORMAT",$fh);
    $self->print_and_delete_from_hash($h,"INFO",$fh) if defined $h->{INFO};
    $self->print_and_delete_col_header($h,$fh);

    unless( !keys %header){
        die $self->error_message("Found some values left in header hash: ".Data::Dumper::Dumper($h));
    }
    return 1;
}

#utility function for printing the header
sub print_and_delete_from_hash {
    my $self = shift;
    my $hash = shift;
    my $key = shift;
    my $fh = shift;

    my $out = $hash->{$key};
    my @output;
    if( ref( $out ) eq "HASH" ){
        for my $val ( sort(keys(%{$out}))){
            print $fh join("=",("##".$key, $out->{$val}))."\n";
        }
        delete $hash->{$key};
        return 1;
    }
    print $fh join("=",("##".$key,$out))."\n";
    delete $hash->{$key};
    return 1;   
}

#utility function for printing the column headers
sub print_and_delete_col_header {
    my $self = shift;
    my $hash = shift;
    my $fh = shift;
    print $fh "#".$hash->{CHROM}."\n";
    delete $hash->{CHROM};
    return 1;   
}

# this will get lines from the two inputs as needed, then reformat and or merge them, and print them to _output_fh
sub process_records {
    my $self = shift;
    my $afh = $self->_vcf_a_fh;
    my $bfh = $self->_vcf_b_fh;
    my $ofh = $self->_output_fh;

    #get the first line from each file a and b
    my $line_a = $afh->getline;
    chomp $line_a;
    my $line_b = $bfh->getline;
    chomp $line_b;
    my ($chrom_a,$pos_a) = split /\t/, $line_a;
    my ($chrom_b,$pos_b) = split /\t/, $line_b;

    my $done=0;

    #loop until one of the files reaches EOF (it sets $done to a or b, denoting which file is done)
    while( ! $done ){
        my $comparison = $self->compare($chrom_a,$pos_a,$chrom_b,$pos_b);
        if($comparison == -1){
            $self->print_record($line_a,'a');
            $done = ($self->get_next($afh,\$line_a,\$chrom_a,\$pos_a)) ? 'a' : 0;
        } elsif ($comparison == 1) {
            $self->print_record($line_b,'b');
            $done = ($self->get_next($bfh,\$line_b,\$chrom_b,\$pos_b)) ? 'b' : 0;
        } else {
            my $merged = $self->merge_records($line_a,$line_b);
            $self->print_merged($merged);
            my $done_a = ($self->get_next($afh,\$line_a,\$chrom_a,\$pos_a)) ? 'a' : 0;
            my $done_b = ($self->get_next($bfh,\$line_b,\$chrom_b,\$pos_b)) ? 'b' : 0;
            $done = $done_a || $done_b;
        }
    }

    #Then dump the remaining records from the non EOF file into the output
    my $remaining_fh = ($done eq 'a') ? $bfh : $afh;
    #Make sure to then switch the filehandle label... 
    $done = ($done eq 'a') ? 'b' : 'a';

    my ($chr,$pos);
    while(my $line = $remaining_fh->getline){
        chomp $line;
        $self->print_record($line,$done);
    }

    return 1;
}

# This pulls the next line from the file handle $fh
# It sets $line, $chrom, and $pos, which should be passed in as refs to scalars,
# and returns 1 if the EOF is encountered.
sub get_next {
    my $self = shift;
    my ( $fh, $line, $chrom,$pos) = @_;
    my $answer = 0;    
    unless($$line = $fh->getline){
        return 1;
    }
    chomp $$line;
    ($$chrom,$$pos) = split /\t/, $$line;
    return $answer;
}

# This prints $line to the _output_fh, after processing the FORMAT fields and SAMPLE fields, 
# adding '.' to any that are specified in _format_fields, but not in $line.
sub print_record {
    my $self = shift;
    my $line = shift;
    my $s = shift;
    my $fh = $self->_output_fh;

    my $source = ($s eq 'a') ? $self->vcf_file_a_source : $self->vcf_file_b_source;
  
    my @cols = split /\t/, $line;
    if($cols[INFO] eq '.'){
        $cols[INFO] = "VC=".$source;
    } else {
        $cols[INFO] .= ";VC=".$source;
    }

    ($cols[FORMAT], $cols[SAMPLE]) = $self->adjust_format_and_sample_string($cols[FORMAT], $cols[SAMPLE]);
    print $fh join("\t",@cols) . "\n";

    return 1;
}

# This method constructs format and sample fields in the proper order, adding dots where fields were missing
sub adjust_format_and_sample_string {
    my $self = shift;
    my $format_in = shift;
    my $sample_in = shift;

    my @formats = split /:/, $format_in;
    my @sample = split /:/, $sample_in;    

    my %form = %{ $self->_format_fields };

    my @output_formats;
    my @output_sample;
    for my $field (0..(scalar(@formats)-1)){
        if(exists($form{$formats[$field]})){
            push @output_formats, $formats[$field];
            push @output_sample, $sample[$field];
            delete $form{$formats[$field]};
        } else {
            die $self->error_message("Could not locte format field definiton for: ".$formats[$field]." in header.");
        }
    }
    unless( ! keys( %form )){
        for my $field (sort(keys(%form))){
            push @output_formats, $field;
            push @output_sample, '.';
            delete $form{$field};
        } 
    }

    my $format_string = join(":",@output_formats);
    my $sample_string = join(":",@output_sample);

    return ($format_string, $sample_string);
}

#This prints a vcf line hash to the _output_fh
sub print_merged {
    my $self = shift;
    my $m = shift;
    my %m = %{ $m };
    my $fh = $self->_output_fh;
    my $line = join("\t",($m{CHROM},$m{POS},$m{ID},$m{REF},$m{ALT},$m{QUAL},$m{FILTER},$m{INFO},$m{FORMAT},$m{SAMPLE}));
    print $fh $line . "\n";

    return 1;
}

#return -1 if $chr_a,$pos_a represents a lower position than $chr_b,$pos_b, 0 if they are the same, and 1 if b is lower
sub compare {
    my $self = shift;
    my ($chr_a,$pos_a,$chr_b,$pos_b) = @_;
    if(($chr_a eq $chr_b) && ($pos_a == $pos_b)){
        return 0;
    }
    if($chr_a eq $chr_b){
        return ($pos_a < $pos_b) ? -1 : 1;
    }
    return ($self->chr_cmp($chr_a,$chr_b)) ? 1 : -1;
}

# return 0 if $chr_a is lower than $chr_b, 1 otherwise
sub chr_cmp {
    my $self = shift;
    my ($chr_a, $chr_b) = @_;
    my @chroms = ($chr_a,$chr_b);
    my @answer = nsort @chroms;
    return ($answer[0] eq $chr_a) ? 0 : 1;
}

#combine two intersecting records
sub merge_records {
    my $self = shift;
    my ($line_a,$line_b) = @_;

    my @cols_a = split /\t/, $line_a;
    my @cols_b = split /\t/, $line_b;
   
    unless(@cols_a == @cols_b){
        die $self->error_message("both records did not have the same number of columns! \n\t\t ".$line_a."\n\t\t".$line_b);
    } 
    
    my %merged;
    my %gt;
    $merged{CHROM} = $cols_a[CHROM];
    $merged{POS} = $cols_a[POS];
    $merged{ID} = $cols_a[ID];
    $merged{REF} = $cols_a[REF];
    my $alt_alleles;
    # We check later on if there is a conflict between the two ALT values... but regardless always trust A
    $merged{ALT} = $cols_a[ALT];

    if($cols_a[QUAL] ne '.'){
        $merged{QUAL} = $cols_a[QUAL];
    } else {
        $merged{QUAL} = $cols_b[QUAL];
    }

    my $a_pass = $cols_a[FILTER] eq 'PASS';
    my $b_pass = $cols_b[FILTER] eq 'PASS';
    if($a_pass && $b_pass){
        $merged{FILTER} = 'PASS';
    } elsif( $a_pass || $b_pass ){
        if($self->intersection){
            $merged{FILTER} = $a_pass ? $cols_b[FILTER] : $cols_a[FILTER];
        } else {
            $merged{FILTER} = 'PASS';
        }
    } else {
        $merged{FILTER} = $cols_a[FILTER];
    }

    my $a_info = $cols_a[INFO] eq '.';
    my $b_info = $cols_b[INFO] eq '.';

    my $source_a = $self->vcf_file_a_source;
    my $source_b = $self->vcf_file_b_source;


    if($a_info && $b_info){
        $merged{INFO} = "VC=".$source_a.",".$source_b;
    } elsif ($a_info || $b_info) {
        $merged{INFO} = $a_info ? $cols_b[INFO] : $cols_a[INFO];
        $merged{INFO} .= ";VC=".$source_a.",".$source_b;
    } else {
        my %info;
        my @a_info = split /\,/, $cols_a[INFO];
        for my $info (@a_info){
            my ($tag,$value) = split /\=/, $info;
            $info{$tag}=$value;
        }
        my @b_info = split /\,/, $cols_b[INFO];
        for my $info (@b_info){
            my ($tag,$value) = split /\=/, $info;
            unless(exists($info{$tag})){
                $info{$tag}=$value; # FIXME this will potentially stomp on values from A... ok for now probably?
            }
        }
        my @info;
        for my $key (sort(keys(%info))){
           push @info, join("=",$key,$info{$key} );
        }
        $merged{INFO} = join(":",@info);
        $merged{INFO} .= ";VC=".$source_a.",".$source_b;
    }

    my @format_keys_a = split ":", $cols_a[FORMAT];
    my @format_values_a= split ":", $cols_a[SAMPLE];
    unless (scalar @format_keys_a == scalar @format_values_a) {
        die $self->error_message("Format keys and values counts do not match. Malformed VCF at lines $line_a, $line_b");
    }
    my $format_a;
    for my $key (@format_keys_a) {
        $format_a->{$key} = shift @format_values_a;
    }

    my @format_keys_b = split ":", $cols_b[FORMAT];
    my @format_values_b= split ":", $cols_b[SAMPLE];
    unless (scalar @format_keys_b == scalar @format_values_b) {
        die $self->error_message("Format keys and values counts do not match. Malformed VCF at lines $line_a, $line_b");
    }
    my $format_b;
    for my $key (@format_keys_b) {
        $format_b->{$key} = shift @format_values_b;
    }

    # FIXME this is somewhat redundant since we have logic that fills in answers from A or . ... we just need to merge in B here... refactor
    my ($format_string, $sample_string) = $self->adjust_format_and_sample_string($cols_a[FORMAT], $cols_a[SAMPLE]);
    my @format_fields = split ":", $format_string;
    my $merged_format_values;
    # Merge the answers from line b in where we had no answer before
    for my $key (@format_fields) {
        if (defined $format_a->{$key} ) {
            $merged_format_values->{$key} = $format_a->{$key};
        } elsif (defined $format_b->{$key} ) {
            $merged_format_values->{$key} = $format_b->{$key};
        } else {
            $merged_format_values->{$key} = ".";
        }
    }

    $self->check_alt_and_gt($line_a, $line_b, $cols_a[ALT], $cols_b[ALT], $format_a->{GT}, $format_b->{GT});

    $merged{FORMAT} = $format_string;
    $merged{SAMPLE} = join (":", map {$merged_format_values->{$_} } @format_fields);

    return \%merged;
}

# This method will check the two inputs alt and gt values and die or warn if they are something worthy of inspection 
sub check_alt_and_gt {
    my $self = shift;
    my ($line_a, $line_b, $alt_a, $alt_b, $gt_a, $gt_b) = @_;

    if ($alt_a eq $alt_b) {
        # If alt and gt match, everything is peachy
        if ($gt_a eq $gt_b) {
            return 1;
        # If Alt is the same but GT does not match, it should be 0/1 and 1/1... warn but proceed
        } else {
            $self->warning_message("ALT values are the same but gt values are different. Trusting line A and proceeding. Lines:\n$line_a\n$line_b");
            return;
        }
    } else {
        # If alt is different but gt is the same... this could be 0/1 0/2, 1/1, 2/2
        if ($gt_a eq $gt_b) {
            $self->warning_message("ALT values are different and GT is the same. Trusting line A and proceeding. Lines:\n$line_a\n$line_b");
            return;
        # Die in the case that there are two homozygous snps that do not agree... this should not happen.
        } elsif ($gt_a eq "1/1" && $gt_b eq "1/1" ) {
            die $self->error_message("Two homozygous snps that do not agree on the genotype. Lines:\n$line_a\n$line_b");
        # If ALT and GT are different, this is 0/1 2/2 or 1/1 1/2
        } else {
            $self->warning_message("ALT and GT values differ.  Trusting line A and proceeding. Lines:\n$line_a\n$line_b");
            return;
        }
    }
}

1;
