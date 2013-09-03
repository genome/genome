package Genome::Model::Tools::Vcf::MultiSampleJoinVcf;

use strict;
use warnings;
use Genome;
use File::stat;
use File::Basename;
use DateTime;
use POSIX;
use Sort::Naturally;
use List::MoreUtils qw/ uniq /;

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

class Genome::Model::Tools::Vcf::MultiSampleJoinVcf {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "Output merged VCF",
        },
        vcf_list => {
            is => 'Text',
            is_input => 1,
            doc => 'Path to a file containing a list of vcfs to merge, one per line, with the source_name following',
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
        _vcf_handles => {
            doc => 'hash ref, keyed by source_name, to filehandles for all vcfs',
        },
        _vcf_list => {
            doc => 'hash ref, keyed by source_name, of all paths',
        },
        _sample_order => {
            doc => 'sorted list of sample column names',
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
        _format_string => {
            doc => 'The string to be used for each format field, perhaps this can vary, but for now it will be the same for each record',
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

    $self->process_input_list;

    $self->merge_headers; 

    $self->add_sample_header_tags;

    $self->set_sample_cols;

    $self->process_format;

    $self->set_format_string;

    my $output = $self->output_file;

    if($self->use_gzip_files){
        $self->_output_fh(Genome::Sys->open_gzip_file_for_writing($output));
    } else {
        $self->_output_fh(Genome::Sys->open_file_for_writing($output));
    }

    $self->print_header;

    $self->process_records;

    $self->close_handles;

    $self->_output_fh->close;

    return 1;
}


sub _process_input_line {
    my($self, $line, $paths, $handles) = @_;

    my ($path,$source_name,$build_id) = split /\s+/,$line;
    if(exists($paths->{$source_name})){
        die $self->error_message("Already have a record for: ".$source_name);
    }
    $paths->{$source_name}{path} = $path;
    if(defined($build_id)){
        $paths->{$source_name}{build_id} = $build_id;
    }
    if($self->use_gzip_files){
        $handles->{$source_name} = Genome::Sys->open_gzip_file_for_reading($path);
    } else {
        $handles->{$source_name} = Genome::Sys->open_file_for_reading($path);
    }
}
 
sub process_input_list {
    my $self = shift;

    my $file = $self->vcf_list;
    my $fh = Genome::Sys->open_file_for_reading($file);

    my %paths;
    my %handles;
    while(my $line = $fh->getline){
        chomp $line;
        $self->_process_input_line($line, \%paths, \%handles);
    }
    $self->_vcf_list(\%paths);
    $self->_vcf_handles(\%handles);

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

sub set_format_string {
    my $self = shift;
    my %form = %{ $self->_format_fields };
    my @answer;
    if(exists($form{GT})){
        push @answer, "GT";
        delete $form{GT};
    }
    if(exists($form{GQ})){
        push @answer, "GQ";
        delete $form{GQ};
    }
    if(exists($form{DP})){
        push @answer, "DP";
        delete $form{DP};
    }

    push @answer, nsort(keys(%form));
    $self->_format_string(join(":",@answer));
}


#### >>>>>>>>>>>>>>>>>>
sub set_sample_cols {
    my $self = shift;

    my $paths = $self->_vcf_list;
    my @samples = nsort keys(%{$paths});
    $self->_sample_order(\@samples);

    my $h = $self->_header;

    $h->{CHROM} = join("\t",($h->{CHROM},@samples));

    return 1;
}
### ===============

sub close_handles {
    my $self = shift;

    my $h = $self->_vcf_handles;

    for my $key (sort(keys(%{$h}))){
        $h->{$key}->close;
    }
    return 1;
}


#smartly combine two headers and stuff a hash of it into $self->_header
sub merge_headers {
    my $self = shift;
    my %header;
    my %filter;
    my %info;
    my %format;
    my $handles = $self->_vcf_handles;
    my @file_handles = map{ $handles->{$_} } keys(%{ $handles });
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
                            unless($header{$tag} =~ m/$data/){
                                $header{$tag} .= ",".$data;
                            }
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
                    my @header = split /\t/, $line;
                    delete $header[-1];
                    unless(exists($header{CHROM})){
                        #unless($header eq $header{CHROM}){
                        #    die $self->error_message("Column headers are not equal!");
                        #}
                    #} else {
                        $header{CHROM} = join("\t",@header);
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

sub add_sample_header_tags {
    my $self = shift;
    my $header = $self->_header;
    my $files = $self->_vcf_list;
    for my $sample (nsort(keys(%{$files}))){
        my $line = "<ID=".$sample.",Build_ID=".$files->{$sample}{build_id}.",Description=\"TGI build_id to sample map\">";
        my $key = $sample.",Number";
        $header->{SAMPLE}{$key} = $line;
    }
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
    $self->print_and_delete_from_hash($h,"SAMPLE",$fh);
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

sub opening_lines {
    my $self = shift;
    my $lines = shift;
    my $handles = $self->_vcf_handles;
    for my $key (sort(keys(%{$handles}))){
        my $line = $handles->{$key}->getline;
        chomp $line;
        my ($chr,$pos) = split /\s+/, $line;
        $lines->{$key}{line} = $line;
        $lines->{$key}{chr} = $chr;
        $lines->{$key}{pos} = $pos;
    }
    return 1;
}

# this will get lines from the two inputs as needed, then reformat and or merge them, and print them to _output_fh
sub process_records {
    my $self = shift;
    my $ofh = $self->_output_fh;

    #my $handles
    my %lines;
    #get the first line from each file
    $self->opening_lines(\%lines);

    #loop until one of the files reaches EOF (it sets $done to a or b, denoting which file is done)
    while( keys(%lines)  ){
        my @answer;
        my @keys = sort(keys(%lines));
        push @answer, shift @keys;
        for my $key (@keys){
            my $cmp = $self->compare(
                $lines{$key}{chr},
                $lines{$key}{pos},
                $lines{$answer[0]}{chr},
                $lines{$answer[0]}{pos},
            );
            if($cmp == 0){
                push @answer, $key;
            } elsif( $cmp == -1){
                @answer = ();
                push @answer, $key;
            }
        }
        $self->_print_processed_records(\@answer, \%lines);

        $self->get_lines(\%lines,\@answer);
    }

    return 1;
}


sub _print_processed_records {
    my($self, $answer, $lines) = @_;

    if (@$answer == 1) {
        $self->print_record($lines->{$answer->[0]}{line}, $answer->[0]);
    } else {
        $self->_print_multiple_processed_records($answer, $lines);
    }
}

sub _print_multiple_processed_records {
    my($self, $answer, $lines) = @_;

    my $line = $self->merge_records($answer,$lines);
    $self->print_merged($line, $answer);
}



# This function takes a hash of sample names, which contain hashes of line, chrom, and pos
# and a list of stale samples, which means the line was used in the output and a new line
# needs to be drawn from the file handle associated with that sample
sub get_lines {
    my $self = shift;
    my $lines = shift;
    my $stale = shift;    

    my $handles = $self->_vcf_handles;

    for my $key (@{$stale}){
        my $line;

        # If the getline doesn't return data (EOF), remove that sample's record
        # from the incoming hash ref, which signals the calling function that
        # the file contains no more data.
        unless($line = $handles->{$key}->getline){
            $handles->{$key}->close;
            delete $handles->{$key};
            delete $lines->{$key};
            next;
        }
        chomp $line;
        my ($chr,$pos) = split /\s+/,$line;
        $lines->{$key}{line} = $line;
        $lines->{$key}{chr} = $chr;
        $lines->{$key}{pos} = $pos;    
    }

    return 1;
}

# This prints $line to the _output_fh, after processing the FORMAT fields and SAMPLE fields, 
# adding '.' to any that are specified in _format_fields, but not in $line.
sub print_record {
    my $self = shift;
    my $line = shift;
    my $source = shift;
    my $fh = $self->_output_fh;

    my @cols = split /\t/, $line;

    $cols[FORMAT] = $self->_format_string;
    $cols[SAMPLE] = $self->adjust_sample_string($cols[FORMAT], $cols[SAMPLE], $source);

    print $fh join("\t",@cols) . "\n";

    return 1;
}

#This prints a vcf line hash to the _output_fh
sub print_merged {
    my $self = shift;
    my $m = shift;
    my $fh = $self->_output_fh;
    my $line = join("\t",($m->{CHROM},$m->{POS},$m->{ID},$m->{REF},$m->{ALT},$m->{QUAL},$m->{FILTER},$m->{INFO},$m->{FORMAT},$m->{SAMPLE}));
    print $fh $line . "\n";

    return 1;
}

sub _inputs_for_merge_records {
    my($self, $answer, $lines) = @_;

    my %inputs;
    for my $key (@{$answer}){
        my @cols = split /\t/, $lines->{$key}{line};
        $inputs{$key} = \@cols;
    }
    return \%inputs;
}

#combine two or more intersecting records
sub merge_records {
    my $self = shift;
    my $answer = shift;
    my $lines = shift;

    my $inputs = $self->_inputs_for_merge_records($answer, $lines);

    my $num_cols=0;
    for my $key (sort(keys(%$inputs))){
        my $nc = scalar( @{$inputs->{$key}});
        if(! $num_cols){
            $num_cols = $nc;
        } else {
            if($nc != $num_cols){
                die $self->error_message("Had differing number of columns in incoming data!".Data::Dumper::Dumper($answer));
            }
        }
    }

    my @samples = nsort(keys(%$inputs));
    my $a = $samples[0];
 
    my %merged;
    $merged{CHROM} = $inputs->{$a}->[CHROM];
    $merged{POS} = $inputs->{$a}->[POS];
    $merged{ID} = $inputs->{$a}->[ID];

    for my $key (@samples){
        unless( $inputs->{$key}->[ID] eq $merged{ID} ){
            die $self->error_message("Found disagreement on ID field at ".$merged{CHROM}." ".$merged{POS});
        }
    }

    $merged{REF} = $inputs->{$a}->[REF];

    for my $key (@samples){
        unless( $inputs->{$key}->[REF] eq $merged{REF} ){
            die $self->error_message("Found disagreement on REF field at ".$merged{CHROM}." ".$merged{POS});
        }
    }

    my  %alts;
    for my $key (nsort(keys(%$inputs))){
        my @alts = split /\,/, $inputs->{$key}->[ALT];
        for my $alt (@alts){
            $alts{$alt} =1;
        }
    }
    $merged{ALT} = join(",", sort(keys(%alts)));

    #FIXME Find a better approach than simply dropping the quality score... some way of merging them?
    $merged{QUAL} = '.';

    $merged{FILTER} = $inputs->{$a}->[FILTER];

    for my $s (@samples){
        my $merged_pass = $merged{FILTER} eq 'PASS';
        my $pass = $inputs->{$s}->[FILTER] eq 'PASS';
        if( $merged_pass xor $pass ){
            if($self->intersection){
                $merged{FILTER} = $merged_pass ? $inputs->{$s}->[FILTER] : $merged{FILTER};
            } else {
                $merged{FILTER} = 'PASS';
            }
        }
    }

    $merged{INFO} = $inputs->{$a}->[INFO];

    my %info;

    for my $sample (@samples){

        my @info = split /\;/, $inputs->{$sample}->[INFO];
        my %tags;
        for my $record (@info){
            my ($t,$v) =  split /\=/, $record;
            unless(exists($info{$t})){
                $info{$t} = $v;
                next;
            }
            $tags{$t} = $v;
        }

        if(exists($tags{VC})&& exists($info{VC})){
            unless($tags{VC} eq $info{VC}){
                my @info_vals = split /\,/, $info{VC};
                my @tag_vals = split /\,/, $tags{VC};        
                my @detectors = uniq (@info_vals,@tag_vals);
                $info{VC} = join(",",@detectors);
            }
        }
        $self->_merge_records_info_for_samples(\%info, \%merged)
    }

    $self->_merge_records_final_processing(\%info, \%merged, $inputs, \@samples);
    $merged{FORMAT} = $self->_format_string;

    return \%merged;
}

sub _merge_records_info_for_samples {
    my($self, $info, $merged) = @_;

    my @finfo;
    for my $key (sort(keys(%$info))){
        push @finfo, $key."=".$info->{$key};
    }
    $merged->{INFO} = join(";",@finfo);
}

sub _merge_records_final_processing {
    my($self, $info, $merged, $inputs, $samples) = @_;

    my @format_fields = split /:/, $merged->{FORMAT};
    my %format_fields;
    @format_fields{@format_fields} = 1;
    my %format;

    for my $sample (@$samples){
        my %format_fields;
        @format_fields{@format_fields} = 1;
        my @tags = split /\:/, $inputs->{$sample}->[FORMAT];
        my @vals = split /\:/, $inputs->{$sample}->[SAMPLE];

        my @sample_fields;
        for my $num (0..(scalar(@tags)-1)){
            if(exists($format_fields{$tags[$num]})){
                $format{$sample}{$tags[$num]} = $vals[$num];
                delete $format_fields{$tags[$num]};
            } else {
                die $self->error_message("Found a format tag that was not defined in the header: ".$tags[$num]." for ".$sample);
            }
        }
        unless( ! keys( %format_fields )){
            for my $field (sort(keys(%format_fields))){
                $format{$sample}{$field} = ".";
                delete $format_fields{$field};
            } 
        }
        $format{$sample}{GT} = Genome::Model::Tools::Vcf::Convert::Base->regenerate_gt( 
            $merged->{REF},               #Reference Base(s)
            $inputs->{$sample}->[ALT],    #Original ALT
            $format{$sample}{GT},       #Original Genotype
            $merged->{ALT},               #New ALT
        );
        my @answer;
        for my $tag (@format_fields){
            push @answer, $format{$sample}{$tag};
        }
        $inputs->{$sample}->[SAMPLE] = join(":",@answer);
    }

    my @gt_fields;
    my @sources;
    for my $sample (@$samples){
        push @gt_fields, $inputs->{$sample}->[SAMPLE];
        push @sources, $sample;
    }

    $merged->{SAMPLE} = $self->get_samples_string(\@gt_fields,\@sources);
}

# This method constructs format and sample fields in the proper order, adding dots where fields were missing
sub adjust_sample_string {
    my $self = shift;
    my $format_in = shift;
    my $sample_in = shift;
    my $source = shift;

    my @formats = split /:/, $format_in;
    my @sample = split /:/, $sample_in;    

    my %formats;

    my @format_fields = split /:/, $self->_format_string;
    my %form;
    @form{@format_fields}=1;

    my @output_formats;
    my @output_sample;
    for my $field (0..(scalar(@formats)-1)){
        if(exists($form{$formats[$field]})){
            $formats{$formats[$field]} = $sample[$field];
            delete $form{$formats[$field]};
        } else {
            die $self->error_message("Could not locte format field definiton for: ".$formats[$field]." in header.");
        }
    }
    unless( ! keys( %form )){
        for my $field (sort(keys(%form))){
            $formats{$field} = '.';
        }
    }
    for my $field (@format_fields){
        push @output_sample, $formats{$field};
    }

    my $sample_string = scalar(@output_sample)>1 ? join(":",@output_sample) : $output_sample[0];

    my @string = ($sample_string);
    my @source = ($source);
    my $full_samples_string = $self->get_samples_string(\@string,\@source);
    return $full_samples_string;
}

sub get_samples_string {
    my $self= shift;
    my $strings = shift;
    my $sources = shift;
    my %ss;
    for my $num (0..(scalar(@{$sources})-1)){
        $ss{$sources->[$num]}=$strings->[$num];
    }

    my $list = $self->_sample_order;
    my @answer;
    for my $s (@{$list}){
        if(exists($ss{$s})){
            push @answer, $ss{$s};
            delete $ss{$s};
        } else {
            push @answer, '.';
        }
    }
    return join("\t",@answer);
}

1;
