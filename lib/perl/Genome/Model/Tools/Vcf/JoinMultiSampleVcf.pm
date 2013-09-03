package Genome::Model::Tools::Vcf::JoinMultiSampleVcf;

use strict;
use warnings;
use Genome;

use Sort::Naturally qw(nsort);

use constant {   # Copied from Genome::Model::Tools::Vcf::MultiSampleJoinVcf
    ALT => 4,
    FORMAT => 8,
    SAMPLE => 9,
};

class Genome::Model::Tools::Vcf::JoinMultiSampleVcf {
    is => 'Genome::Model::Tools::Vcf::MultiSampleJoinVcf',
    has => [
        vcf_list => {
            is => 'Text',
            is_input => 1,
            doc => 'Path to a file containing a list of multi-sample-vcfs to merge, one per line, with the source_name following, followed by build_id',
        },
        per_chrom => {
            is => 'Text',
            doc => 'If this is set, the files listed in the vcf_list will have <vcf_filename>_CHR.vcf.gz appended',
            is_input => 1,
            is_optional => 1,
        }
    ],
};

sub _process_input_line {
    my($self, $line, $paths, $handles) = @_;

    my ($path) = split /\s+/,$line;

    my $d = Digest::MD5->new;
    $d->add($path);
    my $key = $d->hexdigest;
    
    if(exists($paths->{$key})){
        die $self->error_message("Already have a record for: ".$path);
    }
    if($self->per_chrom){
        my @stuff = split /\./, $path;
        my $chr = $self->per_chrom;
        $path = join(".",($stuff[0]."_".$chr,$stuff[1],$stuff[2]));
    }

    $paths->{$key}{path} = $path;

    if($self->use_gzip_files){
        $handles->{$key}{handle} = Genome::Sys->open_gzip_file_for_reading($path);
    } else {
        $handles->{$key}{handle} = Genome::Sys->open_file_for_reading($path);
    }

}

sub set_sample_cols {
    my $self = shift;

    my $paths = $self->_vcf_handles;
    my @samples;
    for my $key (keys(%{$paths})){
        push @samples, split /\,/, $paths->{$key}{samples}
    }
    #my @uniq_samples = nsort uniq @samples;
    my %big;
    @big{@samples}=1;
    my @cols;
    for my $key ( nsort keys(%big)){
        push @cols, $key;
    }
    $self->_sample_order(\@cols);
    print "Sample list: \n";

    for (@cols){ print $_ . "\n";}

    my $h = $self->_header;

    $h->{CHROM} = join("\t",($h->{CHROM},@cols));

    return 1;
}

#smartly combine n headers, possibly containing n samples, and stuff a hash of it into $self->_header
sub merge_headers {
    my $self = shift;
    my %header;
    my %filter;
    my %info;
    my %format;
    my $handles = $self->_vcf_handles;
    my @file_handles = map{ $handles->{$_}{handle} } keys(%{ $handles });

    #By using the FILE label instead of the "last" command, the file handles are left at the first non-header line of each file, after processing the headers

#    FILE: for my $fh (@file_handles){
    FILE: for my $handle_key (keys(%{$handles})){
        while(my $line = $handles->{$handle_key}{handle}->getline){
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

                        } elsif($tag =~ m/SAMPLE/){
                            my $key = $data[1];
                            if(exists($header{SAMPLE}{$key})){
                                #unless($header{SAMPLE}{$key} eq $data){
                                #    die $self->error_message("Already have a SAMPLE record in the header that o'er laps: ".$data."   and    ".$header{SAMPLE}{$key});
                                #}
                            } else {
                                $header{SAMPLE}{$key} = $data;
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
                    } elsif ($tag =~ m/SAMPLE/){
                        $key = $data[1];
                        $header{SAMPLE}{$key} = $data;
                    } else {
                        $header{$tag} = $data;
                    }
                }
            } else {
                if( $line =~ m/^#CHROM/){
                    $line =~ s/^#//;
                    my $header = $line;
                    my @header = split /\t/, $line;
                    my @samples;
                    for my $header (@header[SAMPLE..$#header]){
                        push @samples, $header;
                    }

                    $handles->{$handle_key}{samples} = join(",",@samples);

                    unless(exists($header{CHROM})){
                        delete @header[SAMPLE..$#header];
                        $header{CHROM} = join("\t",@header);
                    }
                    next FILE;
                } else {
                    die $self->error_message("Parsed header but did not find column headers!");
                }
            }
        }
    }
   
    $self->_header(\%header);
    return 1;
}

sub add_sample_header_tags {
    return 1;
}


sub _print_processed_records {
    my($self, $answer, $lines) = @_;
    $self->_print_multiple_processed_records($answer, $lines);
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
        unless($line = $handles->{$key}{handle}->getline){
            $handles->{$key}{handle}->close;
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


sub _inputs_for_merge_records {
    my($self, $answer, $lines) = @_;

    my %inputs;

    for my $key (@{$answer}){
        my @cols = split /\t/, $lines->{$key}{line};
        my $last_idx = scalar(@cols)-1;
        my $sample_string = join("\t",@cols[SAMPLE..$last_idx]);
        my @new_cols = @cols[0..SAMPLE];
        $new_cols[SAMPLE] = $sample_string;
        $inputs{$key} = \@new_cols;
    }

    return \%inputs;
}

sub _merge_records_info_for_samples {} # do nothing for this part

sub _merge_records_final_processing {
    my($self, $info, $merged, $inputs, $samples) = @_;

    my @finfo;
    for my $key (sort(keys(%$info))){
        push @finfo, $key."=".$info->{$key};
    }

    $merged->{INFO} = join(";",@finfo);
    $merged->{SAMPLE} = $self->multi_sample_fields_to_string($merged,$inputs);
}

sub multi_sample_fields_to_string {
    my $self = shift;

    my $merged = shift;
    my $inputs = shift;

    my $handles = $self->_vcf_handles;
    my @format_fields = split /:/, $merged->{FORMAT};
    my %format_fields;
    @format_fields{@format_fields} = 1;
    my %format;
    for my $key ( keys( %{ $inputs } )){
        my @samples = split /\t/, $inputs->{$key}->[SAMPLE];
        my @sample_names = split /\,/, $handles->{$key}{samples}; 
        my @per_sample_data;
        for my $idx (0..$#samples){
            #my %format_fields;
            if( $samples[$idx] ne '.' ){
                @format_fields{@format_fields} = 1;
                my @tags = split /\:/, $inputs->{$key}->[FORMAT];
                my @vals = split /\:/, $samples[$idx];

                my @sample_fields;
                for my $num (0..(scalar(@tags)-1)){
                    if(exists($format_fields{$tags[$num]})){
                        $format{$samples[$idx]}{$tags[$num]} = $vals[$num];
                        delete $format_fields{$tags[$num]};
                    } else {
                        die $self->error_message("Found a format tag that was not defined in the header: ".$tags[$num]." for ".$samples[$idx]);
                    }
                }
                unless( ! keys( %format_fields )){
                    for my $field (sort(keys(%format_fields))){
                        $format{$samples[$idx]}{$field} = ".";
                        delete $format_fields{$field};
                    } 
                }
                $format{$samples[$idx]}{GT} = Genome::Model::Tools::Vcf::Convert::Base->regenerate_gt( 
                    $merged->{REF},               #Reference Base(s)
                    $inputs->{$key}->[ALT],    #Original ALT
                    $format{$samples[$idx]}{GT},       #Original Genotype
                    $merged->{ALT},               #New ALT
                );
                my @answer;
                for my $tag (@format_fields){
                    push @answer, $format{$samples[$idx]}{$tag};
                }
                push @per_sample_data,join(":",@answer);
            } else {
                push @per_sample_data, '.';

                #TODO  For later pipeline integration, add an output to mark all "dots"
                # to allow for quick bamreadcounts to backfill once this is done
            }
        }   
        $inputs->{$key}->[SAMPLE] = join("\t",@per_sample_data);
    }

    my @gt_fields;
    my @sources;
    for my $key (keys(%{$inputs})){
        my @per_sample_fields = split /\t/,$inputs->{$key}->[SAMPLE];
        push @gt_fields, @per_sample_fields;
        my @ordered_sample_names = split /\,/, $handles->{$key}{samples};
        push @sources, @ordered_sample_names;
    }

    return $self->get_samples_string(\@gt_fields,\@sources);
}

1;
