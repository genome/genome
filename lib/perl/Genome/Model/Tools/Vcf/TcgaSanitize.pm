package Genome::Model::Tools::Vcf::TcgaSanitize;

use strict;
use warnings;
use Genome;
use File::Basename;
use List::MoreUtils qw(firstidx);
use Genome::Utility::Vcf ('get_vcf_header', 'parse_vcf_line', 'deparse_vcf_line', 'get_samples_from_header', 'parse_header_format');

class Genome::Model::Tools::Vcf::TcgaSanitize{
    is => 'Command',
    doc => "Fix up a VCF file for TCGA submission",
    has => [
        output_file => {
            is => 'Text',
            doc => "Output sanitized VCF",
        },
        input_file => {
            is => 'Text',
            doc => "Input vcf to be sanitized",
        },
        package_for_tcga => {
            is => 'Boolean',
            default => 0,
            doc => "Instead of just printing out a sanitized vcf, package and md5 it for TCGA validation and submission",
        },
    ],
};

sub help_synopsis {
    <<'HELP';
Fix up a VCF file for TCGA submission
HELP
}

sub help_detail {
    <<'HELP';
Fix up a VCF file for TCGA submission
HELP
}

sub execute {
    my $self = shift;
    my ($ifh, $header);
    my $ofh;
    # If we are packaging for tcga, we don't want the output file zipped, because it will be later
    if ($self->package_for_tcga) {
        # This is dumb but better than being unclear
        unless ($self->output_file =~ m/\.tar\.gz$/) {
            die $self->error_message("If packaging for TCGA, the file must end in .tar.gz (sorry, blame the tcga validator)");
        }
        $ofh = Genome::Sys->open_file_for_writing($self->output_file);
    } else {
        $ofh = Genome::Sys->open_gzip_file_for_writing($self->output_file);
    }

    ($ifh, $ofh, $header) = $self->fix_header($ofh);
    $self->fix_body($ifh, $ofh, $header);
    $ifh->close; $ofh->close;

    # TODO deploy and use the validator
    if ($self->package_for_tcga) {
        $self->package_file_for_tcga;
    }

    return 1;
}

# Fix all problems with the header
sub fix_header {
    my ($self, $ofh) = @_;
    my ($header, $ifh);

    ($header, $ifh) = get_vcf_header($self->input_file);

    $header = $self->fix_reference_and_add_assembly($header);
    $header = $self->fix_header_source($header);
    $header = $self->fix_header_vcf_process_log($header);
    $header = $self->fix_header_sample($header);

    $ofh->print($header."\n");

    return ($ifh, $ofh, $header);
}

# Fix bad characters in the ##reference line and add the ##assembly line so we can fix odd chromosomes
sub fix_reference_and_add_assembly {
    my ($self, $header) = @_;
    my @lines = split("\n", $header);
    for my $i (0..$#lines) {
        my $line = $lines[$i];
        if ($line =~ m/^##reference/) {
            my ($reference_value) = $line =~ m/##reference=(.*)$/;
            unless ($reference_value) {
                die $self->error_message("Malformed reference line?: $line");
            }

            $reference_value =~ s/[\s|=|,|;]/-/g;
            $lines[$i] = "##reference=$reference_value";
            my $assembly_line = "##assembly=$reference_value";
            splice(@lines, $i+1, 0, $assembly_line);
            $header = join "\n", @lines;

            return $header;
        }
    }

    die $self->error_message("Could not find ##reference line in header:\n$header");
}

# In the header, aggregate all "##source=" lines in into one line
sub fix_header_source {
    my ($self, $header) = @_;
    $header = $self->aggregate_header_lines($header, 'source');
    return $header;
}

# In the header, aggregate all "##vcfProcessLog=" lines in into one line
sub fix_header_vcf_process_log {
    my ($self, $header) = @_;
    $header = $self->aggregate_header_lines($header, 'vcfProcessLog');

    #At this point a line like this exists: <InputVCFSource=<Samtools>>-<InputVCFSource=<Sniper>>-<InputVCFSource=<VarscanSomatic>
    # What we want so we pass validation  : <InputVCFSource=<Samtools-Sniper-VarscanSomatic>> 
    my @lines = split("\n", $header);
    my @index_to_fix;
    for my $i (0..$#lines) {
        # Find the vcf process log line and make the repair
        my $line = $lines[$i];
        if ($line =~ m/^##vcfProcessLog=/) {
            push @index_to_fix, $i;
        }
    }
    for my $i (@index_to_fix) {
        my $line = $lines[$i];
        splice(@lines, $i, 1);
        my @values = $line =~ m/<InputVCFSource=<(\w+)>/g;
        my $new_line = "##vcfProcessLog=<InputVCFSource=<" . join("-", @values) . ">>";
        splice(@lines, $i, 0, $new_line);
    }

    $header = join("\n", @lines);

    return $header;
}

# Given the header and a header line type, aggregate all of the lines of that type onto one line
sub aggregate_header_lines {
    my ($self, $header, $type) = @_;

    my $header_base = "##$type=";

    my @lines = split("\n", $header);
    my $original_line_count = scalar(@lines);
    my @index_to_fix;
    for my $i (0..$#lines) {
        my $line = $lines[$i];
        if ($line =~ m/^$header_base/) {
            push @index_to_fix, $i;
        }
    }

    my @values;
    for my $i (@index_to_fix) {
        # Adjust the index we are looking at, since values are getting spliced out
        # Subtract the first index by 0, second index by 1, third index by 2, etc
        $i -= firstidx { $_ == $i } @index_to_fix;
        my $line = $lines[$i];
        my ($value) = $line =~ /^$header_base(.*)/;
        unless ($value) {
            die $self->error_message("Could not obtain value from line $line");
        }
        push @values, $value;
        splice(@lines, $i, 1);
    }

    my $new_header_line = $header_base . join("-", @values);
    # Put the new line where the first old line was
    splice(@lines, $index_to_fix[0], 0, $new_header_line);

    # Make sure things went alright...
    unless (scalar (@lines) + scalar(@index_to_fix) - 1  == $original_line_count) {
        $self->error_message(scalar (@lines) . " " . scalar(@index_to_fix) . " - 1 " .  $original_line_count);
        die $self->error_message("Errors encountered in aggregate_header_lines for type $type");
    }

    $header = join "\n", @lines;

    return $header;
}

# In the header, ensure there is one "##SAMPLE=" line for every sample column (such as the ones created by joinx with the -D option, one column per sample and detector)
sub fix_header_sample {
    my ($self, $header) = @_;
    return $header;
}

# Iterate through each variant line in the vcf and call methods to fix things that are wrong about each
sub fix_body {
    my ($self, $ifh, $ofh, $header) = @_;

    my @header_lines = split("\n", $header);
    my @samples = get_samples_from_header(\@header_lines);
    my $header_format = parse_header_format(\@header_lines);
    while (my $line = $ifh->getline) {
        $line = $self->fix_null_fields($line, \@samples, $header_format);
        $line = $self->fix_ft($line, \@samples);
        $line = $self->fix_odd_chromosomes($line);
        $ofh->print($line);
    }

    return 1;
}

# Every null sample field "." needs to be fixed so it has explicit null values for every FORMAT field, i.e. ".:.:.:." ...etc
# Also make sure that each null field within this has a number of comma,delimited null values according to the header i.e. for NUMBER=4 ".:.,.,.,.:." ... etc
# Remember, Number=. means one per alt
sub fix_null_fields {
    my ($self, $line, $sample_names, $header_format) = @_;
    my $parsed_line = parse_vcf_line($line, $sample_names);

    # For every sample in the line
    my @samples = keys %{$parsed_line->{"sample"}};
    for my $sample (@samples) {
        # For every format field in that sample
        my @fields = keys %{$parsed_line->{"sample"}->{$sample}};
        for my $field (@fields) {
            # If the value is not defined, put a "." there. If there is already a ".", make sure there are the appropriate number of them based upon the header
            if (not defined $parsed_line->{"sample"}->{$sample}->{$field} or $parsed_line->{"sample"}->{$sample}->{$field} eq ".") {
                my $expected_null_values;
                my $number = $header_format->{$field}->{Number};
                if ($number eq ".") {
                    # If the header says number=., it means one per alt (but it seems like for now tcga is fine with just one null value for Number=.)
                    $expected_null_values = 1;
                    #my @number_of_alts = split(",", $parsed_line->{alt});
                    #if (@number_of_alts) {
                    #    $expected_null_values = scalar(@number_of_alts);
                    #} else {
                    #    $expected_null_values = 1;
                    #}
                } else {
                    unless ($number =~ m/^\d+$/) {
                        die $self->error_message("Number value $number for format field $field is not a number or .?");
                    }
                    $expected_null_values = $number;
                }

                my $null_value = join(",", split("", ("." x $expected_null_values)));
                $parsed_line->{"sample"}->{$sample}->{$field} = $null_value;
            }
        }
    }

    $line = deparse_vcf_line($parsed_line, $sample_names);


    return $line;
}

# In our FT fields, we often have multiple filter failures (IntersectionFailure;SnpFilter). The semicolon or commas are not allowed here. 
sub fix_ft {
    my ($self, $line, $sample_names) = @_;
    my $parsed_line = parse_vcf_line($line, $sample_names);
    
    $parsed_line->{"filter"} =~ s/;/-/g;
    for my $sample (@$sample_names) {
        if (defined $parsed_line->{"sample"}->{$sample}->{"FT"}) {
            $parsed_line->{"sample"}->{$sample}->{"FT"} =~ s/;/-/g;
        }
    }
    $line = deparse_vcf_line($parsed_line, $sample_names);

    return $line;
}

# Any chromosome that is not 1-22,X,Y,MT must be enclosed in < >. 
sub fix_odd_chromosomes {
    my ($self, $line) = @_;
    my ($chromosome, @stuff) = split("\t", $line);
    my @acceptable_chromosomes = (1..22, "X", "Y", "MT");
    unless ( grep{uc($_) eq $chromosome} @acceptable_chromosomes ) {
        $chromosome = "<$chromosome>";
    }
    $line = join("\t", ($chromosome, @stuff) );
    return $line;
}

# Instructions from feiyu:
# First TCGA requires the specific archive file name. The format is: '[center]_[disease].[platform].[archive_type].[batch].[revision].[series].tar.gz'
# mkdir genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0
# cd genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0
# md5sum genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0.vcf>   MANIFEST.txt
#   ..
# tar -cvzf genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0.tar.gz genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0/
# md5sum genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0.tar.gz>   genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0.tar.gz.md5

# Package up the vcf in a tarball with md5s as needed to submit to tcga
sub package_file_for_tcga {
    my $self = shift;

    my ($filename, $directory, $suffix) = fileparse($self->output_file, ("\.tar\.gz"));
    my $working_dir = Genome::Sys->create_directory("$directory/$filename");
    my $working_file_basename = "$filename.vcf";
    my $working_file = $working_dir . "/$working_file_basename";

    Genome::Sys->copy_file($self->output_file, $working_file);
    unlink($self->output_file);

    my $md5 = Genome::Sys->md5sum($working_file);
    my $manifest_file = "$working_dir/MANIFEST.txt";
    my $mfh = Genome::Sys->open_file_for_writing($manifest_file);
    $mfh->print($md5 . " " . $working_file_basename);
    $mfh->close;

    my $tar_cmd = "tar -cvzf " . $self->output_file . " $working_dir";
    unless (Genome::Sys->shellcmd(cmd =>$tar_cmd), output_files => [$self->output_file, $self->output_file.".md5"]) {
        die $self->error_message("Failed to run $tar_cmd");
    }

    my $whole_md5 = Genome::Sys->md5sum($self->output_file);
    my $md5_file = $self->output_file . ".md5";
    my $md5_fh = Genome::Sys->open_file_for_writing($md5_file);
    $md5_fh->print($whole_md5);
    $md5_fh->close;

    return 1;
}

1;
