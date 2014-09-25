package Genome::Model::Tools::Annotate::Adaptor::PileupToAnnotation;

use strict;
use warnings;
use IO::File;
use File::Temp;

use Genome;

class Genome::Model::Tools::Annotate::Adaptor::PileupToAnnotation{
    is => ['Genome::Model::Tools::Annotate'], 
    has => [
        snp_file => {
            is => 'Path',
            is_input => '1',
            is_optional => '1',
            doc => 'SNP file in SAMTools format (Pileup format).  This file should be sorted by chromosome and start position',
        },
        indel_file => {
            is => 'Path',
            is_input => '1',
            is_optional => '1',
            doc => 'Indel file in SAMTools format (Pileup format) This file should be sorted by chromosome and start position',
        },
        output_file => {
            is => 'Path',
            is_input => '1',
            is_optional => '0',
            doc => 'Output file in the tab delimited annotator format (chromosome start stop reference variant)',
        },
    ],
};

sub help_synopsis{
    return <<EOS
gmt annotate adaptor pileup-to-annotation --indel-file indelFile.indel --snp-file snpFile.snp --output-file outputFile.out
EOS
}

#TODO: Fill in
sub help_detail{
    return <<EOS

EOS
}

sub execute{

    my $self = shift;
    my $temp_snp_path = "/tmp/snp_file_test"; 
    my $temp_indel_path = "/tmp/indel_file_test";

    unless ($self->snp_file or $self->indel_file){
        $self->error_message("snp_file and/or index_file must be defined");
        return undef;
    }

    my @files_to_merge;

    if ($self->snp_file){
        my @headers = qw/ chromosome_name start reference variant /;
        my $svr = Genome::Utility::IO::SeparatedValueReader->create(
            input=> $self->snp_file,
            headers => \@headers,
            separator => "\t",
            is_regex => 1, 
            ignore_extra_columns => 1,
        );
        my $converted_snp_file = IO::File->new($temp_snp_path, 'w');
        while (my $line = $svr->next){
            $line->{'stop'} = $line->{'start'}; #add the stop position
            @headers = qw/ chromosome_name start stop reference variant /; #throw the new stop header in the headers array so we can make the translation
            my $final_line = join("\t", map { $line->{$_} } @headers) . "\n"; #create a translated line
            $converted_snp_file->print($final_line); #write it to file
        }
        push @files_to_merge, $temp_snp_path;
        undef $converted_snp_file;
    }

    if($self->indel_file){
        my @headers = qw/ chromosome_name start junk_a junk_b junk_c junk_d junk_e junk_f allele_a allele_b /;
        my $svr = Genome::Utility::IO::SeparatedValueReader->create(
            input=> $self->indel_file,
            headers => \@headers,
            separator => "\t",
            is_regex => 1, 
            ignore_extra_columns => 1,
        );
        my $converted_indel_file = IO::File->new($temp_indel_path, 'w');
        while (my $line = $svr->next){
            my @variants;
            #get a stop position
            $line->{'allele_a_stop'} = _determine_indel_stop($line, 'a');
            $line->{'allele_b_stop'} = _determine_indel_stop($line, 'b');
            if (defined $line->{'allele_a_stop'} and defined $line->{'allele_b_stop'}){
                my $variant_a = _create_indel_variant($line, 'a'); 
                my $variant_b = _create_indel_variant($line, 'b'); 
                if ($line->{'allele_a_stop'} <= $line->{'allele_b_stop'}){
                    push @variants, $variant_a;
                    push @variants, $variant_b;
                }
                else{
                    push @variants, $variant_b;
                    push @variants, $variant_a;
                }
                
            }elsif (defined $line->{'allele_a_stop'}){
                my $variant = _create_indel_variant($line, 'a'); 
                push @variants, $variant;
            }elsif (defined $line->{'allele_b_stop'}){
                my $variant = _create_indel_variant($line, 'b');
                push @variants, $variant;
            }
            for my $variant (@variants){
                $converted_indel_file->print($variant."\n");
            }
        }
        push @files_to_merge, $temp_indel_path;
        undef $converted_indel_file;
    }

    if ($self->indel_file and $self->snp_file){
        #Merge two, write to output file in sorted fashion
        my ($temporary_snp_path, $temporary_indel_path) = @files_to_merge;
        my $temporary_indel = IO::File->new($temporary_indel_path, 'r');
        my $temporary_snp = IO::File->new($temporary_snp_path, 'r');
        my $output = new IO::File($self->output_file, "w");

        my $indel = <$temporary_indel>;
        my $snp = <$temporary_snp>;

        while (defined($indel) || defined($snp)){
            if(compare_variants($indel, $snp) eq $indel){
               $output->print($indel);
               $indel = <$temporary_indel>;
            }else{
                $output->print($snp);
                $snp = <$temporary_snp>;
            }
        }
        undef $output;
        undef $temporary_indel;
        undef $temporary_snp;
    }
    else {
        #write sorted file to output file
        my $temporary_file_path = pop @files_to_merge;
        my $temporary_file = new IO::File($temporary_file_path, "r");
        my $output = new IO::File($self->output_file, "w");
        while (my $line = <$temporary_file>){
            $output->print($line);
        }
        undef $output;
        undef $temporary_file;
    }
    unlink $temp_indel_path;
    unlink $temp_snp_path;
    return 1;
}

sub _determine_indel_stop{
    my ($indel, $allele_designation) = @_;
    my $allele = $indel->{"allele_".$allele_designation};
    return undef if ($allele eq '*'); #this allele isn't a variant, therefore no stop
    return $indel->{'start'} + 1 if ($allele =~ /\+/); #insertion: return start + 1
    return $indel->{'start'} + length($allele) - 1; #deletion: return start + string length of indel - 1 (the 1 represents the negative sign in the indel string)
}

#Return the variant that comes first.  Compares the chromosome name, start position, and stop position.  If all of these are equal, compare_variants returns the first value
sub compare_variants{
    my ($first, $second) = @_;
    
    #if either parameter is undefined, return the other
    return $first unless defined $second;
    return $second unless defined $first;
    
    my ($first_chrom, $first_start, $first_stop) = split("\t", $first);
    my ($second_chrom, $second_start, $second_stop) = split("\t", $second);
    my $return_value;

    if($first_chrom lt $second_chrom){
        $return_value =  $first;
    }
    elsif ($second_chrom lt $first_chrom){
        $return_value =  $second;
    }
    elsif($first_start < $second_start){
        $return_value =  $first;
    }
    elsif($second_start < $first_start){
        $return_value =  $second;
    }
    elsif($first_stop < $second_stop){
        $return_value = $first
    }
    elsif($second_stop < $first_stop){
        $return_value = $second;
    }
    else{
        $return_value = $first;
    }
    
    return $return_value;
}

sub _create_indel_variant{
    my ($indel, $allele_designation) = @_;
    my $allele_name = "allele_".$allele_designation;
    my $indel_variant;
    if ($indel->{$allele_name} =~ /\+/ ){
        #This is an insertion.  Build an insertion variant
        $indel_variant = join("\t", $indel->{'chromosome_name'}, $indel->{'start'}, $indel->{$allele_name."_stop"}, "-" , substr($indel->{$allele_name}, 1));
    }
    else{
        #This is a deletion.  Build a deletion variant
        $indel_variant = join("\t", $indel->{'chromosome_name'}, $indel->{'start'}, $indel->{$allele_name."_stop"}, substr($indel->{$allele_name}, 1), "-");
    }
    return $indel_variant;
}

1;
