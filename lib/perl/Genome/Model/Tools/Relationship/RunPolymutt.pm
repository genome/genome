package Genome::Model::Tools::Relationship::RunPolymutt;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $DEFAULT_VERSION = '0.11';
use Cwd;
use File::Basename;
use File::Path;
use Genome::Model::Tools::Relationship::RepairVcf 'fix_alt_and_GT_field';

class Genome::Model::Tools::Relationship::RunPolymutt {
    is => 'Command',
    has => [
       fix_alt_and_gt => {
            is => "Boolean",
            default => 0,
            doc => "If set to true, fix all cases where the REF allele is present in the ALT (this currently happens when we set all_sites or roi_file to force genotype)",
            is_optional => 1,
            is_input => 1,
       },
    ],
    has_input => [
        version => {
            is => 'Text',
            default => $DEFAULT_VERSION,
            doc => "Version to use",
        },
        denovo => {
            is=>'Text',
            is_optional=>1,
            default=>0,
        },
        output_vcf => {
            is=>'Text',
            is_output=>1,
        },
        glf_index => {
            is=>'Text',
        },
        dat_file => {
            is=>'Text',
        },
        ped_file => {
            is=>'Text',
        },
        threads => {
            is=>'Text',
            is_optional=>1,
            default=>4,
        },
        bgzip => {
            is_optional=>1,
            default=>1,
            doc=>'set this to 0 if you prefer uncompressed',
        },
        chr2process=> {
            is_optional=>1,
            default=>undef,
        },
       all_sites => {
            is => "Boolean",
            is_optional => 1,
            default => 0,
            doc => "Output calls on all sites (force genotype) in the glf files.",
       },
       roi_file => {
            is => "Path",
            is_optional => 1,
            doc => "Output calls on all sites (force genotype) in this roi file.",
       },
    ],
    has_param => [
        lsf_resource => {
            is => 'Text',
            default => "-R 'span[hosts=1] rusage[mem=1000] -n 4'",
        },
        lsf_queue => {
            is => 'Text',
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}

my %VERSIONS = (
    '0.02' => '/usr/bin/polymutt0.02',
    '0.10' => '/gscmnt/gc6126/info/medseq/launch_cleft_lip_phenotype_correlation/polymutt_binary/polymutt0.10',
    '0.11' => '/gscmnt/gc6126/info/medseq/launch_cleft_lip_phenotype_correlation/polymutt_binary/polymutt0.11',
    '0.13' => '/usr/bin/polymutt0.13',
    # New, experimental version that optionally runs on vcf and can force-genotype sites
    'vcf.0.01' => '/gscmnt/gc6126/info/medseq/launch_cleft_lip_phenotype_correlation/polymutt_binary/polymuttvcf.0.01',
);

sub path_for_version {
    my $class = shift;
    my $version = shift || $DEFAULT_VERSION;

    unless(exists $VERSIONS{$version}) {
        $class->error_message('No path found for polymutt version ' . $version);
        die $class->error_message;
    }

    return $VERSIONS{$version};
}

sub default_version {
    my $class = shift;

    unless(exists $VERSIONS{$DEFAULT_VERSION}) {
        $class->error_message('Default polymutt version (' . $DEFAULT_VERSION . ') is invalid.');
        die $class->error_message;
    }

    return $DEFAULT_VERSION;
}

sub available_versions {
    return keys(%VERSIONS);
}

sub execute {
    $DB::single=1;
    my $self=shift;
    my $polymutt_cmd= $self->path_for_version($self->version);
    my $ped_file = $self->ped_file;
    my $dat_file = $self->dat_file;
    my $glf_index= $self->glf_index;
    my $threads = $self->threads;
    my $output_vcf = $self->output_vcf;
    my ($temp_output) = Genome::Sys->create_temp_file_path();

    if ( ($self->all_sites or $self->roi_file) and not ($self->fix_alt_and_gt) ) {
        die $self->error_message("all_sites or roi_file set without setting fix_alt_and_gt. This will produce junky output. Is this intentional?");
    }

    my $cmd = $polymutt_cmd;
    $cmd .= " -p $ped_file";
    $cmd .= " -d $dat_file";
    $cmd .= " -g $glf_index";
    $cmd .= " --minMapQuality 1";
    $cmd .= " --nthreads $threads";

    # The output param name changed in the latest version
    if ($self->version eq "0.02" or $self->version eq "0.10" or $self->version eq "0.11") {
        $cmd .= " --vcf $temp_output";
    } else {
        $cmd .= " --out_vcf $temp_output";
    }

    if ($self->chr2process) {
        my $chrs = $self->chr2process;
        $cmd .= " --chr2process $chrs";
    }
    if($self->denovo) {
        $cmd .= " --denovo";
    }
    if($self->all_sites) {
        $cmd .= " --all_sites";
    }
    if($self->roi_file) {
        $cmd .= " --pos " . $self->roi_file;
    }

    print "Running: $cmd\n"; 
    my $rv = Genome::Sys->shellcmd(cmd=> $cmd);
    if($rv != 1) {
        return;
    }

    $self->write_fixed_vcf($temp_output, $output_vcf);

    return 1;
}

# FIXME remove this from merge_and_fix_vcfs
# Fix the various problems with the polymutt header to bring it up to vcf spec
sub write_fixed_vcf {
    my ($self, $input_vcf, $output_vcf) = @_;
    my @info_lines_to_add = (qq|##INFO=<ID=DA,Number=1,Type=Integer,Description="De Novo Mutation Allele">|);
    my @format_lines_to_add = (qq|##FORMAT=<ID=DNGL,Number=10,Type=Integer,Description="Denovo Genotype Likelihoods">|);
    push @format_lines_to_add, qq|##FORMAT=<ID=DNGT,Number=1,Type=String,Description="Genotype">|;
    push @format_lines_to_add, qq|##FORMAT=<ID=DNGQ,Number=1,Type=Integer,Description="Genotype Quality">|;
    ####hack to add in header for bingshan's tag that he hasn't added himself
    my @missing_header_lines = $self->missing_header_lines($input_vcf);

    # Put missing lines in the right places
    for my $missing_line (@missing_header_lines) {
        if ($missing_line =~ m/^##FORMAT/) {
            push @format_lines_to_add, $missing_line;
        } elsif ($missing_line =~ m/^##INFO/) {
            push @info_lines_to_add, $missing_line;
        } else {
            die $self->error_message("Could not figure out if this missing line is header or format: $missing_line");
        }
    }

    ######

    my $ifh = Genome::Sys->open_file_for_reading($input_vcf);
    my $fh;
    if ($self->bgzip) {
        $fh = Genome::Sys->open_gzip_file_for_writing($output_vcf);
    } else {
        $fh = Genome::Sys->open_file_for_writing($output_vcf);
    }

    my ($info_printed, $format_printed)=(0,0);
    while(my $line = $ifh->getline) {
        if($line =~m/^#/) {
            if ($line =~ m/ID=GL.*Type=Unsigned Char/) {
                $line =~ s/Unsigned Char/Integer/;
                $line =~ s/Number=10/Number=3/;
            }

            # Fix incorrect data types
            if($line =~m/ID=PS/) {
                $line =~ s/Integer/Float/;
            }
            #Fix spacing before description
            $line =~s/, Description/,Description/;

            # Fix invalid type labels
            $line =~s/,String/,Type=String/;

            # Fix lines without a "number" tag
            if ( ($line =~ m/^##INFO/ or $line =~ m/^##FORMAT/) and not ($line =~ m/Number=/) ) {
                $line =~ s/,/,Number=1,/;
            }

            if($line =~m/#INFO/ && !$info_printed) {
                for my $info_line (@info_lines_to_add) {
                    $fh->print($info_line ."\n");
                }
                $info_printed=1;
            }
            if($line =~m/#FORMAT/ && !$format_printed) {
                for my $format_line (@format_lines_to_add) {
                    $fh->print($format_line ."\n");
                }
                $format_printed=1
            }
            $fh->print($line);
        }
        else {
            if ($self->fix_alt_and_gt) {
                my $fixed_line = $self->fix_variant_line($line);
                $fh->print("$fixed_line\n");
            } else {
                $fh->print($line);
            }
        }
    }
    $fh->close; $ifh->close;

    return 1;
}

sub fix_variant_line {
    my ($self, $line) = @_;
    chomp $line;
    my($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split("\t", $line);
    my ($fixed_alt, $fixed_format, @fixed_samples) = fix_alt_and_GT_field($ref, $alt, $format, @samples);
    my $new_line = join("\t", ($chrom, $pos, $id, $ref, $fixed_alt, $qual, $filter, $info, $fixed_format, @fixed_samples));
    return $new_line;
}

# FIXME copied... call this in Genome::Model::Tools::Relationship::MergeAndFixVcfs or move to a base class
sub missing_header_lines {
    my ($self, $input_vcf) = @_;

    my %possible_tags = (
        AB => qq|##INFO=<ID=AB,Number=1,Type=Float,Description="Allelic Balance">|,
        BA => qq|##INFO=<ID=BA,Number=1,Type=String,Description="Best Alternative Allele">|,
        DQ => qq|##INFO=<ID=DQ,Number=1,Type=Float,Description="De Novo Mutation Quality">|,
        DS => qq|##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage: Defined As the Expected Alternative Allele Count">|,
    );

    my @header_lines_to_add;
    for my $tag (keys %possible_tags) {
        $DB::single=1;
        chomp(my @number_of_tags = `cat $input_vcf | grep $tag`);
        my $has_header=0;
        for my $line (@number_of_tags) {
            if($line=~m/<ID=$tag,/) {
                # Make sure the tag is of the same type...
                if ( ($line =~ m/FORMAT/ and $possible_tags{$tag} =~ m/FORMAT/ ) or ($line =~ m/INFO/ and $possible_tags{$tag} =~ m/INFO/ ) ) {
                    $has_header=1;
                }
            }
        }
        unless($has_header) {
            push @header_lines_to_add, $possible_tags{$tag};
        }
    }

    return @header_lines_to_add;
}

1;
