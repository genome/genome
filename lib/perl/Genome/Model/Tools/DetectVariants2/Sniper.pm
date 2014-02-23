package Genome::Model::Tools::DetectVariants2::Sniper;

use warnings;
use strict;

use Genome;
use Genome::Info::IUB;
use Workflow;

my $DEFAULT_VERSION = '0.7.3';
my $LEGACY_SNIPER_COMMAND = 'bam-somaticsniper';
my $SNIPER_COMMAND = 'bam-somaticsniper1.0.0';

class Genome::Model::Tools::DetectVariants2::Sniper {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces a list of high confidence somatic snps and indels.",
# TODO ... make sure this works without old default snv and indel params default => '-q 1 -Q 15',
    # Make workflow choose 64 bit blades
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>100000] span[hosts=1]',
        },
    ],
};

my %SNIPER_VERSIONS = (
    '0.7'   => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7/'   . $LEGACY_SNIPER_COMMAND,
    '0.7.1' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7.1/' . $LEGACY_SNIPER_COMMAND,
    '0.7.2' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7.2/' . $LEGACY_SNIPER_COMMAND,
    '0.7.3' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7.3/' . $LEGACY_SNIPER_COMMAND,
    '1.0.0' => '/usr/bin/' . $SNIPER_COMMAND,
    '1.0.1' => '/usr/bin/bam-somaticsniper1.0.1',
    '1.0.2' => '/usr/bin/bam-somaticsniper1.0.2',
);

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic sniper --aligned-reads-input tumor.bam --control-aligned-reads-input normal.bam --output-directory sniper
gmt somatic sniper --aligned-reads tumor.bam --control normal.bam --out sniper --quality 25
EOS
}

sub help_detail {                           
    return <<EOS 
    Provide a tumor and normal BAM file and get a list of somatic snps.  
EOS
}

sub _detect_variants {
    my $self = shift;
    $self->debug_message("beginning execute");

    #sniper1.0.0 does not have indel output anymore. We will always
    #run sniper1.0.0 with vcf as output format by giving "-F vcf" via
    #strategy detector parameter
    my $params = $self->params;
    my $indel_output = $self->_indel_staging_output;
    my $snp_output;

    if ($params =~ /\-F\s/) { #vcf output
        $snp_output = $self->_temp_staging_directory .'/'.$self->_snv_base_name.'.raw.vcf';
    }
    else {
        $snp_output = $self->_snv_staging_output;
    }
    
    my $cmd = $self->sniper_path . " " . $params . " -f ".$self->reference_sequence_input." ".$self->aligned_reads_input." ".$self->control_aligned_reads_input ." " . $snp_output . " " . $indel_output;
    my $result = Genome::Sys->shellcmd( cmd=>$cmd, input_files=>[$self->aligned_reads_input,$self->control_aligned_reads_input], output_files=>[$snp_output], skip_if_output_is_present=>0, allow_zero_size_output_files => 1, );

    #For now need convert sniper vcf output to classic output to make
    #downstream process easier.
    $self->_convert_vcf_to_classic($snp_output) if $params =~ /\-F\s/;

    #Manually check for $self->_indel_staging_output as there might not be any indels and shellcmd()
    # chokes unless either all are present or all are empty.
    #(This means shellcmd() can check for the SNPs file on its own and still work given an empty result.)
    #Varied the warning text slightly so this message can be disambiguated from shellcmd() output in future debugging
    unless(-s $self->_indel_staging_output) {
        #Touch the file to make sure it exists
        my $fh = Genome::Sys->open_file_for_writing($self->_indel_staging_output);
        unless ($fh) {
            $self->error_message("failed to touch " . $self->_indel_staging_output . "!: " . Genome::Sys->error_message);
            die;
        }
        $fh->close;
        
        $self->warning_message("ALLOWING zero size output file " . $self->_indel_staging_output);
    }

    $self->debug_message("ending execute");
    return $result; 
}

#MQ gets 1 value and means total mapping quality.
#The BQ and AMQ conversion is tricky. Always look at genotype (0/2,1/2)
#first then convert above to bases (G/C, A/C) then sort bases (CG, AC)
#note the order in vcf BQ,AMQ,BCOUNT always follow the order (A,C,G,T);
#1/2 means you get 2 vars so "37,37" is tumor var bq, "30,37" is tumor
#var mq, "2,3" (A,C from BCOUNT) is tumor var read depth;

#2	179190596	.	G	A,C	.	.	.	GT:IGT:DP:DP4:BCOUNT:GQ:JGQ:VAQ:BQ:MQ:AMQ:SS:SSC	0/2:0/2:5:2,0,3,0:1,2,2,0:34:.:34:38,38:35:37,37:1:.	1/2:1/2:5:0,0,5,0:2,3,0,0:15:.:51:37,37:34:30,37:2:21
#2	179190596	G	M	S	21	15	51	34	34	34	35	5	5	0	0	0	37,37	30,37	2,3	38	37	2	38	37	2


sub _convert_vcf_to_classic {
    my ($self, $vcf_output) = @_;
    my $snv_output = $self->_snv_staging_output;

    my $in_fh  = Genome::Sys->open_file_for_reading($vcf_output) or die "Failed to open $vcf_output for read\n";
    my $out_fh = Genome::Sys->open_file_for_writing($snv_output) or die "Failed to open $snv_output for writing\n";

    while (my $line = $in_fh->getline) {
        next if $line =~ /^\#/;
        chomp $line;
        my @columns = split /\t/, $line;
        my @n_data  = split /\:/, $columns[9];
        my @t_data  = split /\:/, $columns[10];
        my @bases   = split /\,/, $columns[4];

        unshift @bases, $columns[3];
        my $num = 0;
        my %map;
        for my $base (@bases) {
            $map{$num} = $base;
            $num++;
        }

        my @n_gt_bases = map{$map{$_}}(split /\//, $n_data[0]);
        my @t_gt_bases = map{$map{$_}}(split /\//, $t_data[0]);

        my $n_gt_str = join '', sort @n_gt_bases;
        my $t_gt_str = join '', sort @t_gt_bases;
        my $n_gt_iub = Genome::Info::IUB->string_to_iub($n_gt_str);
        my $t_gt_iub = Genome::Info::IUB->string_to_iub($t_gt_str);

        my @n_outs = $self->_resolve_outs($map{0}, \@n_gt_bases, \@n_data);
        my @t_outs = $self->_resolve_outs($map{0}, \@t_gt_bases, \@t_data);

        $out_fh->print(join "\t", $columns[0], $columns[1], $columns[3], $t_gt_iub, $n_gt_iub, $t_data[12], $t_data[5], $t_data[7], $t_data[9], $n_data[5], $n_data[7], $n_data[9], $t_data[2], $n_data[2], @t_outs, @n_outs);
        $out_fh->print("\n");
    }
    $in_fh->close;
    $out_fh->close;

    return 1;
}


sub _resolve_outs {
    my ($self, $ref, $gt_bases, $data) = @_;

    my @bq = split /\,/, $data->[8];
    my @mq = split /\,/, $data->[10];  #AMQ
    my %uniq_base;
    map{$uniq_base{$_} = 1}@$gt_bases;
    my @sort_uniq_bases = sort keys %uniq_base;
    
    my %bq_base;
    map{$bq_base{$sort_uniq_bases[$_]}=$bq[$_]}(0..$#bq);
    
    my %mq_base;
    map{$mq_base{$sort_uniq_bases[$_]}=$mq[$_]}(0..$#bq);

    my @bases = qw(A C G T);  #must be in this sorted order
    my %b_ct;
    my @b_cts = split /\,/, $data->[4];
    map{$b_ct{$bases[$_]} = $b_cts[$_]}(0..3);

    my %gt_type;
    for my $base (@sort_uniq_bases) {
        if ($base eq $ref) {
            $gt_type{ref} = $base;
        }
        else {
            push @{$gt_type{var}}, $base;
        }
    }

    my $bq_ref = defined $gt_type{ref} ? $bq_base{$gt_type{ref}} : 0;
    my $mq_ref = defined $gt_type{ref} ? $mq_base{$gt_type{ref}} : 0;
    my $dp_ref = $b_ct{$ref};

    my (@bq_vars, @mq_vars, @dp_vars);
    for my $var_base (@{$gt_type{var}}) {
        push @bq_vars, $bq_base{$var_base};
        push @mq_vars, $mq_base{$var_base};
        push @dp_vars, $b_ct{$var_base};
    }
    
    my $bq_var = @bq_vars ? join ',', @bq_vars : 0;
    my $mq_var = @mq_vars ? join ',', @mq_vars : 0;
    my $dp_var = @dp_vars ? join ',', @dp_vars : 0;

    return ($bq_ref, $mq_ref, $dp_ref, $bq_var, $mq_var, $dp_var);
}


sub sniper_path {
    my $self = $_[0];
    return $self->path_for_sniper_version($self->version);
}

sub available_sniper_versions {
    my $self = shift;
    return keys %SNIPER_VERSIONS;
}

sub path_for_sniper_version {
    my $class = shift;
    my $version = shift;

    if (defined $SNIPER_VERSIONS{$version}) {
        return $SNIPER_VERSIONS{$version};
    }
    die('No path for bam-somaticsniper version '. $version);
}

sub default_sniper_version {
    die "default bam-somaticsniper version: $DEFAULT_VERSION is not valid" unless $SNIPER_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    if(exists($SNIPER_VERSIONS{$version})){
        return 1;
    }
    return 0;
}

1;
