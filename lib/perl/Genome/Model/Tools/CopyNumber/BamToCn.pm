package Genome::Model::Tools::CopyNumber::BamToCn;

use strict;
use warnings;
use Getopt::Std;
use Statistics::Descriptive;

class Genome::Model::Tools::CopyNumber::BamToCn {
    is => 'Command',
    has_input => [
        aligned_reads_input => {
            is => 'String',
            doc => 'BAM upon which to operate',
        },
        output_file => {
            is => 'String',
            doc => 'File to send output to',
            is_output => 1,
        }
    ],
    has_optional_input => [
        downsampling_ratio => {
            is => 'String',
            doc => 'Downsampling ratio',
            default => '1',
        },
        quality_cutoff => {
            is => 'String',
            doc => 'MAQ mapping quality cutoff',
            default => '0',
        },
        window_size => {
            is => 'String',
            doc => 'Copy Number Window Size',
            default => '10000',
        },
        neutral_ratio => {
            is => 'String',
            default => '0.25',
            doc => 'ratio diverged from median, used to find copy number neutral region',
        },
        chromosome_number => {
            is => 'String',
            doc => 'Chromosome to operate on',
        },
        properly_matched_reads => {
            is => 'Boolean',
            doc => "use only properly matched reads",
        },
    ],
    has_param => [
        lsf_queue => {
            default => "long",
        },
    ],

};


sub help_detail {
    "Output read depth based copy number in a single genome in one or more map/bam files"
}

# This script generates read counts and copy number from a single bam file,
# which can be used for copy number segmentation

sub execute {
    my $self = shift;

    my $version="BAM2CN-0.0.1r2";
    my %opts1;
    my %opts;

    my $options='';

    $opts{'c'} = $self->chromosome_number;
    $opts{'n'} = $self->neutral_ratio;
    $opts{'w'} = $self->window_size;
    $opts{'q'} = $self->quality_cutoff;
    $opts{'r'} = $self->downsampling_ratio;
    $opts{'p'} = $self->properly_matched_reads;
    
    my $output = Genome::Sys->open_file_for_writing($self->output_file);

    ######################## Read In Configuration file #######################
    my $fbam=$self->aligned_reads_input;
    
    ######################## Compute read counts in sliding windows ########################
    my %data;
    my %Coor;
    my %Chrs;

    my $cmd="samtools view $fbam";
    $cmd.=" $opts{c}" if(defined $opts{c});
    open(MAP,"$cmd |") || die "unable to open $fbam\n";
    my ($ppos,$nread_buf);
    my $pchr=0;
    my $idx=-1;
    while(<MAP>){
      chomp;
      my $rr=rand();
      next if($rr>$opts{r});  #down sampling reads
      my @t=split;
      my ($chr,$pos,$mapqual,$dist,$flag);
      $flag=0;
      $mapqual=$t[4];
      $chr=$t[2];
      $pos=$t[3];
      $flag=1 if(!defined $opts{p} || ($t[1]=~/P/ || $t[1] & 0x0002));

      next unless ($mapqual=~/^\d+$/ && $mapqual>=$opts{'q'} && $flag && $chr=~/\S+/ && $pos=~/^\d+$/);
      next if(defined $opts{'c'} && $opts{'c'} ne $chr);

      if($chr ne $pchr){  # a new chromosome
        # reset
        $Chrs{$chr}=1;
        $pchr=$chr;
        $nread_buf=0;
        $ppos=0;
        $idx=-1;
      }
      $nread_buf++;
      if($pos>$ppos+$opts{'w'}){
        do{
          $ppos+=$opts{'w'};
        } until($pos<$ppos+$opts{'w'});

        #register
        $idx++;
        ${$data{$chr}}[$idx]+=$nread_buf;
        ${$Coor{$chr}}[$idx]=$ppos;
        $nread_buf=0;
      }
    }
    close(MAP);


    my @chrs=sort keys %Chrs;
    die "No reads passed QC.\n" if($#chrs<0);

    my $median=Statistics::Descriptive::Full->new();
    foreach my $chr(@chrs){
      my $tmpchr=$chr; $tmpchr=~s/chr//i;
      next unless($tmpchr=~/^\d+/);  #skip non-autosome
      my $chr_median=&Get_Median($data{$chr});
      #print $output "#Chr${chr}_Median:$chr_median\n";
      $median->add_data($chr_median);
    }
    my $md=$median->median();

    print $output "#WholeGenome_Median:$md\n";
    my $num_CN_neutral_pos=0;
    my $NReads_CN_neutral=0;
    foreach my $chr(@chrs){
      next unless (defined $data{$chr});
      my $Nwi=$#{$data{$chr}};

      for(my $i=0;$i<=$Nwi;$i++){
        my $f2x=1;
        next unless (defined ${$data{$chr}}[$i]);

        $f2x=0 if(${$data{$chr}}[$i]<$md*(1-$opts{n}) || ${$data{$chr}}[$i]>$md*(1+$opts{n}));
        next if(! $f2x);
        $num_CN_neutral_pos++;
        $NReads_CN_neutral+=${$data{$chr}}[$i];
      }
    }

    #subtract the normal from the tumor
    my $depth2x=($num_CN_neutral_pos>10)?$NReads_CN_neutral/$num_CN_neutral_pos:$md;
    printf $output "#2xReadCount:%.2f in %d bp\n",$depth2x,$opts{w};
    print $output "CHR\tPOS\tReadCount\tCopyNumber\n";

    foreach my $chr(@chrs){
      next unless (defined $data{$chr});
      my $Nwi=$#{$data{$chr}};

      for(my $i=0;$i<=$Nwi;$i++){
        next unless (defined ${$data{$chr}}[$i]);
        my $cn=${$data{$chr}}[$i]*2/$depth2x;
        my $poschr=${$Coor{$chr}}[$i];
        printf $output "%s\t%d\t%.6f\t%.6f\n",$chr,$poschr,${$data{$chr}}[$i],${$data{$chr}}[$i]*2/$depth2x;
      }
    }
    $output->close;
    return 1;
}
sub Get_Median {
  my $rpole = shift;
  my @pole = @$rpole;
  my $ret;

  @pole=sort {$a<=>$b} @pole;
  if( (@pole % 2) == 1 ) {
    $ret = $pole[((@pole+1) / 2)-1];
  } else {
    $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
  }
  return $ret;
}

1;
