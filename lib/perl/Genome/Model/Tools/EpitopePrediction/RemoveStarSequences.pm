package Genome::Model::Tools::EpitopePrediction::RemoveStarSequences;

use strict;
use warnings;
use Bio::SeqIO;



class Genome::Model::Tools::EpitopePrediction::RemoveStarSequences {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Input Fasta format file',
        },
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file after removing star sequences',
        },
    ],
};


sub help_brief {
    "Outputs a FASTA file after removing *s from the input Fasta sequence",
}



sub execute {
    my $self = shift;
  #  my $input_fh = Genome::Sys->open_file_for_reading($self->input_file);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);

    my $input_filename = $self->input_file ;
    my $output_filename= $self->output_file;

 #   while (my $line = $input_fh->getline) {
#        chomp $line;
 #}

        my $in  = Bio::SeqIO->new(-file => $input_filename ,
                      -format => 'Fasta');

        my $out = Bio::SeqIO->$output_fh(-file => ">".$output_filename ,
                       -format => 'Fasta');

     #     my $out = Bio::SeqIO->new(-fh => $output_fh,
      #                          -format => 'fasta');


        while ( my $seq = $in->next_seq() )
        {
            my $seq_string = $seq->seq;
            print $seq_string."\n";

#               if ( $seq_string =~ /[^A-Z]/ )
     #          {
#                   print "Skipping Sequence:"."\t".$seq_string."\n";

               }
#            else
  #             {
 #                  $out->write_seq($seq);
  #             }
    #    }

1;
}
return 1;
__END__

 my ($temp_fh_name, $temp_name) = Genome::Sys->create_temp_file();

    my $temp_fh      =  Genome::Sys->open_file_for_writing($temp_name);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);




# READ THE SEQ FOF
my $in_file = $ARGV[0];
open (IN, "$in_file") or die 'could not open in file';
while (<IN>) {
    chomp;
    my $parent_dir='/gscmnt/sata141/techd/jhundal/Schreiber_mouse/SchreiberProject/';
    my $input_filename = $parent_dir.$_."_17_mer.fa";
    my $output_filename = $parent_dir."testing/".$_."_17mer_NoStar.fa";
  #  my $netmhc_output_filename = $parent_dir.$_.".netmhc.out";
    print "Input File is".$input_filename."\n";
    print "Output File is".$output_filename."\n";



##########################################################################


    my $in  = Bio::SeqIO->new(-file => $input_filename ,
                      -format => 'Fasta');

    my $out = Bio::SeqIO->new(-file => ">".$output_filename ,
                        -format => 'Fasta');

    while ( my $seq = $in->next_seq() ) {
    my $seq_string = $seq->seq;
       if ( $seq_string =~ /[^A-Z]/ )
       {
           print $seq_string."\n";

       }
       else
       {

       $out->write_seq($seq);
       }
    }
}
__END__
################CALCULATING NUMBER OF LINES IN INPUT FASTA FILE & SPLITTING#################

    my $lines;
    open (FILE, "$output_filename") or die 'Could not open output file';
    $lines++ while (<FILE>);
    close FILE;
    #print "$lines\n";
    print "Number of lines in ".$output_filename." is ".$lines."\n";
               if ($lines <= 10000)
               {
                   #do nothing#
               }
               else
               {
                   my $command = 'split -l 8000 '.$output_filename." ".$output_filename."." ;
                   #print $command."\n";
                   system ( "$command" );
               }

}
