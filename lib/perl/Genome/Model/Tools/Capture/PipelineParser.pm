#####################################################################################
# perl module to be used to process capture csv file format.
#####################################################################################

package Genome::Model::Tools::Capture::PipelineParser;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use warnings;
use strict;
use MG::IO::Parse::CSV;
use FileHandle;

#our @ISA = qw( MG::Transform::Process );

#use Data::Dumper;
##__SET THE STYLE OF OUTPUT WHEN USING THE SCRIPT IN "CHECK" MODE
#   $Data::Dumper::Indent = 1;

#########
#  NEW  #
#########

sub new {
   my ($class, %arg) = @_;
   my $self = {
      _source => $arg{source} || 'pipeline_parser',
      _processor => $arg{processor} || {} ,
   };
   bless($self, $class || ref($class));
   return $self;
}

###########
#  PARSE  #
###########

sub Parse {
   my ($self,$fh,$file,%args) = @_;

   # Setup this class to be called (Process method) for secondary processing.
   my $parser;
   $args{version} ||= 0;

#PARSE FORMAT STYLE THAT COMES DIRECTLY FROM GENOME TOOLS 'FROM_GT'
######   header_fields => 'chromosome_name,start,stop,reference,variant,type,gene_name,transcript_name,transcript_source,transcript_version,strand,transcript_status,trv_type,c_position,amino_acid_change,ucsc_cons,domain',

    $parser = MG::IO::Parse::CSV->new(
         source => $self->{_source},
         processor => $self,
         keyfields => 'HUGO_SYMBOL:file_line_num',
         header_translation => {
            'chromosome_name' => 'CHROMOSOME',
            'start' => 'START_POSITION',
            'stop' => 'END_POSITION',
            'reference' => 'REFERENCE_ALLELE',
            'variant' => 'TUMOR_SEQ_ALLELE1',
            'type' => 'VARIANT_TYPE',
            'gene_name' => 'HUGO_SYMBOL',
            'transcript_name' => 'TRANSCRIPT',
	    'transcript_source' => 'SOURCE',
	    'transcript_version' => 'GENOME',
            'strand' => 'STRAND',
            'transcript_status' => 'MUTATION_STATUS',
	    'trv_type' => 'TYPE',
	    'c_position' => 'C_POSITION',
            'amino_acid_change' => 'AA_CHANGE',
	    'ucsc_cons' => 'UCSC',
	    'domain' => 'DOMAIN',
         }
      );
   $self->{_no_process} = $args{no_process};
   $args{no_process} = 0;
   $self->{_check} = $args{check};
   $args{check} = 0;
   return $parser->Parse($fh,$file,%args);
}

#############
#  PROCESS  #
#############

sub Process {
   my ($self,$input,%args) = @_;
   my $processor = $self->{_processor};
   my $check = $self->{_check};
   my $do_process = $self->{_no_process};
   my ($output) = {};

   foreach my $hugo (keys (%{$input})) {
      foreach my $line_num (keys (%{$input->{$hugo}})) {
         $output->{$hugo}->{$line_num} = $input->{$hugo}->{$line_num};
      }
   }
	

   return ($output);
}

=head1 NAME

 Genome::Model::Tools::Capture::Pipeline_Parser -- process csv input for Mutation file format.

=head1 SYNOPSIS

 use Genome::Model::Tools::Capture::Pipeline_Parser;

 # Create a new parser object and process
#__SET SOME PARSING PARAMETERS -- UNSURE OF MEANING OF ORIGINAL COMMENTS (MCW)
   my %parse_args = (

   #__PROCESS EVERYTHING INTO A SINGLE STRUCTURE (AND THEN PROCESSED)
      'all' => 1,

#   #__HAVE EVERYTHING PROCESSED INTO A SINGLE STRUCTURE (AND NOT PROCESSED)
#      'no_process' => 1,
   );

my $fh1 = new FileHandle;
my $fh2 = new FileHandle;
# open "file.annotated"
   unless ($fh1->open (qq{$annotation_file})) {
      die "Could not open mutation project file '$annotation_file' for reading";
   }
# output file MAF
   unless (open($fh2,">$out_file")) {
      die "Could not open mutation project file '$out_file' for reading";
   }

my $parser = Genome::Model::Tools::Capture::PipelineParser->new();
my $annotation = $parser->Parse ($fh1, $annotation_file, %parse_args);
foreach my $hugo (keys %{$annotation}) {
# gene name = $hugo
    foreach my $line_num (keys %{$annotation->{$hugo}}) {
        print STDOUT ".";   #report that we are starting a sample (For commandline user feedback)

        my ($line, $aa_change,$transcript,$mstatus,$Variant_Type,$Chromosome,$Start_position,$End_position,$Reference_Allele,$Tumor_Seq_Allele1,$source,$genome,$strand,$trv_type,$c_position,$ucsc_cons,$domain) =
            (
                $annotation->{$hugo}->{$line_num}->{file_line},
                $annotation->{$hugo}->{$line_num}->{AA_CHANGE},
                $annotation->{$hugo}->{$line_num}->{TRANSCRIPT},
                $annotation->{$hugo}->{$line_num}->{MUTATION_STATUS},
                $annotation->{$hugo}->{$line_num}->{VARIANT_TYPE},
                $annotation->{$hugo}->{$line_num}->{CHROMOSOME},
                $annotation->{$hugo}->{$line_num}->{START_POSITION},
                $annotation->{$hugo}->{$line_num}->{END_POSITION},
                $annotation->{$hugo}->{$line_num}->{REFERENCE_ALLELE},
                $annotation->{$hugo}->{$line_num}->{TUMOR_SEQ_ALLELE1},
		$annotation->{$hugo}->{$line_num}->{SOURCE},
		$annotation->{$hugo}->{$line_num}->{GENOME},
		$annotation->{$hugo}->{$line_num}->{STRAND},
		$annotation->{$hugo}->{$line_num}->{TYPE},
		$annotation->{$hugo}->{$line_num}->{C_POSITION},
		$annotation->{$hugo}->{$line_num}->{UCSC},
		$annotation->{$hugo}->{$line_num}->{DOMAIN},
            );
	$fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$Variant_Type\t$source\t$genome\t$strand\t$trv_type\t$c_position\t$ucsc_cons\t$domain\n");
    }
}
exit 1;

=head1 DESCRIPTION

This module is for parsing the annotation output, using the csv parser at MG::IO::Parse::CSV.

=head1 Constructor and Initialization

=item (object) new (arguments)

=head2 arguments (optional)

=head2 Methods

=item (object) Parse (arguments)

=head2 arguments

=item data structure record

 Hash reference to data structure to process.

=head2 EXPORT

None by default.

=head1 SEE ALSO

=head1 FILES

None

=head1 BUGS

I'm sure that you'll find some. Let me know.

=head1 AUTHOR

William Schierding, E<lt>wschierd@watson.wustl.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by William Schierding.  All Rights Reserved.

=cut

__END__

1;

# $Header$
