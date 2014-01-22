package utility;
require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(&parseXmlTree);

%EXPORT_TAGS = (
     all => [qw(&parseXmlTree)]
);

use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor qw(:constants);



#######################################################################################################################################################################
#Utility function for grabbing values from hash/array structures                                                                                                      #
#######################################################################################################################################################################
sub parseXmlTree{
  my %args = @_;
  my $ref = $args{'-ref'};
  my $value_name = $args{'-value_name'};

  #If there was only one record for this data type in the XML, a reference to a key-value HASH will be returned, otherwise a reference to an ARRAY will be returned
  #Figure out which is the case ...
  my $ref_test = ref($ref);
  my @values;
  if ($ref_test eq "HASH"){
    my $value = $ref->{$value_name};

    #If you still have an array reference... dereference, join into a string, and push onto an array
    if (ref($value) eq "ARRAY"){
      #print YELLOW, "\nDebug: $value", RESET;
      my @tmp = @{$value};
      $value = join(", ", @tmp);
    }

    if (defined($value)){
      push(@values, $value);
    }else{
      push(@values, "na");
    }
  }else{
    foreach my $x (@{$ref}){
      my $value = $x->{$value_name};
      if (defined($value)){
        push(@values, $value);
      }else{
        push(@values, "na");
      }
    }
  }
  return(\@values);
}







1;



