# review gsanders
# This seems to be dated... it still uses finfo::SVR instead of Genome::Utility::IO::SVR so at the least update that.
# This seems to, in general, unnecessarily wrap the separated value reader. 
# This module is still used in manual review stuff but perhaps can be dumped in favor of using G:U:IO:SVR directly

package Genome::Utility::VariantReviewListReader;

use strict;
use warnings;

use Genome;

sub new{
    my ($class, $list, $separation_character) = @_;
    $separation_character ||= "|";
    my $separated_value_reader = Genome::Utility::IO::SeparatedValueReader->new(
        input => $list,
        separator => $separation_character,
    );
    die "Can't get SeparatedValueReader for $list!" unless $separated_value_reader;
    my $self = bless({separated_value_reader=>$separated_value_reader, separation_character=>$separation_character}, $class);
    return $self;
}

sub list_columns{
    my ($self, $column_name) = @_;;
    my @columns = (qw/
        chromosome
        begin_position
        end_position
        variant_type
        variant_length
        delete_sequence
        insert_sequence
        genes
        supporting_samples
        supporting_dbs
        finisher_manual_review
        pass_manual_review
        finisher_3730_review
        manual_genotype_normal
        manual_genotype_tumor
        manual_genotype_relapse
        somatic_status
        notes
        /);
    return @columns unless $column_name;
    my $counter=0;
    while ($columns[$counter] ne $column_name){
        $counter++;
    }
    return $counter;
}

sub db_columns{
    #BUILD_ID, REVIEW_TYPE, and DATE are read from the review info file.
    my @columns = ( 
#        'id',
#        'detail_id',
#        'build_id',,
#        'dump_date',
#        'review_type',
    #read only, from VARIANT_REVIEW_DETAIL        
        'chromosome',
        'position',
        'variant_type',
        'delete_sequence',
        'insert_sequence_allele1',
        'insert_sequence_allele2',
        'genes',
        'supporting_samples',
        'supporting_dbs',
    #columns are editable, and used to update SNV_MANUAL_REVIEW
        'notes',
        'genotype_iub_code',#not null
        'pass_manual_review',#not null
        'manual_genotype_iub_normal',#not null
        'manual_genotype_iub_tumor',#not null
        'manual_genotype_iub_relapse',#not null
        'somatic_status',
        'data_needed',
        
        );
    return @columns;
}

sub db_old_columns{
    my @columns = ( qw/
        chromosome
        begin_position
        end_position
        variant_type
        variant_length
        delete_sequence
        insert_sequence_allele1
        insert_sequence_allele2
        genes
        supporting_samples
        supporting_dbs
        finisher_manual_review
        pass_manual_review
        finisher_3730_review
        manual_genotype_normal
        manual_genotype_tumor
        manual_genotype_relapse
        somatic_status
        notes
        /);
    return @columns;
}

sub set_line_hash{
    my ($self, $data) = @_;
    my %hash;
    foreach my $col_name (@{$self->{separated_value_reader}->headers}){#changed from db_column above
        if ($col_name eq 'insert_sequence' && !exists $hash{'insert_sequence_allele1'}){
            my ($insert_sequence_allele1, $insert_sequence_allele2) ;
            ($insert_sequence_allele1, $insert_sequence_allele2) = split (/\//, $data->{$col_name}) if $data->{$col_name};
            $hash{'insert_sequence_allele1'} = $insert_sequence_allele1;
            $hash{'insert_sequence_allele2'} = $insert_sequence_allele2;
        }
        $hash{$col_name} = $data->{$col_name};        
    }
    $self->{hash} = \%hash;
}

sub line_hash{
    my ($self)=@_;
    unless ($self->{hash}){
        die "no current line data present!";
    }
    return $self->{hash};
}

sub next_line_data{
    my ($self) = @_;
    my $line_hash = $self->{separated_value_reader}->next;
    return unless $line_hash;
    foreach (keys %$line_hash){
        $_ =~ s/"'//g;          #strip quotation marks
        undef $_ if $_=~/^-$/;
    }

    $self->set_line_hash($line_hash);
    return $self->line_hash;
}

1;

=pod

=head1 Name

Genome::Utility::VariantReviewListReader

=head1 Synopsis

Processes a variant review list line by line, splits on the given separator, and builds a hash whose keys are the items from the header and values are the data from the processed row.     

=head1 Usage

  use Genome::Utility::VariantRevewListReader;

  my $reader = Genome::Utility::VariantReviewListReader->new(
     list => 'VariantReviewList.csv', # req; file
     separation_character => ',' , # opt; default is '|'
  );

  while (my $line_hash = $reader->next_line_data){
      print "Chromosome: " . $line_hash->{'chromosome'} . " and Begin Position: " . $line_hash->{'begin_position'};
  }

=head1 Methods

=head2 next_line_data

  my $line_data = $reader->next_line_data

=over

=item I<Synopsis>   Gets the next line in hash form from the list

=item I<Params>     none

=item I<Returns>    hashref with the following keys: chromosome, begin_position, end_position, variant_type, variant_length, delete_sequence, insert_sequence_allele1, insert_sequence_allele2, genes, supporting_samples, supporting_dbs, finisher_manual_review, pass_manual_review, finisher_3730_review, manual_genotype_normal, manual_genotype_tumor, manual_genotype_relapse, somatic_status, notes

=back

=head1 See Also

I<Genome::Utility::IO::SeparatedValueReader>

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

Adam Dukes  <adukes@watson.wustl.edu>,
Jim Weible  <jweible@watson.wustl.edu>

=cut
