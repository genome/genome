#----------------------------------
# $Author: dlarson $ 
# $Date: 2008-08-28 15:15:41 -0500 (Thu, 28 Aug 2008) $
# $Revision: 37990 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram/Domain.pm $
#----------------------------------
package Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Domain; 
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
use warnings;
use SVG;

# If you know about and want to use inheritance:
#use Base::Class;
#our @ISA = qw( Base::Class );

sub new {
 	my ($class, %arg) = @_;
	my $self = {
							_backbone => $arg{backbone},
                            _style => $arg{style},
                            _start_amino_acid => $arg{start_aa},
                            _stop_amino_acid => $arg{stop_aa},
                            _text => $arg{text},
                            _id => $arg{id},
						 };
    #calculate the origin and scaling of the backbone                     
    my $backbone = $self->{_backbone};
    $self->{_x_start} = $backbone->x_coordinate_for_amino_acid_pos($self->{_start_amino_acid});
    $self->{_x_stop} = $backbone->x_coordinate_for_amino_acid_pos($self->{_stop_amino_acid});
    
    $self->{_feature_length} = ($self->{_x_stop} - $self->{_x_start});
    bless($self, ref($class) || $class);
    return $self;
}
#----------------------------------
sub length {
    my ($self,) = @_;
    return $self->{_feature_length};
}
#----------------------------------
sub draw {
    my ($self,) =  @_;
    my $backbone = $self->{_backbone};
    my $svg = $backbone->svg;
    $svg->rectangle(
        id => $self->{_id}, 
        x => $self->{_x_start},
        y => $backbone->dimensions()->{y},
        width => $self->{_feature_length},
        height => $backbone->dimensions->{height},
        style => $self->{_style}, 
    );
}
#----------------------------------
# TESTING CODE
#----------------------------------
sub normalize_bbox {
    my @box = @_;
    my $y_shift = 0 - $box[7];
    for(my $i = 1; $i < scalar(@box); $i += 2) {
        $box[$i] += $y_shift;
    }
    return @box;
}

=head1 AUTHOR

David Larson, E<lt>dlarson@watson.wustl.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 by David Larson

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8 or,
at your option, any later version of Perl 5 you may have available.


=cut
1;

