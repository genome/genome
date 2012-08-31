#!/usr/bin/env genome-perl

#############################################################

package TestStd;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use Finfo::Std;

    use Data::Dumper;

    my %i1 :name(i1:r) :isa('int');
    my %i2 :name(i2:o) :isa('int') :default(4);
    my %str_rw :name(str_rw:r) :isa(string) :default(default);
    my %str_ro :name(str_ro:o) :isa(string) :access(ro);
    my %_str :name(_str:p) :isa(string);

    sub START
    {
        my $self = shift;

        return $self->unlink_log_file;
    }

    sub log_file
    {
        return 'log.txt';
    }

    sub unlink_log_file
    {
        my $self = shift;

        unlink $self->log_file if -e $self->log_file;

        return 1;
    }
}

#############################################################

package TestStdInheritance;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use Data::Dumper;

    use base 'TestStd';

    my %i3 :name(i3:r) :isa('int <=> 100 110');
}

#############################################################

package TestObject;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use Data::Dumper;

    use base 'Finfo::Object';

    sub _attrs
    {
        return
        {
            'i1:r' => 
            {
                type => 'integer',
            },
            'i2:r' => 
            {
                type => 'integer', 
                default => 4, 
            },
            'str_rw:o' => 
            {
                type => 'defined',
            },
            'str_ro:o' => 
            {
                type => 'defined',
                access => 'ro',
            },
            '_str:p' => 
            {
                type => 'defined',
            },
        };
    }

    sub _init
    {
        my $self = shift;

        $self->_str('private!')
            or return;

        return 1;
    }
}

#############################################################

package TestObjectInheritance;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use Data::Dumper;

    use base 'TestObject';

    sub _attrs
    {
        my $self = shift;

        my $attrs = $self->SUPER::_attrs;

        $attrs->{'i3:r'} = { type => 'int_between', options => [qw/ 100 110/], };

        return $attrs;
    }
}

#############################################################

package TestOldObject;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use Data::Dumper;

    use base 'Finfo::Object';

    sub _reqs
    {
        return
        {
            i1 => [qw/ integer /],
            str => [qw/ defined /],
        };
    }

    sub _opts
    {
        return
        {
            i2 => [qw/ integer 4 /],
        };
    }
}

#############################################################

package Person;
{
    use Finfo::Std;

    my %name :name(name:o) :isa(string) :default(Ed);
    my %height :name(height:r) :isa('real pos');
    my %weight :name(weight:r) :isa('int pos');
    my %age :name(age:r) :isa('int pos');
}

#############################################################

package Peep;
{
    use Finfo::Std;

    my %name :name(name:o) :isa(string) :default(Ed);
}

#############################################################

package TestReader;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use Data::Dumper;

    use base 'Finfo::Reader';

    my $count = 0;

    sub _return_class
    {
        return 'Person';
    }

    sub _next
    {
        my $self = shift;

        return if ++$count == 4;

        return 
        {
            height => 6.25,
            weight => 215,
            age => 31,
        };
    }
}

#############################################################

package TestWriter;

{
    use strict;
    use warnings;
    no warnings 'reserved';

    use base 'Finfo::Writer';

    sub _write_one
    {
        my $self = shift;

        return 1 if $self->io;
    }
}

#############################################################

package TestSingleton;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use base 'Finfo::Singleton';

    use Data::Dumper;

    my %color :name(color:r) :isa('in_list red green purple');

    sub test_enforce_instance
    {
        my $self = shift;

        $self->_enforce_instance
            or return;

        return 1;
    }
}

#############################################################

package FinfoTest;

use strict;
use warnings;
no warnings 'reserved';

use base qw(Test::Class);

use Data::Dumper;
use Finfo::Iterator;
use Finfo::Msg;
use Finfo::Validate;
use IO::String;
use Test::More;

require Finfo::Logging;
require Storable;

use constant TYPE => 0;
use constant GOOD => 1;
use constant BAD => 2;
use constant OPTIONS => 3;
use constant REF_TYPE => 4;

sub test_logging_msgs
{
    my ($self, $obj) = @_;

    foreach my $level ( Finfo::Msg->levels )
    {
        next if $level eq 'fatal';
        my $method = $level . '_msg';
        $obj->$method("Test level ($level)");
    }

    return 1;
}

sub test01_std : Tests
{
    my $self = shift;

    my $class = 'TestStd';
    print "\nTesting $class\n";

    ok
    (
        Finfo::Logging::add_appender
        (
            class => $class,
            type => 'file',
            params => { filename => $class->log_file, mode => 'clobber' },
            pattern => '%m (%d)%n',
            level => 'DEBUG',
        ), "Added file appender to $class"
    );

    my %p = 
    (
        i1 => 55,
        str_rw => 'simba',
        str_ro => 'nala',
    );

    my $obj = $class->new(%p);
    ok($obj, "create $class");

    # save to use later for validating
    $self->{val_obj} = $obj;

    foreach my $attr ( keys %p )
    {
        my $val = $obj->$attr;
        ok($val eq $p{$attr}, "Check init attr $attr, $val eq $p{$attr}");
        die unless $val;
    }

    ok($obj->i2 eq 4, "Check default set for i2 ".$obj->i2);

    # undef attr
    ok( $obj->undef_attribute('str_rw'), 'Undef attr "simba"');
    $obj->str_rw('simba');
    
    my $set_priv;
    eval
    {
        $set_priv = $obj->_str('set private attr');
    };

    is($set_priv, undef, "Unable to set private attr");
    
    my $set_ro;
    eval
    {
        $set_ro = $obj->str_ro('set read only attr');
    };

    is($set_ro, undef, "Unable to set read only attr");

    ok
    (
        Finfo::Logging::set_screen_logger_level('debug'),
        'Set screen logger to debug, will print all debug msgs',
    );

    $self->test_logging_msgs($obj);
    
    # Test another object writing to log file
    my $obj2 = $class->new(%p);
    ok($obj2, "Created object 2");

    $obj2->info_msg("str_rw is now " . $obj2->str_rw);

    #my $fh = IO::File->new('<' . $obj->log_file);
    #die $obj->log_file ." $!\n" unless $fh;
    #my $log_string = join('', $fh->getlines);
    #$fh->close;

    #ok($log_string, "Logged msgs from $class:\n$log_string");

    # Test dump
    my $dump1 = $obj->DUMP;
    ok($dump1, "Dumped obj1 from $class\:\n$dump1");
    
    # Test storable
    my $serialized = Storable::freeze($obj);
    ok($serialized, "Froze $class");

    my $clone = Storable::thaw($serialized);
    ok($clone, "Thawed $class");

    my $dump2 = $clone->DUMP;
    ok($dump2, "Dumped clone of obj1 from $class\:\n$dump2");

    is_deeply($dump1, $dump2, "Compared obj1 dump to obj1 clone dump");
        
    #$obj->unlink_log_file;

    return 1;
}

sub test02_std_inheritance : Tests
{
    my $self = shift;

    my $class = 'TestStdInheritance';
    print "\nTesting $class\n";

    my %p = 
    (
        i1 => 55,
        i3 => 105,
        str_rw => 'simba',
        str_ro => 'nala',
    );

    my $obj = $class->new(%p);
    ok($obj, "create $class");

    die unless $obj;
    
    foreach my $attr ( keys %p )
    {
        my $val = $obj->$attr;
        ok($val eq $p{$attr}, "check attr $attr, $val eq $p{$attr}");
    }

    ok($obj->i2 eq 4, "Check default set for i2");

    # undef attr
    ok( $obj->undef_attribute('i3'), 'Undef attr "simba"');
    $obj->str_rw(400);
    
    my $get_priv;
    eval
    {
        $get_priv = $obj->_str('get private attr');
    };

    is($get_priv, undef, "Unable to get private attr");
 
    my $set_priv;
    eval
    {
        $set_priv = $obj->_str('set private attr');
    };

    is($set_priv, undef, "Unable to set private attr");
    
    my $set_ro;
    eval
    {
        $set_ro = $obj->str_ro('set read only attr');
    };

    is($set_ro, undef, "Unable to set read only attr: $@");

    $self->test_logging_msgs($obj);

    return 1;
}

sub test03_object : Tests
{
    my $self = shift;

    my $class = 'TestObject';
    print "\nTesting $class\n";

    my %p = 
    (
        i1 => 55,
        str_rw => 'simba',
        str_ro => 'nala',
    );

    my $obj = $class->new(%p);
    ok($obj, "create $class");

    foreach my $attr ( keys %p )
    {
        my $val = $obj->$attr;
        ok($val eq $p{$attr}, "check attr $attr, $val eq $p{$attr}");
    }

    ok($obj->i2 eq 4, "Check default set for i2");

    ok
    ( 
        Finfo::Logging::set_msg_detail_level($class, 'more'), 
        "Set $class msg detail level to more"
    );
    
    my $get_priv;
    eval
    {
        $get_priv = $obj->_str('get private attr');
    };

    is($get_priv, undef, "Unable to get private attr");
 
    my $set_priv;
    eval
    {
        $set_priv = $obj->_str('set private attr');
    };

    is($set_priv, undef, "Unable to set private attr");
    
    my $set_ro;
    eval
    {
        $set_ro = $obj->str_ro('set read only attr');
    };

    is($set_ro, undef, "Unable to set read only attr: $@");

    $self->test_logging_msgs($obj);

    return 1;
}

sub test04_object_inheritance : Tests
{
    my $self = shift;

    my $class = 'TestObjectInheritance';
    print "\nTesting $class\n";

    my %p = 
    (
        i1 => 55,
        i3 => 105,
        str_rw => 'simba',
        str_ro => 'nala',
    );

    my $obj = $class->new(%p);
    ok($obj, "create $class");

    die unless $obj;
    
    foreach my $attr ( keys %p )
    {
        my $val = $obj->$attr;
        ok($val eq $p{$attr}, "check attr $attr, $val eq $p{$attr}");
    }

    ok($obj->i2 eq 4, "Check default set for i2");

    my $get_priv;
    eval
    {
        $get_priv = $obj->_str('get private attr');
    };

    is($get_priv, undef, "Unable to get private attr");
 
    my $set_priv;
    eval
    {
        $set_priv = $obj->_str('set private attr');
    };

    is($set_priv, undef, "Unable to set private attr");
    
    my $set_ro;
    eval
    {
        $set_ro = $obj->str_ro('set read only attr');
    };

    is($set_ro, undef, "Unable to set read only attr: $@");

    $self->test_logging_msgs($obj);

    return 1;
}

sub test04_oldobject : Tests
{
    my $self = shift;

    my $class = 'TestOldObject';
    print "\nTesting $class\n";

    my %p = 
    (
        i1 => 55,
        str => 'simba',
    );

    my $obj = $class->new(%p);
    ok($obj, "create $class");

    foreach my $attr ( keys %p )
    {
        my $val = $obj->$attr;
        ok($val eq $p{$attr}, "check attr $attr, $val eq $p{$attr}");
    }

    ok($obj->i2 eq 4, "Check default set for i2");

    ok
    (
        Finfo::Logging::set_screen_logger_level('info'),
        'Set screen logger to info, will NOT print debug msgs',
    );

    # works - prints a lot...
    #ok
    #( 
    #    Finfo::Logging::set_msg_detail_level($class, 'stack'), 
    #    "Set $class msg detail level to stack"
    #);

    $self->test_logging_msgs($obj);

    return 1;
}

sub test05_validations : Tests
{
    my $self = shift;

    #return 1;
    my @validation_types = 
    (
        [ 'defined', 0, undef],
        [qw/ positive_integer 1 0 /], 
        [qw/ non_negative_integer 0 -1/],
        [qw/ integer 1 t /], 
        [qw/ integer_gt 8 0 /, [ 5 ] ], 
        [qw/ integer_gte 8 0 /, [ 5 ] ], 
        [qw/ integer_lt 0 11 /, [ 5 ] ], 
        [qw/ integer_lte 0 8 /, [ 5 ] ], 
        [qw/ integer_between 8 0 /, [ 5, 10 ] ], 
        [qw/ pos_int 1 0 /], 
        [qw/ non_neg_int 0 -1/],
        [qw/ int 1 1.1 /], 
        [qw/ int_gt 8 0 /, [ 5 ] ], 
        [qw/ int_gte 8 0 /, [ 5 ] ], 
        [qw/ int_lt 0 11 /, [ 5 ] ], 
        [qw/ int_lte 0 8 /, [ 5 ] ], 
        [qw/ int_between 8 0/, [ 5, 10 ] ], 
        [qw/ real 1.3 t /], 
        [qw/ positive_real 13 -22 /],
        [qw/ non_negative_real 0.25 -44 /],
        [qw/ real_gt 8.3 0.4 /, [ 5 ] ], 
        [qw/ real_gte 5.1 0.6 /, [ 5 ] ], 
        [qw/ real_lt 3.5 6.7 /, [ 5 ] ], 
        [qw/ real_lte 3.7 9.5 /, [ 5 ] ], 
        [qw/ real_between 0 -1 /, [ 0, 10 ] ], 
        [qw/ float 1.3 t /], 
        [qw/ positive_float 13 -22 /],
        [qw/ non_negative_float 0.25 -44 /],
        [qw/ float_gt 8.3 0.4 /, [ 5 ] ], 
        [qw/ float_gte 5.1 0.6 /, [ 5 ] ], 
        [qw/ float_lt 3.5 6.7 /, [ 5 ] ], 
        [qw/ float_lte 3.7 9.5 /, [ 5 ] ], 
        [qw/ float_between 0 -1/, [ 0, 10 ] ], 
        [ 'aryref', [], 1 ],
        [ 'aryref', [qw/ 1 2 3 /], [qw/ 1 4.3 4 /], undef, 'int' ],
        [ 'non_empty_aryref',['test'], [] ],
        [ 'non_empty_aryref', [qw/ 1 2 3 /], [qw/ 1 4 3 ff /], undef, 'int' ],
        [ 'hashref', {}, 1 ],
        [ 'non_empty_hashref', { hash => 'ref' }, {} ],
        [ 'non_empty_hashref', { hash => 6 }, { hash => 6, re => 66 }, [qw/ 5 7 /], 'integer_between' ],
        [qw| file /gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt . |],
        [qw| input_file  /gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt nowaythisexists |],
        [qw| output_file jjj /gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt |],
        [qw/ path . hoopla_dir /],
        [qw/ input_path . hoopla_dir /],
        [qw/ output_path . hoopla_dir /],
        [qw/ io_path . hoopla_dir /],
        [qw/ dir . hoopla_dir /],
        [qw/ input_dir . hoopla_dir /],
        [qw/ output_dir . hoopla_dir /],
        [qw/ io_dir . hoopla_dir /],
        [qw/ in_list yes NO /, [qw/ y yes YES/] ],
        [qw/ regex good not_good /, [qw/ ^go ^super /] ],
        [ 'regex_aryref', [qw/ good goodest /], [qw/ bad baddest not_bad /], [qw/ ^good ^bad /] ],
        [ 'code',  sub{ print }, 'not code' ],
        [qw/ is_validation_type y_or_n not_a_type /],
        [qw/ y_or_n no hya /],
        [ 'not_blank', 'not', '' ],
        [ 'string', 'a string', { not_a => 'string' } ],
        [ 'non_empty_aryref', ['a really long string'], ['short'], [qw/ int_gt 7 /], 'string' ],
        #[qw/ project_name M_BB0392D19 not_a_project /],
    );

    my $class = 'Finfo::Validate';
    print "\nTesting $class\n";
    my $val_obj = $self->{val_obj};

    foreach my $validation_info ( @validation_types )
    {
        my $type = $validation_info->[TYPE];
        my $good_val = $validation_info->[GOOD];
        my $bad_val = $validation_info->[BAD];
        my $opts = $validation_info->[OPTIONS];
        my $ref_type = $validation_info->[REF_TYPE];

        $val_obj->info_msg("Testing validation - $type");

        ok
        (
            $class->validate
            (
                attr => $type,
                value => $good_val,
                type => $type,
                options => $opts,
                ref_type => $ref_type,
                err_cb => $val_obj,
            ),
            sprintf
            (
                '%s w/ good val (%s)',
                $type,
                ( ref($good_val) ) ? Dumper($good_val) : $good_val,
            )
        );

        my $bad_result;
        eval
        {
            $bad_result = $class->validate
            (
                attr => $type,
                value => $bad_val,
                type => $type,
                options => $opts,
                ref_type => $ref_type,
                err_cb => $val_obj,
            );
        };

        is
        (
            $bad_result,
            undef,
            sprintf
            (
                '%s w/ bad value (%s)',
                $type,
                ( defined $bad_val ) 
                ? ( ref($bad_val) ) ? Dumper($bad_val) : $bad_val 
                : 'not defined'
            )
        );
    }
}

sub test05_validations2 : Tests
{
    my $self = shift;

    #return 1;
    my @validation_infos = 
    (
        {
            isa => 'defined',
            g => 0,
            b => undef,
        },
        {
            isa => 'int',
            g => 1,
            b => 0.55,
        },
        {
            isa => 'int pos',
            g => 1,
            b => 0,
        },
        {
            isa => 'int neg',
            g => -1,
            b => 1,
        },
        {
            isa => 'int non_neg',
            g => 0,
            b => -1,
        },
        {
            isa => 'int gt 6',
            g => 7,
            b => 0,
        },
        {
            isa => 'int gte 6',
            g => 6,
            b => -5,
        },
        {
            isa => 'int lt 6',
            g => -5,
            b => 7,
        },
        {
            isa => 'int lte 6',
            g => -6,
            b => 85,
        },
        {
            isa => 'int between 6 8',
            g => 7,
            b => -5,
        },
        {
            isa => 'real',
            g => 1.111,
            b => 'g',
        },
        {
            isa => 'real pos',
            g => 1.6,
            b => 0,
        },
        {
            isa => 'real neg',
            g => -1.0202,
            b => 1.5,
        },
        {
            isa => 'real non_neg',
            g => 0.555,
            b => -1.77,
        },
        {
            isa => 'real gt 6',
            g => 7.64323,
            b => 0,
        },
        {
            isa => 'real gte 6',
            g => 6.36,
            b => -5.88,
        },
        {
            isa => 'real lt 6',
            g => -5.55,
            b => 7,
        },
        {
            isa => 'real lte 6',
            g => -6,
            b => 55.4,
        },
        {
            isa => 'real between 6 8',
            g => 7.1,
            b => -5,
        },
        {
            ds => 'aryref',
            empty_ok => 1,
            g => [],
            b => 1,
        },
        {
            ds => 'aryref',
            g => [qw/ 1 2 3 /],
            b => [],
        },
        {
            ds => 'aryref',
            isa => 'int',
            g => [qw/ 1 2 3 /],
            b => [qw/ 1 2 f /],
        },
        {
            ds => 'hashref', 
            empty_ok => 1,
            g => {},
            b => 1,
        },
        {
            ds => 'hashref', 
            g => { key => 1 },
            b => {},
        },
        {
            ds => 'hashref', 
            isa => 'int between 3 4',
            g => { key1 => 3, key2 => 4 },
            b => { key3 => 3, key4 => 'y' },
        },
        {
            isa => 'file',
            g => '/gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt',
            b => '.',
        },
        {
            isa => 'file_r',
            g => '/gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt',
            b => 'nowaythisexists',
        },
        {
            isa => 'file_w',
            g => 'nowaythisexists',
            b => '/gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt',
        },
        {
            isa => 'file_rw',
            g => '/gsc/var/cache/testsuite/data/Genome-Utility-Filesystem/existing_file.txt',
            b => '/lost+found/nowaythisexists',
        },
        {
            isa => 'dir',
            g => '.',
            b => 'hoopla_dir',
        },
        {
            isa => 'dir_r',
            g => '.',
            b => 'hoopla_dir',
        },
        {
            isa => 'dir_w',
            g => '.',
            b => '/',
        },
        {
            isa => 'dir_rw',
            g => '.',
            b => '/',
        },
        {
            isa => 'in_list foo bar',
            g => 'foo',
            b => 'yes',
        },
        {
            ds => 'aryref',
            isa => 'in_list foo bar',
            g => [qw/ foo bar /],
            b => [qw/ foo yes YES/],
        },
        {
            isa => 'regex ^go',
            g => 'good',
            b => 'bad',
        },
        {
            ds => 'aryref',
            isa => 'regex ^good',
            g => [qw/ goody /],
            b => [qw/ good bad /],
        },
        {
            isa => 'code',
            g => sub{},
            b => 'not code',
        },
        {
            isa => 'y_or_n',
            g => 'no',
            b => 'h',
        },
        {
            isa => 'not_blank',
            g => 'not',
            b => '',
        },
        {
            isa => 'string',
            g => 'a string', 
            b => { not_a => 'string' },
        },
        {
            ds => 'aryref',
            isa => 'string > 7',
            g => ['a really long string'],
            b => [ 'short' ],
        },
        {
            isa => 'object',
            g => $self->{val_obj},
            b => 'not an object',
        },
        {
            isa => 'object TestStd',
            g => $self->{val_obj},
            b => IO::String->new(),
        },
    );

    my $class = 'Finfo::Validate';
    print "\nTesting $class\n";
    my $val_obj = $self->{val_obj};

    foreach my $val_info ( @validation_infos )
    {
        my $isa = Dumper($val_info->{isa});
        if ( $isa =~ /undef/ )
        {
            $isa = 'no isa';
        }
        else
        {
            $isa =~ s/\n//g;
            $isa =~ s/\$VAR1 = //;
            $isa =~ s/;//g;
            $isa =~ s/\s+/ /g;
        }

        $isa .= " ds " . $val_info->{ds} if $val_info->{ds};

        $val_info->{attr} = $isa;
        $val_info->{msg} = 'fatal';
        $val_info->{obj} = $self->{val_obj},
        
        my $gv = delete $val_info->{g};
        my $bv = delete $val_info->{b};

        my %gvh = %$val_info;
        $gvh{value} = $gv;

        my %bvh = %$val_info;
        $bvh{value} = $bv;  

        $val_obj->info_msg("Testing validation - $isa");

        ok
        (
            $class->validate(%gvh),
            #sprintf
            #(
                'w/ good val',
                #'w/ good val (%s)',
                #( ref($gv) ) ? Dumper($gv) : $gv,
                #)
        );

        my $bad_result;
        eval
        {
            $bad_result = $class->validate(%bvh);
        };

        is
        (
            $bad_result,
            undef,
            #sprintf
            #(
                'w/ bad value',
                #'w/ bad value (%s)',
                #   ( defined $bv ) 
                #? ( ref($bv) ) ? Dumper($bv) : $bv 
                #: 'undef'
                #)
        );
    }

    # test is_isa and ds
    my $ds_and_isa = 
    {
        is_isa =>
        {
            good => 'y_or_n',
            bad => 'not_an_isa',
        },
        is_ds =>
        {
            good => 'hashref', 
            bad => 'big_hashref',
        },
    };

    foreach my $method ( keys %$ds_and_isa )
    {
        $val_obj->info_msg("Testing $method");
        my $good = $ds_and_isa->{$method}->{good};
        my $bad = $ds_and_isa->{$method}->{bad};
        ok
        (
            $class->$method
            (
                attr => "testing good $method",
                value => $good,
                obj => $val_obj,
                msg => 'fatal',
            ),
            "w/ good $method val ($good)",
        );

        my $bad_result;
        eval
        {
            $bad_result = $class->$method
            (
                attr => "testing bad $method",
                value => $bad,
                obj => $val_obj,
                msg => 'fatal',
            );
        };

        is
        (
            $bad_result,
            undef,
            "w/ bad $method value ($bad)",
        );
    }

    return 1;
}

sub test06_reader : Tests
{
    my $self = shift;

    my $class = 'TestReader';
    print "\nTesting $class\n";
    my $obj = $class->new(io => IO::String->new());
    ok($obj, "Created $class");
    die unless $obj;

    my @refs = $obj->all;
    ok(@refs == 3, sprintf("Got %d refs", scalar(@refs)));
    is(ref($refs[0]), 'HASH', "Got hashrefs");
    
    $obj = undef;
    $obj = $class->new
    (
        io => IO::String->new(),
        return_as_objs => 1,
    );
    ok($obj, "Created reader");
    
    my $ed = $obj->next;
    is(ref($ed), 'Person', "Got Ed");

    return 1;
}

sub test07_writer : Tests
{
    my $self = shift;

    my $class = 'TestWriter';
    print "\nTesting $class\n";
    my $obj = $class->new(io => IO::String->new());
    ok($obj, "Created $class");
    die unless $obj;

    ok($obj->write_one(1), "Wrote one");
    is($obj->write_many([1, 2, 3]), 3, "Wrote many - 3");

    return 1;
}

sub test08_singleton : Tests
{
    my $self = shift;

    my $class = 'TestSingleton';
    print "\nTesting $class\n";
    my $obj1 = $class->instance(color => 'red');
    ok($obj1, "Created $class");
    die unless $obj1;

    ok($obj1->test_enforce_instance, "Enforce instance");
    
    $obj1->info_msg("color is " . $obj1->color);
    $obj1->color('green');
    ok($obj1->color eq 'green', "set color to green");
    $obj1->info_msg("color is " . $obj1->color);
    
    my $instance_with_params;
    eval
    {
        $instance_with_params = $class->instance(color => "purple");
    };
    is($instance_with_params, undef, "Unable to get singleton once created with instance (params)");
 
    my $obj2 = $class->instance;
    ok($obj2, "Got singleton instance");

    ok($obj1->color('purple'), "Set color to purple on obj1");
    ok($obj2->color eq 'purple', "obj2 color is purple, set on obj1");

    my $bad_color;
    eval
    {
        $bad_color = $obj1->color("yellow");
    };
    is($bad_color, undef, "Unable to set color to 'yellow'");

    return 1;
}

sub test09_clo : Tests
{
    my $self = shift;
    # TODO
    return 1;
}

sub test10_iterator : Tests
{
    my $self = shift;
    # count reset first last next all find search

    print "\nTesting Finfo::Iterator\n";
    
    my $names = [qw/ Sally Jon Kyung Ed Neha Edward Adam /];
    my $iterator = Finfo::Iterator->new
    (
        ids => $names,
        cb => sub{ return Peep->new(name => $_[0]); },
        id_method => 'name',
    );
    is($iterator->count, scalar @$names, sprintf('Iterator has %d values', scalar @$names));
    is($iterator->first->name, 'Sally', 'Sally is first');
    is($iterator->last->name, 'Adam', 'Adam is last');

    my $i = 0;
    while ( my $person = $iterator->next )
    {
        ok($person->name eq $names->[$i], "Got " . $names->[$i]);
        $i++;
    }
    ok($iterator->reset, "Reset iterator");
    is_deeply([ map { $_->name } $iterator->all ], $names, "Got all");
    my $ed = $iterator->find({ name => 'Ed' });
    ok($ed && $ed->name eq 'Ed', 'Found Ed');
    my $eds = $iterator->search({ name => 'ed' }, { name => 'like' });
    is($eds->count, 2, "Searched for eds");
    is_deeply([ map { $_->name } $eds->all ], [qw/ Ed Edward /], "Got the Eds");
    
    return 1;
}

#################################################
#################################################

package main;

use strict;
use warnings;
no warnings 'reserved';

Test::Class->runtests('FinfoTest');

exit 0;


