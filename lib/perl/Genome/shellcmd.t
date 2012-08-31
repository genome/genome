

use above Genome;
use Test::More tests => 6;


my $dir = '/tmp';
my $in = ["$dir/shellcmd1.test","$dir/shellcmd2.test"];
my $val = 'ian says its broked';


create($in->[0], $val);
create($in->[1], $val);
ok(cat(), 'two good input files');

unlink($in->[0]);
create($in->[1], $val);
ok(! cat(), 'one good input file, one missing file');

create($in->[0], '');
create($in->[1], '');
ok(cat(1), 'two empty input files (allow_zero_size_input_files = 1)');

create($in->[0], '');
create($in->[1], '');
ok(! cat(), 'two empty input files (allow_zero_size_input_files = 0)');

unlink($in->[0]);
unlink($in->[1]);
ok(! cat(), 'no input files (allow_zero_size_input_files = 0)');

unlink($in->[0]);
unlink($in->[1]);
ok(! cat(1), 'no input files (allow_zero_size_input_files = 1)');


exit;



sub create {
    my ($fn, $v) = @_;
    open(my $fh, ">$fn");
    print $fh $v;
    close($fh); 
}

sub cat {
    my ($allow_zero) = @_;

    eval {
        Genome::Sys->shellcmd(
            cmd => 'cat ' . join(' ',@$in) . '> /dev/null',
            input_files => $in,
            allow_zero_size_input_files => $allow_zero || 0
        );
    };
    return $@ eq '' ? 1 : undef;
}


