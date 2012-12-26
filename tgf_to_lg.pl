#!/usr/bin/perl -w
use strict;
use warnings;

print "t # 0\n";

my $prefix = 'v ';

my $line;
while (defined($line = <STDIN>))
{
    my @values = split(' ', $line);
    if (@values)
    {
	if ($values[0] =~ "^#.*")
	{
	    $prefix = 'e ';
	}
	else
	{
	    print $prefix . $values[0] . " " . $values[1]. "\n";
	}
    }
}
