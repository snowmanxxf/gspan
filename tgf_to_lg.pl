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
	    if ($prefix =~ 'e')
	    {
		print $prefix . ($values[0]-1) . " " . ($values[1]-1) . "\n";
	    }
	    else
	    {
		print $prefix . ($values[0]-1) . " " . $values[1] . "\n";
	    }
	}
    }
}
