#!/usr/bin/perl
# wangsen 20150212
#

my $dir=$ARGV[0];
my @file;
my %hash;
opendir DIR,"$dir" or die $!; 
for my $file(readdir DIR) { 
        next unless $file=~/^\w+/;
        open IN,"$dir/$file" or die $!; 
        push @file,$file; 
        while (<IN>) {
                chomp;
                my @array=split;
                #next unless ($array[1]=~/^\d+/);
                $hash{$array[0]}{$file}=$array[4]; 
        }
        close IN;
}
closedir DIR;
for my $i(keys %hash) {
        print "$i";
        for my $f(@file) {
                $c=($hash{$i}{$f})?$hash{$i}{$f}:0;
                print "\t$c";
        }
        print "\n";
}


