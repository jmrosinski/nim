#!/usr/bin/perl -p
s/\b(real)\(rt\)/$1/gi; s/([\d.])_rt\b/$1/gi;
