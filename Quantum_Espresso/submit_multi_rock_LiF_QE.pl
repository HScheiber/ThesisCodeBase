#!/usr/bin/perl

use POSIX qw/floor/;

######################
$QE_ver=0;
$hours0=3;
$mins0=00;

@l1i=(0,1,2,3,4,5,6,7,8,9);
@l1a=(7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0);

$l1l='ALAT';

$nMols_per_Task=-1; # -1 to fix the number of cores
$nCores=32;
$nTasks_per_Node=32;


$workDir=`pwd`;
chomp($workDir);

foreach $l1 (@l1i) {

 $l1s=sprintf("%.6f",$l1a[$l1]);
 $l1d="${l1l}_${l1s}";

 # create/change dir
 if (!(-d $l1d)) {
  mkdir($l1d);
 }
 chdir("$workDir/$l1d");

 if ($nMols_per_Task < 0) {
  $nTasks=$nCores;
  $hours_calc=$hours0;
  $mins_calc=$mins0;
  #$mins_calc=floor(((2**$l1)**$scaling)*($hours0*60+$mins0));
 } else {
  $nTasks=$l1a[$l1]/$nMols_per_Task;
  $hours_calc=$hours0;
  $mins_calc=$mins0;
 }

 $ids="$l1d";
 $sTasks=sprintf("%05d",$nTasks);

 $path = "$workDir";
 $path =~ s/\//\\\//g;

 system "sed -e \'s/##$l1l##/$l1s/g\' $workDir/LiF_rock.template > input.$l1s.inp";

 # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
 system "submit_QE.pl -1 $nTasks $nTasks_per_Node -1 -1 $hours_calc $mins_calc LiF-R-$l1s-$nCores input.$l1s.inp output.$l1s.out 1 $QE_ver";
  
 chdir("$workDir");

} #level1

