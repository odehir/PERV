#!/usr/bin/perl -w 

use strict;

my $file = shift;


my %list;
$list{"5ltr"}= "0-564";
#$list{"5ltr"}= "1-564";
$list{"3ltr"}= "7518-8200";
#$list{"3ltr"}= "7518-7976";

$list{gag}= "585-2159";
#$list{pol}= "2307-5744";
$list{pol}= "2307-5747";
$list{env}= "5620-7536";

my $list;
foreach my $i (keys %list){
  my @tmp = split /\-/, $list{$i};
  $list->{$i}->{s} = $tmp[0];
  $list->{$i}->{l} = $tmp[1];
}

my $seq;
my $n = 0;

if($file =~ /.gz$/){
  open FILE, "gunzip -cd $file";
}else{
  open FILE, "$file";
}
while(<FILE>){
  s/\A\s+//;
  s/\s+\Z//;
  my $line = $_;

  if(/^\>([\W\w]+)$/){
    my $name = $1;
    $n++;
    $name =~ s/\,//g;
    $seq->{$n}->{n} = $name;
    $seq->{$n}->{s} = "";
  }else{
    $seq->{$n}->{s} = $seq->{$n}->{s} . uc($line); 
  }
}
close FILE;

printf "#sequence ID;REF,";
printf "length(5'-LTR),length(3'-LTR),";

printf "#C+m(5'-LTR),#C+m(3'-LTR),";
printf "#C+h(5'-LTR),#C+h(3'-LTR),";
printf "NA sequence(5'-LTR),";
printf "NA sequence(3'-LTR),";
print "\n";

my @ref = split //, $seq->{1}->{s};
for(my $i = 2; $i <= $n; $i++){
  my %len;
  my %numM;
  my %numH;
  my %ss;
  my @tmp = split //, $seq->{$i}->{s};
  my $j = 0;
  for(my $i2 = 0; $i2 <= $#ref; $i2++){
    if($ref[$i2] =~ /[A-Z]/){
      $j++;
    }

    foreach my $id ("5ltr", "3ltr"){
      if(! exists $len{$id}){ $len{$id} = 0; }
      if(! exists $ss{$id}){ $ss{$id} = ""; }
      if(! exists $numM{$id}){ $numM{$id} = 0; }
      if(! exists $numH{$id}){ $numH{$id} = 0; }

      if($j >= $list->{$id}->{s} && $j <= $list->{$id}->{l}){
        if($tmp[$i2] =~ /[A-Z]/){
          $len{$id} = $len{$id} + 1;
          $ss{$id} = $ss{$id} . $tmp[$i2];

          if($tmp[$i2] =~ /M/){ $numM{$id} = $numM{$id} + 1; }
          if($tmp[$i2] =~ /H/){ $numH{$id} = $numH{$id} + 1; }
        }
      }
    }
  }

  printf "%s,", $seq->{$i}->{n};

  printf "%d,%d,", $len{"5ltr"}, $len{"3ltr"};
  printf "%d,%d,", $numM{"5ltr"}, $numM{"3ltr"};
  printf "%d,%d,", $numH{"5ltr"}, $numH{"3ltr"};

  printf "%s,", $ss{"5ltr"};
  printf "%s,", $ss{"3ltr"};
  print "\n";
}



