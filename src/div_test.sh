p=$GRP/ukb; cd $p
s=dat/cal/001/{01..22}
o=$SCR/ukb/cal/bat/bvl

i=dat/phe/cnb.eid
j=dat/phe/smk.eid

u=dat/cal/002/{01..04}
t=dat/cal/002/[0-2][1-3]

sh src/div.sh -i $i -i $j -b 2048 -o $o -v --retain $t $u
sh src/div.sh -i $i -b 2048 -o $o -v --retain $t $u
sh src/div.sh -i $j -b 2048 -o $o -v --retain $t $u



p=$GRP/ukb; cd $p
s=dat/cal/002/{01..22}
o=dat/b2k/dbg
i='dat/phe/bvl/*.eid'
sh src/div.sh -i "$i" -b 2048 -o $o -v --retain $s
