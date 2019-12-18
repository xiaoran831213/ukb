#!/bin/sh
# saner programming env: these switches turn some bugs into errors
set -o errexit
set -o pipefail
# set -o noclobber # when set, bash does not overite with >, >> or <>
set -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi

# basic input
declare -a _i                   # individual IDs?
declare -a _g                   # genotype (multiple)
_o=                             # output
_b=1024                         # batch size
_w=                             # working directory
_v=                             # verbose?
_e=                             # erase existing files?
retain=n                        # retain temp files?

SOPT=i:g::b:o:evc
LOPT=iid:,genotype::,batchsize:,out:,workdir,erase,verbose,retain

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$SOPT --longoptions=$LOPT --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi

# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

# now enjoy the options in order and nicely split until we see --
echo "Parsing command line:"
while true; do
    # echo "$1 $2"
    case "$1" in
        -g|--gno)          s=2; _g+=("$2") ;;
        -i|--iid)          s=2; _i+=("$2") ;;
        -b|--batchsize)    s=2; _b="$2" ;;
        -e|--erase)        s=1; _e=y ;;
        -v|--verbose)      s=1; _v=y ;;
        -o|--output)       s=2; _o="$2" ;;
        --retain)          s=1; retain=y ;;
        --)                shift; break ;;
        * ) echo "error parsing arguments"; exit 3 ;;
    esac
    echo -n "  $1"
    if [ $s -eq 2 ]; then echo -n " $2"; fi
    shift $s
    echo
done
# display positional arguments: the genotype words
echo "  -- $@"
echo

# process positional arguments: the genotype words
if [[ $# -lt 1 ]]; then
    echo "$0: require at least one genotype word."
    exit 4
fi
_g=("$@")                       # array of genotype words
_m=${#_g[@]}                    # number of words
echo "$# genotype words specified."
[ $_v ] && for _ in ${_g[@]}; do echo '  '$_; done
echo
echo
# -------------------- done arguement parsing -------------------- #


# --------------------    helper functions    -------------------- #
# sneak peeking a text file
peek() {
    [ $# -gt 0 ] || return 255  # file unspecified
    [ ! -e $#  ] || return 254  # file not found
    
    local n=$(cat $1 | wc -l)   # line count

    local l=5                   # line limit
    if [ $# -gt 1 ]; then l=$2; fi

    # peek now
    if [ $n -gt $[l*2] ]; then
        head -n$l "$1"
        echo ...
        tail -n$l "$1"
    else
        cat $1
    fi
}

# manage output and working directories
echo MSG: create output and working  directories
_o=${_o:-"$_o".}                # output
_w=${_w:-"$_o"/wrk}             # working
mkdir -p $_o $_w                # create directories
echo out="$_o"
echo wrk="$_w"
echo


# ------------------------- working now -------------------------- #

# collect genotype files
for i in $(seq 0 $[$_m-1]); do
    eval "ls ${_g[$i]}*.bed"
done | sed 's/[.]bed$//' | sort -u > "$_w/gls"
_m=($(wc -l "$_w/gls"))
echo "$_w/gls": $_m genotype files.
peek "$_w/gls"
echo

# list individuals appeared in all genotypes
echo "list IID shared by all genotype:"
if [ -s "$_w/gid" ]; then       # exist/erase?
    s=$([ $_e ] && echo ERASE || echo EXIST)
else
    s=CREATE                     # or new?
fi
if [ $s = ERASE -o $s = CREATE ]; then
    while IFS= read -r g; do
        awk <"$g.fam" '{print $2}'
    done < "$_w/gls" \
        | sort | uniq -c \
        | awk -v m=$_m '$1==m {print $2}' > "$_w/gid"
    r=${PIPESTATUS[0]}
else
    r=PASS
fi
echo "$s    $r    $_w/gid"      # report, exit on error
[ $r = "0" ] || [ $r = PASS ] || exit $r
echo

_N=($(wc -l "$_w/gid"))         # report
echo \""$_w/gid"\": $_N IID saved.
[ $_v ] && peek "$_w/gid"
echo


# gether individuals retain some individuals
echo MSG: filter IID if necessary
if [ -s "$_w/iid" ]; then       # exist/erase?
    s=$([ $_e ] && echo ERASE || echo EXIST)
else
    s=CREATE                    # or new?
fi
if [ $s = ERASE -o $s = CREATE ]; then
    if [ "$_i" ]; then          # -i/-iid: some individuals
        echo ${#_i[@]} IID words specified by --iid/-i,
        for i in "${_i[@]}"; do
            [ $_v ] && echo '  '$i
            fs+=($(eval "ls $i"))
        done
        echo which expands to ${#fs[@]} IID files.
        [ $_v ] && for f in "${fs[@]}"; do echo "  $f"; done
        echo

        echo MSG: sorting and merging IID files
        sort -u "${fs[@]}" > "$_w/kid"
        r=$?; [ ! $r -eq 0 ] && echo ERR: $r && exit $r

        echo MSG: joining sorted IID with genotyped IID
        join "$_w/kid" "$_w/gid" > "$_w/iid"
        r=$?; [ ! $r -eq 0 ] && echo ERR: $r && exit $r
    else                        # allow all individuals
        cp "$_w/gid" "$_w/iid"
        r=COPY
    fi
else
    r=PASS
fi
echo "$s    $r    $_w/iid"
_n=($(wc -l $_w/iid))
echo "$_w/iid:" $_n IID retained.
[ $_v ] && peek "$_w/iid"
echo


# shuffle and duplicate it as FID
echo "MSG: shuffle IID, and duplicate as FID:"
sort -R "$_w/iid" | awk '{print $1,$1}'> "$_w/cid"
echo "$_w/cid": $_n shuffled FID and IID saved.
[ $_v ] && peek "$_w/cid"
echo


# batch division, also get the number of batches
echo MSG: divide $_n IID into batches of $_b each
_q=$[ $_n / $_b ]
if [ -d "$_w/div" ] && \
       [ $(find "$_w/div" -name "*.cid" | wc -l) -eq $_q ] && \
       [ $(cat "$_w/div/"*.cid | wc -l) -eq $_n ]; then
    s=$([ $_e ] && echo ERASE || echo EXIST)
else
    mkdir -p "$_w/div"
    s=CREATE                    # or new?
fi
if [ $s = CREATE -o $s = ERASE ]; then
    split -n l/$_q -d -a3 "$_w/cid" "$_w/div/" --additional-suffix ".cid"
    r=$?
else
    r=PASS
fi
echo "$s    $r    $_w/div"
[ $r = "0" -o $r = PASS ] || exit $r
echo "$_w/div/": $_q batches
wc -l "$_w/div"/*.cid | head -n-1 | awk '{print $1,$2}' > "$_w/bsz"
[ $_v ] && peek "$_w/bsz"
echo


echo "MSG: divide genotype accordingly"
while IFS= read -r f            # genotype files
do
    g=${f//[ \/]/_}             # genotype id
    t="$_w/div/gno.$g"          # genotype extract
    echo GNO: $f $t
    if [ -s "$t.bed" -a -s "$t.bim" -a -s "$t.fam" ]; then
        s=$([ $_e ] && echo ERASE || echo EXIST)
    else
        s=CREAT; fi
    
    if [ $s = ERASE -o $s = CREAT ]; then
        plink --bfile "$f" --keep "$_w/cid" --keep-allele-order \
          --make-bed --out "$t" &> /dev/null
        r=$?
    else
        r=PASS; fi
    echo "$s    $r    $t"   # report, exit on error
    [ $r = "0" ] || [ $r = PASS ] || exit $r

    # divide the extracted genotype
    for b in "$_w/div/"*.cid    # batch numbers
    do
        o="${b%.cid}.${g}"
        if [ -s "$o.bed" -a -s "$o.bim" -a -s "$o.fam" ]
        then
            s=$([ $_e ] && echo ERASE || echo EXIST)
        else
            s=CREAT
        fi

        # commence batch extraction?
        r=PASS
        if [ $s = ERASE -o $s = CREAT ]; then
            # preserve allele ordering
            plink --bfile "$t" --keep "$b" --keep-allele-order \
                  --make-bed --out "$o" &> /dev/null
            r=${PIPESTATUS[0]}
        fi
        echo "$s    $r    $o"   # report, exit on error
        [ $r = "0" ] || [ $r = PASS ] || exit $r
    done
    # rm "$t".*                   # remove extract
done < "$_w/gls"
echo


# merge into one genome for each batch
echo MSG: merge genotype genome-wise:
for b in $(seq -w 000 $[$_q-1])
do
    o="$_o"/$b
    for f in "$_w/div/"$b.*.bed
    do
        echo ${f%.bed}
    done > "$o.lst"

    pass=n
    if [ -s "$o.bed" -a -s "$o.bim" -a -s "$o.fam" ]
    then
        if [ -z $_e ]
        then
            echo -ne EXIST\\t
            pass=y
        else
            echo -ne ERASE\\t
        fi
    else
        echo -ne CREAT\\t
    fi

    # do not pass the batch extraction
    if [ $pass = y ]
    then
        echo -ne PASS\\t
    else
        plink --merge-list "$o.lst" --keep-allele-order --make-bed --out "$o" \
            &>/dev/null
        echo -ne $?\\t
    fi
    echo "$o"
done
echo "result is written to" "$_o"


# reaching this point, the sample division is successful,
# therefore we can safely delete intermediate files.
if [ $retain = y ]
then
    [ $_v ] && echo MSG: "$_w" - retained.
else
    [ $_v ] && echo MSG: "$_w" - clean up.
    rm -rf "$_w" "$_o"/*.{log,lst}
fi
