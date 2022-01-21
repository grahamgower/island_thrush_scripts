# shell functions and settings

# Make exit status of programs propagate through pipes.
set -o pipefail

# So lets hide tmp files elsewhere.
export TMPDIR=/willerslev/scratch/grg/tmp
mkdir -p -m 1775 $TMPDIR

# Make java shit play slightly nicer...
ulimit -n $(ulimit -Hn)
java="java -Xmx8G -XX:ParallelGCThreads=2"
picard=/willerslev/software/picard/picard.jar
gatk=/home/srx907/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
export java picard gatk

# Print error message and exit.
die() {
        echo Error: $@
        exit 1
}
export -f die

# Check file existence and basic bam integrity.
bam_check() {
        if [ $1 == "-i" ]; then
                check_idx=1
                shift
        else
                check_idx=0
        fi
        bam=$1

        if [ ! -f "${bam}" ]; then
                return 1
        fi
        if [ $check_idx == "1" ] && \
	   [ ! -f "${bam%.bam}.bai" -a ! -f "${bam}.bai" ]; then
                return 2
        fi
        if ! samtools quickcheck "${bam}"; then
                return 3
        fi
        return 0
}
export -f bam_check
