#!/bin/bash

######################################### 0.init
function printHelp(){
    echo "Version: 1.0.0"
    echo
    echo
    echo "Required Arguments:"
    echo
    echo "--sample=File                      Sample file, separated by \",\",  (Default: sample.txt)"
    echo "--paired=Boolean                   Whether paired end, Possible values: {true or false} (Default: true)"
    echo
    echo
    echo "Optional Arguments:"
    echo
    echo "--nomodel=Boolean                  It is generally set to --nomodel for histone marks and ATAC, not for TF.  (Default: true)"
    echo "--broad=Boolean                    It is generally set to --broad for histone modifications,"
    echo "                                   but narrow peaks for transcription factors and ATAC-seq data. (Default: true)"
    echo "--species=character                'hs' for human,'mm' for mouse and 'dm' for drosophila. (Default: hs)"
    echo "--suffix=character                 Suffix of fastq files, 'fastq.gz' or 'fq.gz'. (Default: fastq.gz)"
    echo "--split=character                  Rule of how to name paired reads, a.1.fq.gz or a_1.fq.gz, '.' or '_'. (Default: _)"
    echo "--extend=numeric                   If single end reads, reads need to be extended, Default: 200)"
}


for p in "$@"
do
    if [[ "$p" = *"help"* || "$p" = *"-h" ]]; then
        printHelp
        exit 88
    fi

    if [[ "$p" = "--sample"* ]]; then
        # skip argument to prevent propogating and get its value
        SAMPLE=`echo $p | cut -d "=" -f 2`
        #echo $SAMPLE
    fi

    if [[ "$p" = "--paired"* ]]; then
        PAIRED=`echo $p | cut -d "=" -f 2`
        #echo $PAIRED
    fi

    if [[ "$p" = "--suffix"* ]]; then
       suffix=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--split"* ]]; then
      split=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--nomodel"* ]]; then
        NOMODEL=`echo $p | cut -d "=" -f 2` 
    fi

    if [[ "$p" = "--broad"* ]]; then
        BROAD=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--broad-cutoff"* ]]; then
        BROADCUTOFF=`echo $p | cut -d "=" -f 2`
    fi	

    if [[ "$p" = "--qvalue"* ]]; then
        qvalue=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--species"* ]]; then
        species=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--extend"* ]]; then
        extend=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--thread"* ]]; then
        nThread=`echo $p | cut -d "=" -f 2`
    fi
done
######################################### Helper ends here



######################################### 1. Default variant
if [ ! -n "$SAMPLE" ]; then
    SAMPLE="sample.txt"
fi
if [ ! -n "$nThread" ]; then
    nThread=16
fi
if [ ! -n "$suffix" ]; then
    suffix=".fastq.gz"
fi
if [ ! -n "$split" ]; then
    split="_"
fi
if [ ! -n "$PAIRED" ]; then
   PAIRED="true"
fi
if [ ! -n "$species" ]; then
    species="hs"
fi
if [ ! -n "$NOMODEL" ]; then
    NOMODEL="true"
fi
if [ ! -n "$BROAD" ]; then
    BROAD="true"
fi
if [ ! -n "$BROADCUTOFF" ]; then
    BROADCUTOFF=0.1
fi
if [ ! -n "$qvalue" ]; then
    qvalue=0.05
fi
if [ ! -n "$extend" ]; then
    extend=200
fi

if [ "$species" = "hs" ]; then
    echo "species output: "$species
    bowtie2Index='/home/reference/Hsapiens/UCSC/hg38/Sequence/Bowtie2Index/hg38'
    chromSize='/public/software/ucsc-tools/hg38.chrom.sizes'
    wuzhong='hs'
    echo $bowtie2Index
fi
if [ "$species" = "mm" ]; then
    echo "species output: "$species
    bowtie2Index='/home/reference/Mmusculus/UCSC/mm9/Sequence/Bowtie2Index/mm9'
    chromSize='/public/software/ucsc-tools/mm9.chrom.sizes'
    echo $bowtie2Index
    wuzhong='mm'
fi
if [ "$species" = "dm" ]; then
    echo "species output: "$species
    bowtie2Index='/home/reference/Dmelanogaster/UCSC/dm6/Sequence/Bowtie2Index/dm6'
    chromSize='/public/software/ucsc-tools/dm6.chrom.sizes'
    echo $bowtie2Index
    wuzhong='dm'
fi
######################################### Default variant ends here

MACS2_ARGS=""
if [[ "$NOMODEL" = "true" ]]; then
    MACS2_ARGS+="--nomodel "
fi

if [[ "$BROAD" = "true" ]]; then
    MACS2_ARGS+="--broad --broad-cutoff $BROADCUTOFF"
else
    MACS2_ARGS+="-q $qvalue"
fi

######################################### 2.Test sample infomation
if [ ! -e $SAMPLE ];then
    printHelp
    echo "Erro: Sample information file is required, default file is sample.txt, default format is"
    echo "      \"treat,input\", with input, if file name are treat.fq.gz and input.fq.gz"
    echo "      \"treat\", without input, if file name is treat_1.fastq.gz"
    exit 99
else
    sampleLine=`cat $SAMPLE | wc -l`
    if [[ "$sampleLine" = "0" ]]; then
        echo "$SAMPLE is empty!"
                exit 100
        fi
fi
######################################### Test sample infomation  ends here



######################################### Main
echo
echo Start time: `date`
echo Working directory: `pwd`
dir=`pwd`
fastqDir=${dir}/Rawdata
echo Fastq files directory: $fastqDir
ls -lh $fastqDir/*$suffix

for line in $(cat $SAMPLE)
do
    var=(${line//,/ })
    for sample in {${var[0]},${var[1]}}
    do
        if [[ "$PAIRED" = "true" ]]; then
            if [ ! -e $fastqDir/${sample}${split}1${suffix} ];then
                echo "$fastqDir/${sample}${split}1${suffix} do not exist!"
                exit 110
            fi
        else
            if [ ! -e $fastqDir/${sample}${suffix} ];then
                echo "$fastqDir/${sample}${suffix} do not exist!"
                exit 110
            fi            
        fi
    done
done

echo
echo "Input Sample Information:"
printf "%-20s%-20s\n" "Treat" "Input"
for line in $(cat $SAMPLE)
do
    var=(${line//,/ })
    printf "%-20s%-20s\n" ${var[0]} ${var[1]} 
done

echo
echo "Data Processing..."


################################3.mkdir
#if [ ! -d fastqc ];then
#        mkdir -p fastqc
#fi

if [ ! -d trim ];then
    echo "mkdir trim"
    mkdir -p trim
fi
if [ ! -d bam ];then
    echo "mkdir bam"
    mkdir -p bam
fi
if [ ! -d bw ];then
    echo "mkdir bw"
    mkdir -p bw
fi
if [ ! -d peak ];then
    echo "mkdir peak"
    mkdir -p peak
fi

#################################

for line in $(cat sample.txt)
do
    var=(${line//,/ })
    ######################################## Operation for all samples: trimming, mapping and rmDuplicates
    for sample in {${var[0]},${var[1]}}
    do
        echo "Start analysising for ${sample}..."

        TRIM_ARGS=""
        MAP_G_ARGS=""

        if [[ "$PAIRED" = "true" ]]; then
            TRIM_ARGS+="--paired $fastqDir/${sample}${split}1${suffix} $fastqDir/${sample}${split}2${suffix}"
            MAP_G_ARGS+="-1 trim/${sample}${split}1_val_1.fq.gz -2 trim/${sample}${split}2_val_2.fq.gz"

        else
            TRIM_ARGS+="$fastqDir/${sample}${suffix}"
            MAP_G_ARGS+="-U trim/${sample}_trimmed.fq.gz"

        fi

        if ([[ "$TRIM_ARGS" = "--paired"* ]] && [ ! -e trim/${sample}${split}1_val_1.fq.gz ] && [ ! -e trim/${sample}${split}2_val_2.fq.gz ]) || ([[ "$TRIM_ARGS" != "--paired"* ]] && [ ! -e trim/${sample}_trimmed.fq.gz ]); then
            echo "Trimming adaptor and low-quality reads..."
        #    trim_galore $TRIM_ARGS --gzip -o trim -j 8
        fi

        if [ ! -e bam/${sample}.bam ];then
            echo "Mapping to reference genome..."
            echo $bowtie2Index
            bowtie2 -x $bowtie2Index -p $nThread $MAP_G_ARGS | samtools sort -@ $nThread -O bam -o bam/${sample}.bam
        fi


        if [ ! -e bam/${sample}.rmdup.bam ];then
            echo "rmDuplicates..."
            gatk MarkDuplicates -I bam/${sample}.bam --ADD_PG_TAG_TO_READS false --REMOVE_SEQUENCING_DUPLICATES true --CREATE_INDEX true -O bam/${sample}.rmdup.bam -M bam/${sample}.rmdup.matrix.txt
        fi
    done
    ######################################## Trimming, mapping and rmDuplicates end here

    
    ######################################## With or without input, Operation: Creating bw file
    if [ -n "${var[1]}" ]; then
        echo "Treat:${var[0]}  Input:${var[1]}   Creating bw file... "
        # If single end, extending is needed
        if ([[ "$PAIRED" = "true" ]] && [ ! -e bw/${var[0]}.compareInput.bedgraph ]); then
            bamCompare -p $nThread -b1 bam/${var[0]}.rmdup.bam -b2 bam/${var[1]}.rmdup.bam --skipNAs --scaleFactorsMethod readCount --operation subtract --outFileFormat bedgraph -o bw/${var[0]}.compareInput.bedgraph
	fi
        if ([[ "$PAIRED" = "false" ]] && [ ! -e bw/${var[0]}.compareInput.bedgraph ]); then
            bamCompare -p $nThread -b1 bam/${var[0]}.rmdup.bam -b2 bam/${var[1]}.rmdup.bam --skipNAs --scaleFactorsMethod readCount --operation subtract --outFileFormat bedgraph -o bw/${var[0]}.compareInput.bedgraph --extendReads $extend
        fi
	
	if [ ! -e bw/${var[0]}.move0.bedgraph ];then
        	awk '{if($4<0){$4=0};print}' bw/${var[0]}.compareInput.bedgraph >bw/${var[0]}.move0.bedgraph
	fi

	if [ ! -e bw/${var[0]}.rpm.bedgraph ];then
        	totalreads=`awk '{sum+=$4}END{print sum}' bw/${var[0]}.move0.bedgraph`
        	echo 'reads left after remove duplication: '$totalreads
        	awk -v totalreads=$totalreads '{$4=$4/totalreads*1000000;print}' bw/${var[0]}.move0.bedgraph >bw/${var[0]}.rpm.bedgraph
	fi

	if [ ! -e bw/${var[0]}.sort.bedgraph ];then
        	sort -k1,1 -k2,2n bw/${var[0]}.rpm.bedgraph >bw/${var[0]}.sort.bedgraph
	fi

	if [ ! -e bw/${var[0]}.rpm.bw ];then
        	bedGraphToBigWig bw/${var[0]}.sort.bedgraph $chromSize bw/${var[0]}.rpm.bw
	fi
	
	################### call peak
        echo "Treat:${var[0]}  Input:${var[1]}   Calling peak... "
        # If single end, extending is needed
        if ([[ "$PAIRED" = "true" ]] && [ ! -e peak/${var[0]}_peaks.xls ]); then
            macs2 callpeak -t bam/${var[0]}.rmdup.bam -c bam/${var[1]}.rmdup.bam -g $wuzhong -f BAMPE -n ${var[0]} --outdir peak $MACS2_ARGS 2>${var[0]}.masc2.log
	fi

        if ([[ "$PAIRED" = "false" ]] && [ ! -e peak/${var[0]}_peaks.xls ]); then
            macs2 callpeak -t bam/${var[0]}.rmdup.bam -c bam/${var[1]}.rmdup.bam -g $wuzhong -f BAM -n ${var[0]} --outdir peak $MACS2_ARGS --extsize $extend 2>${var[0]}.masc2.log
        fi

    else

        echo "Treat:${var[0]}  No input   Creating bw file... "
        # If single end, extending is needed
        if ([[ "$PAIRED" = "true" ]] && [ ! -e bw/${sample}.CPM.bw ]); then
            bamCoverage -p $nThread -b bam/${sample}.rmdup.bam --skipNAs --normalizeUsing CPM -o bw/${sample}.CPM.bw
        fi

	if ([[ "$PAIRED" = "false" ]] && [ ! -e bw/${sample}.CPM.bw ]); then
            bamCoverage -p $nThread -b bam/${sample}.rmdup.bam --skipNAs --normalizeUsing CPM -o bw/${sample}.CPM.bw --extendReads $extend
        fi
    
        echo "Treat:${var[0]}  No input   Calling peak... "
        # If single end, extending is needed
        if ([[ "$PAIRED" = "true" ]] && [ ! -e peak/${var[0]}_peaks.xls ]); then
            macs2 callpeak -t bam/${var[0]}.rmdup.bam -g $wuzhong -n ${var[0]} -f BAMPE --outdir peak 2>${var[0]}.masc2.log $MACS2_ARGS
	fi

        if ([[ "$PAIRED" = "false" ]] && [ ! -e peak/${var[0]}_peaks.xls ]); then
            macs2 callpeak -t bam/${var[0]}.rmdup.bam -g $wuzhong -n ${var[0]} -f BAM --outdir peak 2>${var[0]}.masc2.log $MACS2_ARGS --extsize $extend
        fi

    fi
    ######################################## Ends here

done

echo End time: `date`
