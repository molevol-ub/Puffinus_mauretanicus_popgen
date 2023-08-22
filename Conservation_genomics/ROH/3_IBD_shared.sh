#!/bin/bash
#$ -cwd
#$ -R y
#$ -e roh_plink.err
#$ -o roh_plink.out
#$ -q h13.q
#$ -pe ompi255h13 12
#$ -V                    #export environment var
#$ -N ROH_plink             #name Job

# List of names
inds=("ALT78" "CZA11" "ILA13" "ILA2" "PORQ" "TZE1" "M8" "M1" "M5" "M20" "M4" "G3" "M2" "M3" "M18" "M13" "M19" "M14" "M12" "M11" "M17" "M10" "M21" "G12" "G9" "G10" "G11" "G14" "G15" "M16" "G4" "M6" "COP1" "LT2")
VCF=/users-d3/jferrer/gizquierdo/TFM/conservation_genomics/ROHs/plink/proves/Puffinus_ROHs.vcf.gz

# Iterate through combinations of ind1 and ind2
for ((i=0; i<${#inds[@]}; i++)); do
    for ((j=i+1; j<${#inds[@]}; j++)); do
        ind1="${inds[i]}"
        ind2="${inds[j]}"
        output_file="${ind2}.${ind1}.intersect.bed"

        # Check if the output file already exists, if not, run the command
        if [ ! -f "$output_file" ]; then
            echo "Running bedtools intersect for $ind1 and $ind2"
            bedtools intersect -a $ind1.bed -b $ind2.bed > "$output_file"
            
            vcftools --gzvcf $VCF --bed $output_file --indv $ind1 --indv $ind2 --recode --out prova.$ind1.$ind2
            
            # Run the custom script that generates hybrid VCF files and calculates ROH
            grep -v "#" prova.$ind1.$ind2.recode.vcf > prova.vcf
            
            while IFS= read -r line; do
                fields=($line)  # Split the line into fields
            
                col10="${fields[9]}"
                col11="${fields[10]}"
            
                if [[ $col10 == "0/1"* && $col11 == "0/1"* ]]; then
                    col10="0/0${col10:3}"
                elif [[ $col10 == "0/1"* && $col11 == "0/0"* ]]; then
                    col10="0/0${col10:3}"
                elif [[ $col10 == "0/1"* && $col11 == "1/1"* ]]; then
                    :
                elif [[ $col10 == "0/0"* && $col11 == "0/0"* ]]; then
                    :
                elif [[ $col10 == "0/0"* && $col11 == "0/1"* ]]; then
                    :
                elif [[ $col10 == "0/0"* && $col11 == "1/1"* ]]; then
                    col10="0/1${col10:3}"
                elif [[ $col10 == "1/1"* && $col11 == "0/0"* ]]; then
                    col10="0/1${col10:3}"
                elif [[ $col10 == "1/1"* && $col11 == "0/1"* ]]; then
                    col10="0/1${col10:3}"
                elif [[ $col10 == "1/1"* && $col11 == "1/1"* ]]; then
                    :
                fi
            
                echo "$col10"
            done < prova.vcf > prova.2.vcf

            grep -v "#" prova.$ind1.$ind2.recode.vcf > prova.vcf
            grep "#" prova.$ind1.$ind2.recode.vcf > $ind1.$ind2.vcf
            cut -f 1,2,3,4,5,6,7,8,9 prova.vcf > prova
            mv prova prova.vcf
            paste prova.vcf prova.2.vcf >> $ind1.$ind2.vcf
            cut -f 1,2,3,4,5,6,7,8,9,10 $ind1.$ind2.vcf >prova.vcf
            mv prova.vcf $ind1.$ind2.vcf
            
            VCF=$ind1.$ind2.vcf
            plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# -indep-pairwise 50 10 0.1 --homozyg-window-snp 100 --homozyg-window-het 1 --homozyg-window-threshold 0.05 --homozyg-snp 25 --homozyg-kb 100 --homozyg-density 50 --homozyg-gap 1000 --homozyg-het 750 --out $ind1.$ind2
            
        else
            echo "Skipping $ind1 and $ind2, output file already exists: $output_file"
        fi
    done
done

echo "Intersect commands completed!"

rm *intersect*

#--------------------------------------------------------------
grep -v "#" prova.$ind1.$ind2.recode.vcf > prova.vcf

while IFS= read -r line; do
    fields=($line)  # Split the line into fields

    col10="${fields[9]}"
    col11="${fields[10]}"

    if [[ $col10 == "0/1"* && $col11 == "0/1"* ]]; then
        col10="0/0${col10:3}"
    elif [[ $col10 == "0/1"* && $col11 == "0/0"* ]]; then
        col10="0/0${col10:3}"
    elif [[ $col10 == "0/1"* && $col11 == "1/1"* ]]; then
        :
    elif [[ $col10 == "0/0"* && $col11 == "0/0"* ]]; then
        :
    elif [[ $col10 == "0/0"* && $col11 == "0/1"* ]]; then
        :
    elif [[ $col10 == "0/0"* && $col11 == "1/1"* ]]; then
        col10="0/1${col10:3}"
    elif [[ $col10 == "1/1"* && $col11 == "0/0"* ]]; then
        col10="0/1${col10:3}"
    elif [[ $col10 == "1/1"* && $col11 == "0/1"* ]]; then
        col10="0/1${col10:3}"
    elif [[ $col10 == "1/1"* && $col11 == "1/1"* ]]; then
        :
    fi

    echo "$col10"
done < prova.vcf > prova.2.vcf

grep -v "#" prova.$ind1.$ind2.recode.vcf > prova.vcf
grep "#" prova.$ind1.$ind2.recode.vcf > $ind1.$ind2.vcf
cut -f 1,2,3,4,5,6,7,8,9 prova.vcf > prova
mv prova prova.vcf
paste prova.vcf prova.2.vcf >> $ind1.$ind2.vcf
cut -f 1,2,3,4,5,6,7,8,9,10 $ind1.$ind2.vcf >prova.vcf
mv prova.vcf $ind1.$ind2.vcf

VCF=$ind1.$ind2.vcf
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# -indep-pairwise 50 10 0.1 --homozyg-window-snp 100 --homozyg-window-het 1 --homozyg-window-threshold 0.05 --homozyg-snp 25 --homozyg-kb 100 --homozyg-density 50 --homozyg-gap 1000 --homozyg-het 750 --out $ind1.$ind2
