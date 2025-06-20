### this script was used to search for homologous elements of ERVs and Errantiviruses

### software used
### Edirect
https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started
### ncbi-blast-2.9.0+-src
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/
### bowtie-1.2.3-linux-x86_64
https://github.com/BenLangmead/bowtie/releases/tag/v1.2.3
### bedtools/2.28.0
https://github.com/arq5x/bedtools2/releases/tag/v2.28.0
### samtools/1.10
https://github.com/samtools/samtools/releases/tag/1.10
### EMBOSS-6.6.0
https://emboss.sourceforge.net/download/
### mafft-7.505
https://mafft.cbrc.jp/alignment/software/source.html
### iqtree-1.6.12-Linux
http://www.iqtree.org/release/v1.6.12
### hmmer 3.3.1
http://hmmer.org/
### HHpred
https://github.com/soedinglab/hh-suite
### colabfold/1.4.0
https://github.com/YoshitakaMo/localcolabfold/releases/tag/v1.4.0

### define the variables
DIRECTORY="<directory where the genome sequence file and the tblastn results are kept>"
CLASS="<species phylum, order, class>"
references="<directory where files of tblastn bait sequences and hmmscan are kept>"
mkdir -p ${DIRECTORY}/${CLASS} # directory where genome assembly files and other results files are kept per CLASS
mkdir -p ${DIRECTORY}/${CLASS}/blastdb # directory where blast databases are kept
mkdir -p ${DIRECTORY}/${CLASS}/tblastn_results # directory where tblastn results are kept
mkdir -p ${DIRECTORY}/${CLASS}/fasta # directory where fasta files of ERV genomic sequences and the putative PBS sequences are kept
mkdir -p ${DIRECTORY}/${CLASS}/indices # directory where bowtie indices are kept
mkdir -p ${DIRECTORY}/${CLASS}/tRNA_PBS # directory where bowtie mapping results for tRNA PBS search are kept
mkdir -p ${DIRECTORY}/${CLASS}/blastn_results # directory where blastn results are kept

### step 1: preparating for the tblastn search

### step 1-1: make a master table per "CLASS"

### list of "CLASS" used for the study
Diptera
Lepidoptera
Coleoptera
Hymenoptera
Hemiptera
Polyneoptera
Palaeoptera
Arthropoda # outside insects
Protostomia # outside Arthropoda
Cnidaria
Ctenophora
Porifera
Deuterostomia # outside Vertebrata

### baits used for tblastn
Dipteran_errantivirus.first-set_ERV.fasta
Dipteran_errantivirus.first-set_GAG.fasta
Dipteran_errantivirus.first-set_POL.fasta
Dipteran_HSV-ENV.fasta
Errantivirus_and_vERV.second-set_ENV.fasta
Errantivirus_and_vERV.second-set_GAG.fasta
Errantivirus_and_vERV.second-set_POL.fasta
Errantivirus_and_vERV.third-set_ENV.fasta
Errantivirus_and_vERV.third-set_GAG.fasta
Errantivirus_and_vERV.third-set_POL.fasta
Errantivirus_and_vERV.fourth-set_ENV.fasta
Errantivirus_and_vERV.fourth-set_GAG.fasta
Errantivirus_and_vERV.fourth-set_POL.fasta
Errantivirus_and_vERV.fifth-set_ENV.fasta
Errantivirus_and_vERV.fifth-set_GAG.fasta
Errantivirus_and_vERV.fifth-set_POL.fasta

### fetch assembly names using edirect
export PATH=${PATH}:${HOME}/edirect
esearch -db Assembly -query "${CLASS} [ORGN]" |\
esummary | xtract -pattern DocumentSummary -element FtpPath_Stats_rpt,AssemblyAccession,Taxid,Organism > ${DIRECTORY}/${CLASS}_genome_assemblies.txt
### <ftp path for download> <genome assembly accession> <taxonomy ID> <species name>

### print out genome assembly and species names
awk '{split($1,a,"/"); print a[10],$4"_"$5}' ${DIRECTORY}/${CLASS}_genome_assemblies.txt > ${DIRECTORY}/${CLASS}_genome_assemblies.names.txt
### <assembly name> <species name>

### printing out the URL for download by wget
awk '{split($1,a,":"); n=split(a[2],b,"/"); print "https://ftp.ncbi.nlm.nih.gov/genomes/all/"b[n-5]"/"b[n-4]"/"b[n-3]"/"b[n-2]"/"b[n-1]"/"b[n-1]"_genomic.fna.gz"
}' ${DIRECTORY}/${CLASS}_genome_assemblies.txt > ${DIRECTORY}/${CLASS}_genome_assemblies.https_URL.txt
### <URL for download>

### fetch the information of reference genomes, assembly level and length from NCBI
### here is the page, for example, for Diptera
https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7147

### assembly a master table that contains genome assembly name, species names, family that the species belongs to, and the length of the assembly
${DIRECTORY}/${CLASS}/${CLASS}_reference-genomes.assembly.species.family.size.txt
### <assembly name> <species name> <family name> <assembly length>
### below is an example
GCA_932526495.1_idBomMajo1.1 Bombylius_major Bombyliidae 304311880


### common part per genome assembly --- START ---

### step 1-2: download the genome sequence using wget
if [[ ! -f "${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna.gz" ]]; then
cat ${DIRECTORY}/${CLASS}_genome_assemblies.https_URL.txt | grep -e ${ASSEMBLY} | while read URL; do
wget --directory-prefix=${DIRECTORY}/${CLASS} $URL
done
fi

### step 1-3: make blast database
if [[ ! -f "${DIRECTORY}/${CLASS}/blastdb/${ASSEMBLY}_NCBI_blast.nsi" ]]; then
zcat ${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna.gz > ${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna
fasta="${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna"
makeblastdb -in ${fasta} -parse_seqids -dbtype nucl -title "${ASSEMBLY}_NCBI_assembly" \
-out ${DIRECTORY}/${CLASS}/blastdb/${ASSEMBLY}_NCBI_blast
#rm ${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna
fi

### step 1-4: running tblastn
### below is an example of tblastn using Dipteran_errantivirus.first-set_${ORF}.fasta as bait sequences "ORF" is ENV, GAG or POL
### bait sequences from vertebrate retroviruses often contain both gag and pol in the same ORFs. We treated them as "POL" and looked for regions that have both "POL" and "ENV" ORFs.
tblastn -db ${DIRECTORY}/${CLASS}/blastdb/${ASSEMBLY}_NCBI_blast -outfmt 10 \
-query ${references}/Dipteran_errantivirus.first-set_${ORF}.fasta \
-out ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_tblastn_results.out.tab
### column headers:
### qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

### step 1-5: extracting high-score hits of ORFs and merge them

### step 1-5-1: for POL, only consider hits with >50 bitscore, merge them allowing a maximum of 300 nt gaps, and extract insertions that are longer than 1500nt with multiple hits
ORF="POL"
bitscore="50"
min_length="1500"
### extract ORFs tblastn hits
awk -v BITSCORE=${bitscore} 'BEGIN{FS=","} {if ($9<$10 && $12>BITSCORE) print $2,$9-1,$10,$1":"$7"-"$8,$5":"$6":"$11":"$12,"+";
else if ($9>$10 && $12>BITSCORE) print $2,$10-1,$9,$1":"$7"-"$8,$5":"$6":"$11":"$12,"-"
}' ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_tblastn_results.out.tab |\
tr ' ' '\t' | sort -k1,1 -k2,2n > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.bed

bedtools merge -d 300 -s -c 5,6 -o count,distinct -i ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.bed |\
awk -v LENGTH=${min_length} '{if ($3-$2>LENGTH && $4>1) print $1,$2,$3,$4,".",$5}' | tr ' ' '\t' > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.high-score.bed

### step 1-5-2: for GAG and ENV, only consider hits with >30 bitscore, merge them allowing a maximum of 300 nt gaps, and extract insertions that are longer than 300nt with multiple hits
ORF="GAG"
ORF="ENV"
bitscore="30"
min_length="300"
### extract ORFs tblastn hits
awk -v BITSCORE=${bitscore} 'BEGIN{FS=","} {if ($9<$10 && $12>BITSCORE) print $2,$9-1,$10,$1":"$7"-"$8,$5":"$6":"$11":"$12,"+";
else if ($9>$10 && $12>BITSCORE) print $2,$10-1,$9,$1":"$7"-"$8,$5":"$6":"$11":"$12,"-"
}' ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_tblastn_results.out.tab |\
tr ' ' '\t' | sort -k1,1 -k2,2n > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.bed

bedtools merge -d 300 -s -c 5,6 -o count,distinct -i ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.bed |\
awk -v LENGTH=${min_length} '{if ($3-$2>LENGTH && $4>1) print $1,$2,$3,$4,".",$5}' | tr ' ' '\t' > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.high-score.bed

### step 1-5-3: identify insertions that carry all three ORFs of GAG, POL and ENV
### allowing gaps of 500 nt between ORFs for insect elements and 1000 nt for all the other elements
gap="500"
### We set 4000 and 10000 nt for minimum and maximum insertion size, respectively, when all three ORFs of GAG, POL and ENV are merged.
### We set 3000 and 10000 nt for minimum and maximum insertion size, respectively, when POL and ENV are merged.
min_size="4000"
max_size="10000"
awk -v GAP=${gap} '{if($6=="+") print $1,$2,$3+GAP,$4,$5,$6; else print $1,$2-GAP,$3,$4,$5,$6}' ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_GAG_preset_tblastn_results.out.high-score.bed |\
awk '{if($2<0) $2=0; print}' | tr ' ' '\t' |\
bedtools intersect -s -wao -a - -b ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_POL_preset_tblastn_results.out.high-score.bed |\
awk '{if($6=="+" && $NF!=0) print $1,$2,$9+500,$4":"$10,".",$6; else if($6=="-" && $NF!=0) print $1,$8-500,$3,$4":"$10,".",$6}' |\
awk '{if($2<0) $2=0; print}' | tr ' ' '\t' |\
bedtools intersect -s -wao -a - -b ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_ENV_preset_tblastn_results.out.high-score.bed |\
awk '{if($6=="+" && $NF!=0) print $1,$2,$9,$4":"$10,".",$6; else if($6=="-" && $NF!=0) print $1,$8,$3,$4":"$10,".",$6}' |\
awk -v MIN_SIZE=${min_size} -v MAX_SIZE=${max_size} '{if($3-$2>MIN_SIZE && $3-$2<MAX_SIZE) print $1,$2-300,$3+300,$1":"$2":"$3":"$4,".",$6}' |\
awk '{if($2<0) $2=0; print}' | tr ' ' '\t' > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset.bed
### We extended the genomic coordinate of the insertions to upstream and downstream to predict ORFs. See methods for details.

### step 1-6: obtain amino acid sequences of uninterrupted ORFs
### intersect the high-score ORFs with the predicted GAG-POL-ENV insertions
bedtools intersect -wo -s -a ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset.bed \
-b ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_errantivirus.first-set_${ORF}_preset_tblastn_results.out.high-score.bed |\
awk -v ORF=${ORF} '{print $7,$8,$9,$4":"ORF,".",$12}' | tr ' ' '\t' |\
bedtools getfasta -name -s -fi ${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna -bed - | fasta_formatter - -t | tr ':' '@' |\
awk '{print ">"$1"\n"toupper($2)}' > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset_${ORF}.fasta

### use transeq from EMBOSS to convert the nucleotide sequence to an amino acid sequence
### uninterrupted ORFs without frame shifts and with fewer than five stop codons (marked by *) are taken for further analyses
transeq ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset_${ORF}.fasta ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset_${ORF}.pep2
cat ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset_${ORF}.pep2 | tr '@' ':' | fasta_formatter - -t |\
awk '{if (gsub("*","*",$2)<5) print ">"$1"\n"$2}' > ${DIRECTORY}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_Dipteran_GAG-POL-ENV_preset_${ORF}.pep

### we then used Clustal omega (https://www.ebi.ac.uk/Tools/msa/clustalo/) to align predicted peptide sequences.

### common part per genome assembly --- END ---


### phylogenetic analysis
### mafft-7.505 to make multiple sequence alignements
mafft --auto ${input_FASTA} > ${mafft_output_alignment_file}

### iqtree-1.6.12-Linux to generate phylogenetic trees
iqtree -s ${mafft_output_alignment_file} -m rtREV+R4 -bb 1000 -nt 4


### running hmm scan using hmmer 3.3.1
### make a database for GyDB entries
### HMM profiles of retrotransposons/retroviruses ORFs were downloaded from https://gydb.org/index.php?title=Collection_HMM
mkdir -p ${references}/GyDB_collection/database
mkdir -p ${references}/GyDB_collection/profiles
cat ${references}/GyDB_collection/profiles/*.hmm > ${references}/GyDB_collection/database/GyDB
hmmpress ${references}/GyDB_collection/database/GyDB

### run a hmmscan
hmmscan -o ${output_hhr_file} ${references}/GyDB_collection/database/GyDB ${query_fasta_file}


### running HHpred
mkdir -p ${references}/HHpred
### step 1: make a multiple alingment (a3m file) using the database UniRef30_2022_02 (https://uniclust.mmseqs.com/)
hhblits -i ${input_fasta} -d ${references}/HHpred/UniRef30_2022_02/UniRef30_2022_02 -oa3m ${output_a3m_from_hhblits} -n 1

### step 2: search for homologous profiles from pdb70 database (https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)
### same parameter as HHpred web search (https://toolkit.tuebingen.mpg.de/tools/hhpred)
hhsearch -i ${output_a3m_from_hhblits} \
-o ${output_hhr_of_hhsearch} \
-d ${references}/HHpred/pdb70 \
-p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 \
-contxt ${references}/HHpred/context_data.crf


### identifying tRNA primer binding sequences (PBS) --- START ---
### this was done per elements

### making a blast database
input_fasta="${DIRECTORY}/${CLASS}/fasta/${ERV_NAME}.fasta"
makeblastdb -in ${input_fasta} -parse_seqids -dbtype nucl -title "${ERV_NAME}" \
-out ${DIRECTORY}/${CLASS}/blastdb/${ERV_NAME}_blast

### run a blastn search against its own sequence
blastn -db ${DIRECTORY}/${CLASS}/blastdb/${ERV_NAME}_blast -outfmt 10 -query ${fasta} \
-out ${DIRECTORY}/${CLASS}/blastn_results/${ERV_NAME}_blastn_results.out.tab
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

### take hits of a bitscore between 200 and 5000 in order to select hits of long sequence, but not the whole sequence itself
### we only allowed two LTRs whose distances are longer than 5000 and shorter than 14000 bp.
awk 'BEGIN{FS=","} {if( ($12>200 && $12<5000 && $3>99) && ($9-$7>5000 || $7-$9>5000) && ($9-$7<14000 || $7-$9<14000) ) print
}' ${DIRECTORY}/${CLASS}/blastn_results/${ERV_NAME}_blastn_results.out.tab > ${DIRECTORY}/${CLASS}/blastn_results/${ERV_NAME}_putative-LTRs_blastn_results.out

### fetch -50 to +50 of the 3prime end of the 5prime LTRs
awk '{if($9>$7) print $2,$8-50,$8+50,$2":"$7-1":"$8":"$9-1":"$10,". +"
}' ${DIRECTORY}/${CLASS}/blastn_results/${ERV_NAME}_putative-LTRs_blastn_results.out | tr ' ' '\t' |\
### make sure to select one entry per element (this process often excluded LTRs that gave two split blastn hits)
sort -k2,2nr | awk '!seen[$1]++' |\
bedtools getfasta -name -s -fi ${input_fasta} -bed - |\
fasta_formatter - -t | awk '{print ">"$1"\n"toupper($2)}' | fastx_reverse_complement > ${DIRECTORY}/${CLASS}/fasta/${ERV_NAME}_putative-LTRs_blastn_results.putative-PBS.revComp.fasta
### <NAME>:<3prime_LTR_START>:<3prime_LTR_END>:<5prime_LTR_START>:<5prime_LTR_END> in 0 based coordinate

### build bowtie index for the putative-PBS regions
bowtie-build ${DIRECTORY}/${CLASS}/fasta/${ERV_NAME}_putative-LTRs_blastn_results.putative-PBS.revComp.fasta ${DIRECTORY}/${CLASS}/fasta/${ERV_NAME}_putative-LTRs_blastn_results.putative-PBS.revComp

### run bowtie to identify tRNA PBS allowing 0 to 2 mismatches
PBS_index="${DIRECTORY}/${CLASS}/fasta/${ERV_NAME}_putative-LTRs_blastn_results.putative-PBS.revComp"
### the most 12 nucleotides of tRNA sequences were collected per CLASS, see methods for details: ${DIRECTORY}/${CLASS}/${CLASS}_tRNAs.last12nt.fasta
for MMs in 0 1 2; do
bowtie -f -v ${MMs} -a -S ${PBS_index} ${DIRECTORY}/${CLASS}/${CLASS}_tRNAs.last12nt.fasta |\
samtools view -bS - | \
bamToBed -i - | awk -v MMs=${MMs} '{if($6=="+") print $0,MMs}' >> ${DIRECTORY}/${CLASS}/tRNA_PBS/${ERV_NAME}_putative-PBS.last12nt.bed
done

### identifying tRNA primer binding sequences (PBS) --- END ---


### measuring the copy number of ERVs in their respective genome --- START ---
ORFs_genomic_fasta="${DIRECTORY}/${CLASS}/fasta/${ERV_NAME}_ORFs.fasta" ### This fasta file contains the genome sequence that covers the GAG, POL and ENV open reading frames. Also note that it does not include LTRs.
### blastn of ORFs against respective genomes
blastn -db ${DIRECTORY}/${CLASS}/blastdb/${ASSEMBLY}_NCBI_blast -outfmt 10 \
-query ${ORFs_genomic_fasta} \
-out ${DIRECTORY}/${CLASS}/blastn_results/${ASSEMBLY}_${SPECIES}_${ERV_NAME}_ORFs.blastn_results.out.tab

### measure the copy number of elements with >98% identity of >98% region
SIZE=`(fasta_formatter -i ${ORFs_genomic_fasta} -t | awk '{print length($2)}')`
cat ${DIRECTORY}/${CLASS}/blastn_results/${ASSEMBLY}_${SPECIES}_${ERV_NAME}_ORFs.blastn_results.out.tab | tr ',' ' ' |\
awk -v SIZE=${SIZE} -v ERV_NAME=${ERV_NAME} '{if($3>98 && $4>SIZE*0.98) count++} END {print ERV_NAME,count}' >> ${DIRECTORY}/${CLASS}/${CLASS}.copy_number.list

### measuring the copy number of ERVs in their respective genome --- END ---


### Alphafold2 --- START ---

### monomer structure prediction: used for the POL bridge region and GP1
### database search and protein sequence clustering
colabfold_search --db-load-mode 1 --threads 12 ${input_fasta} $COLABFOLDDIR/database ${output_directory}
### colabfold_search generates a multiple sequence alignment in the a3m file format.

### structural prediction
colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax ${a3m_file} ${predicted_dir}
### colabfold_batch generates pdb files as outputs.

### trimer structure prediction: used for the F and HSV/gB Ectodomains
### step 1. database search and protein sequence clustering
colabfold_search --db-load-mode 1 --threads 12 ${input_fasta} $COLABFOLDDIR/database ${output_directory}

### step 2. merging multiple a3m files.
cd ${output_directory}
cp ${output_directory}/0.a3m ${output_directory}/1.a3m
cp ${output_directory}/0.a3m ${output_directory}/2.a3m

$COLABFOLDDIR/a3m2multi.sh # this code merges three a3m files into one file

### step 3. structural prediction for the multier
colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax --model-type AlphaFold2-multimer-v2 ${a3m_file} ${predicted_dir}
### colabfold_batch generates pdb files as outputs.

### Alphafold2 --- END ---
