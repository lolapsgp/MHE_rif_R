cd /fast/AG_Forslund/Lola/Secuencias_INCLIVA_2024/Assembly/Segatella_copri/results/prokka/NR

# Header
echo -e "MAG\tcontigs\tbases\tCDS\trRNA\ttRNA" > prokka_table.txt

# Loop over all MAG folders
for d in */; do
    TXT="$d$(basename $d).txt"
    if [ -f "$TXT" ]; then
        MAG=$(basename $d)
        contigs=$(grep "^contigs:" "$TXT" | awk '{print $2}')
        bases=$(grep "^bases:" "$TXT" | awk '{print $2}')
        CDS=$(grep "^CDS:" "$TXT" | awk '{print $2}')
        rRNA=$(grep "^rRNA:" "$TXT" | awk '{print $2}')
        tRNA=$(grep "^tRNA:" "$TXT" | awk '{print $2}')
        echo -e "${MAG}\t${contigs:-NA}\t${bases:-NA}\t${CDS:-NA}\t${rRNA:-NA}\t${tRNA:-NA}" >> prokka_table.txt
    fi
done
