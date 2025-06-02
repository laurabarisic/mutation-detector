# 🔬 Pronalazak mutacija pomoću treće generacije sekvenciranja

## Opis projekta

Ovaj projekt razvijen je u sklopu kolegija **Bioinformatika 1** na FER-u, s ciljem analize mutacija u DNA sekvencama korištenjem podataka dobivenih metodama sekvenciranja treće generacije.

### Glavni zadaci uključuju:

- poravnanje očitanja (reads) na referentni genom pomoću alata **minimap2**
- identifikaciju mutacija poput **supstitucija**, **insercija** i **delecija**
- usporedbu rezultata s izlazom alata **FreeBayes** i datotekom lambda_mutated.csv
- generiranje pregledne **CSV datoteke** s detektiranim mutacijama

---

## ⚙️ Instalacija

### ✅ Preduvjeti

Za rad na ovom projektu potrebni su sljedeći alati:

- C++ kompajler (preporuka: `g++`)
- Minimap2
- FreeBayes

### 📦 Instalacija alata

#### Minimap2

    git clone https://github.com/lh3/minimap2
    cd minimap2
    make```bash

#### FreeBayes

    git clone --recursive https://github.com/freebayes/freebayes.git
    cd freebayes
    make

### 👣Koraci:
    # Koristimo datoteke
    reference="lambda.fasta"
    reads="lambda_simulated_reads.fasta"
    samOutput="lambda.sam"
    bamOutput="lambda.bam"
    sortedBam="lambda_sorted.bam"
    markedBam="lambda_marked.bam"
    csvOutput="lambda_mutated.csv"
    vcfOutput="lambda_freebayes.vcf"
    freebayesOutput="lambda_freebayes_mutations.csv"
    index="lambda.fasta.fai"
    
    # Step 1: Run Minimap2
    minimap2 -ax map-ont $reference $reads > $samOutput
    
    # Step 2: Call detector to detect mutations
    ./bioinf
    
    # Step 3: Generate index for the reference genome
    samtools faidx $reference
    
    # Step 4: Convert SAM to BAM
    samtools view -bS $samOutput > $bamOutput
    
    # Step 5: Sort BAM file
    samtools sort -o $sortedBam $bamOutput
    
    # Step 6: Index marked BAM file
    samtools index $markedBam
    
    # Step 7: Run FreeBayes for evaluation 
    freebayes -f $reference $markedBam > $vcfOutput

### ▶️ Pokretanje programa

    g++ bioinf.cpp -o bioinf
    ./bioinf

Nakon pokretanja, generira se mutations.csv s popisom detektiranih mutacija.
### 📄 Ulazni i izlazni podaci
    Ulaz:
    
        lambda.fasta, lambda.sam
    
    Izlaz:

       CSV tablica u formatu:
        
        Position,Type,ALT
        X,261,G
        X,627,T
        X,726,A
        D,1043,-
        ...

    Legenda tipova:
    
        X – supstitucija (zamjena jedne baze drugom)
    
        I – insercija (umetanje baze)
    
        D – delecija (brisanje baze)

## 👩‍🔬 Autori

    - Laura Barišić

    - Mia Nazor

## 🎓 Napomena

Ovaj projekt služi isključivo za edukativne svrhe. Potencijalne nadogradnje uključuju:

    - analizu učestalosti mutacija po regijama

    - dodatne vizualizacije rezultata


