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

### 👣 Koraci pripreme:

    # Koristimo datoteke
    reference="lambda.fasta"
    reads="lambda_simulated_reads.fasta"
    samOutput="lambda.sam"
    bamOutput="lambda.bam"
    sortedBam="lambda_sorted.bam"
    markedBam="lambda_marked.bam"
    csvOutput="lambda_mutated.csv"
    vcfOutput="lambda_freebayes.vcf"
    vcfOutputFiltered="lambda_freebayes_filtered.vcf"  
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

    # Step 6:  Copy sorted BAM to create marked BAM
    cp $sortedBam $markedBam
    
    # Step 7: Index marked BAM file
    samtools index $markedBam
    
    # Step 8: Run FreeBayes for evaluation 
    freebayes -f $reference $markedBam > $vcfOutput

    # Step 9: Run FreeBayes with filters
    freebayes -f $reference --min-alternate-count 3 --min-alternate-fraction 0.4 $markedBam > $vcfOutputFiltered

    # Step 10: Create .csv file from .vcf using converter.cpp
    cd ~/mutation-detector/src
    g++ converter.cpp -o converter
    ./converter lambda_freebayes.vcf lambda_freebayes_mutations.csv
    ./converter lambda_freebayes_filtered.vcf lambda_freebayes_filtered_mutations.csv

    
### ▶️ Pokretanje programa

    cd src
    g++ bioinf.cpp -o bioinf
    *ukoliko pokrećemo program nad lambda podacima:
        ./bioinf lambda.fasta lambda.sam lambda_mutations.csv
    *ukoliko pokrećemo program nad ecoli podacima:
        ./bioinf ecoli.fasta ecoli.sam ecoli_mutations.csv    
    

Nakon pokretanja, generira se.csv datoteka s popisom detektiranih mutacija.

### 📄 Ulazni i izlazni podaci

    Ulaz:
    
        lambda.fasta/ecoli.fasta, lambda.sam/ecoli.sam (ulazne datoteke se nalaze u /data direktoriju projekta)
    
    Izlaz:

        - CSV tablica (lambda_mutations.csv/ecoli_mutations.csv) u formatu:
        
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

        - matching.txt datoteka: sadrži ispis uspoređivanja svake baze referentnog genoma s njom poravnatom bazom očitanja
        - voting.txt datoteka: sadrži ispis broja glasova svih mutacija koje su zabiježene na poravnatim očitanjima za svaku poziciju referentnog genoma

        (izlane datoteke se kreiraju u /data direktoriju projekta)

## 📊 Usporedba i točnost rezultata

    - točnost detektiranih mutacija bodujemo sa accuracy.cpp skriptom u kojoj se dobivene mutacije iz .csv datoteke 
    (čije ime proslijedimo prilikom poziva programa) uspoređuju sa .csv datotekom dobivenom u sklopu projekta koju uzimamo kao referentne rezultate
    
        g++ accuracy.cpp -o accuracy
        ./accuracy lambda_mutations.csv lambda_mutated.csv 

    - aaliza usporedbe rezultata dobivenih mutacija ovog programa, FreeBayesovih rezltata te u sklopu projekta dobivenih referentnih rezultata
    može se dobiti pokretanjem skripte compare_mutations.cpp koja ispisuje:
            
            * koliko ukupno varijanti ima svaka datoteka
            * koliko zajedničkih varijanti ima između svake dvije datoteke
            * koliko varijanti ima svaka datoteka, a da ih nema u drugim datotekama

        g++ compare_mutations.cpp -o compare_mutations
        ./compare_mutations


## 👩‍🔬 Autori

    - Laura Barišić

    - Mia Nazor

## 🎓 Napomena

Ovaj projekt služi isključivo za edukativne svrhe. Potencijalne nadogradnje uključuju:

    - analizu učestalosti mutacija po regijama

    - dodatne vizualizacije rezultata


