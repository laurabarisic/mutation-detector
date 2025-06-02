# 🔬 Pronalazak mutacija pomoću treće generacije sekvenciranja

## 🧭 Opis projekta

Ovaj projekt razvijen je u sklopu kolegija **Bioinformatika 1** na FER-u, s ciljem analize mutacija u DNA sekvencama korištenjem podataka dobivenih metodama sekvenciranja treće generacije.

### Glavni zadaci uključuju:

- Poravnanje očitanja (reads) na referentni genom pomoću alata **minimap2**
- Identifikaciju mutacija poput **supstitucija**, **insercija** i **delecija**
- Usporedbu rezultata s izlazom alata **FreeBayes**
- Generiranje pregledne **CSV datoteke** s detektiranim mutacijama

---

## ⚙️ Instalacija

### ✅ Preduvjeti

Za rad na ovom projektu potrebni su sljedeći alati:

- C++ kompajler (preporuka: `g++`)
- Minimap2
- FreeBayes

### 📦 Instalacija alata

#### Minimap2

```bash
git clone https://github.com/lh3/minimap2
cd minimap2
make

#### FreeBayes

git clone --recursive https://github.com/freebayes/freebayes.git
cd freebayes
make

▶️ Korištenje
🔧 Pokretanje programa

g++ bioinf.cpp -o bioinf
./bioinf

Nakon pokretanja, generira se freebayes_mutations.csv s popisom detektiranih mutacija.
📄 Ulazni i izlazni podaci
Ulaz

    freebayes.vcf: datoteka s varijacijama dobivena iz alata FreeBayes

Izlaz

    freebayes_mutations.csv: CSV tablica u formatu:

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

👩‍🔬 Autori

    Laura Barišić

    Mia Nazor

🎓 Napomena

Ovaj projekt služi isključivo za edukativne svrhe. Potencijalne nadogradnje uključuju:

    Analizu učestalosti mutacija po regijama

    Dodatne vizualizacije rezultata

    Integraciju dodatnih bioinformatičkih alata
