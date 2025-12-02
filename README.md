# Casanovo meets Peptide-PRISM

---
A Snakemake-based workflow for de novo sequencing of MS peptidomics data with [Casanovo](https://github.com/Noble-Lab/casanovo) and subsequent annotation of the peptides to the reference genome with [Peptide-PRISM](https://pubmed.ncbi.nlm.nih.gov/32561536/)
---

If you have Casanovo and Peptide-PRISM already installed and you don't want to use Snakemake and Singularity you can simply download workflow/shell/format_conversion_for_PRISM.py and apply it to your Casanovo output. This just an adaptor script to convert the file to the [PEAKS](https://pubmed.ncbi.nlm.nih.gov/14558135/)-like format.

---

## 🛠️ Installation
1. Install snakemake:
   ```bash
   conda create -n snakemake_env
   conda activate snakemake_env
   conda install bioconda::snakemake
   ```
2. Install singularity according to the tutorial [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).
3. Clone the repository:
   ```bash
   git clone https://github.com/mkkoshkina/wf_denovo_sequencing.git
   cd wf_denovo_sequencing
   ```
4. Build the images:
   ```bash
   cd workflow/containers
   sudo sungularity build image/python.3.12.2.debian.sif def/python.3.12.2.debian.def
   sudo sungularity build image/python.3.10.15.gpu.ubuntu.sif def/python.3.10.15.gpu.ubuntu.def
   sudo sungularity build image/peptide.prism.1.1.0.ubuntu.env.sif def/peptide.prism.1.1.0.ubuntu.env.def
   cd ../../
   ```
   The python.3.10.15.gpu.ubuntu.sif image for Casanovo execution is built with cuda:12.4.0. This may need to be adjusted depending on your GPU.
5. Download Peptide-PRISM from [here](https://www.uni-wuerzburg.de/sft/erfindungen-patente-und-lizenzen-jmu-und-ukw/download-software-for-scientific-purposes/). Extract the Peptide-PRISM.zip to your workflow/shell folder:
   ```bash
   unzip ~/Downloads/Peptide-PRISM.zip -d workflow/shell/
   ```
6. Download a reference genome and put it in a folder in genomes. For example:
   ```bash
   cd genomes
   mkdir GRCh38.p14
   cd GRCh38.p14
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
   gzip -d gencode.v47.annotation.gtf.gz
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.p14.genome.fa.gz
   gzip -d GRCh38.p14.genome.fa.gz
   cd ../../
   ```
7. Index the genome:
   ```bash
   singularity run workflow/containers/image/peptide.prism.1.1.0.ubuntu.env.sif
   java -jar workflow/shell/Peptide-PRISM/IndexGenome.jar -s genomes/GRCh38.p14/GRCh38.p14.genome.fa -a genomes/GRCh38.p14/gencode.v47.annotation.gtf -nomapping -o genomes/GRCh38.p14/GRCh38.p14
   exit
   ```
8. Download the Casanovo model weights from [here](https://github.com/Noble-Lab/casanovo/releases). For example:
   ```bash
   cd model
   wget https://github.com/Noble-Lab/casanovo/releases/download/v4.0.0/casanovo_nontryptic.ckpt
   cd ..
   ```
   Modify model/config.yaml if you wish. At least increase/decrease predict_batch_size according to the VRAM size you have available.
9. Modify workflow/profiles/example/config.yaml and workflow/profiles/example/todo.py according to the instructions in them.
10. Run:
    ```bash
    snakemake --profile workflow/profiles/example
    ```

   

