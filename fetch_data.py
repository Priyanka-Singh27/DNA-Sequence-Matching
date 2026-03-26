import os
import urllib.request
import time

FILES = {
    "human/hbb.fasta": "NM_000518",
    "human/brca1.fasta": "NM_007294",
    "human/opn1lw.fasta": "NM_020061",
    "mitochondrial/human_mt.fasta": "NC_012920",
    "mitochondrial/neanderthal_mt.fasta": "NC_011137",
    "mitochondrial/chimp_mt.fasta": "NC_001643",
    "mitochondrial/gorilla_mt.fasta": "NC_011120",
    "mitochondrial/mouse_mt.fasta": "NC_005089",
    "mitochondrial/zebrafish_mt.fasta": "NC_002333",
    "pathogens/influenza_h1n1.fasta": "NC_026433",
    "pathogens/sars_cov2.fasta": "NC_045512",
    "pathogens/ecoli_k12.fasta": "NC_000913",
}

def fetch_fasta(accession, relative_path):
    # Base URL for NCBI E-utilities
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text"
    url = base_url.format(accession)
    
    # E. coli restriction: strictly first 50,000 bp
    if accession == "NC_000913":
        url += "&seq_start=1&seq_stop=50000"
        
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "data", relative_path)
    
    print(f"Downloading {accession} -> {filepath} ...")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # To prevent blocks, NCBI limits requests without API keys to ~3/sec. 
    # Just space them out gently.
    time.sleep(1)
    
    req = urllib.request.Request(url, headers={'User-Agent': 'GenomeX-Tool/1.0'})
    try:
        with urllib.request.urlopen(req) as response:
            content = response.read().decode('utf-8')
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
        print(" [OK]")
        return True
    except Exception as e:
        print(f" [FAILED] ({e})")
        return False

if __name__ == "__main__":
    print("GenomeX - Fetching Sequence Databases")
    print("=" * 40)
    success = 0
    for path, acc in FILES.items():
        if fetch_fasta(acc, path):
            success += 1
            
    print(f"\nDone. Successfully fetched {success}/{len(FILES)} databases.")
