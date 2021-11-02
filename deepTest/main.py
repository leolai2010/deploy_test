from Bio import SeqIO
from fastapi import FastAPI, UploadFile, File
import time
import os 

app = FastAPI()

@app.post("/")
async def root(file: UploadFile = File(...)):
    protein_sequences = {}
    fasta_file = await file.read() 
    dir_path = os.path.dirname(os.path.realpath(__file__))
    filename = f'{dir_path}\\{time.time()}-{file.filename}'
    f = open(f'{filename}', 'wb')
    f.write(fasta_file)
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        protein_sequences[name] = sequence  
    print(protein_sequences)
    return {"Working":"OK"}