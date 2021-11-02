##!/usr/bin/env python
# Import necessary libraries
from Bio import SeqIO
import re
import numpy as np
import pandas as pd
from fastapi import FastAPI
from typing import Optional

# Function for parsing sequence in FASTA file into a Python dictionary relying on BioPython library
def parseFASTA(sequenceFile):
    protein_sequences = {}
    try:
        fasta_sequences = SeqIO.parse(open(sequenceFile), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            protein_sequences[name] = sequence  
        return protein_sequences
    except EOFError:
        print("Unable to parse FASTA file") 
        return False

# Function for parsing peptide names from TSV file into a Python dictionary
def parseTSV(inputFile):
    peptide_unimod_sequences = []
    peptide_sequences = {}    
    position_regex = "(\(UniMod:\d*\))"
    try:
        with open(inputFile) as tsv_file:
            next(tsv_file)    
            for ptm_unimod in tsv_file:
                peptide_unimod = re.findall(position_regex, ptm_unimod)
                if peptide_unimod:
                    peptide_unimod_sequences.append(ptm_unimod)
        peptide_unimod_sequences = list(dict.fromkeys(peptide_unimod_sequences))
        for ptm in peptide_unimod_sequences:
            ptm = ptm.strip()
            peptide_unimod = re.findall(position_regex, ptm)
            if peptide_unimod:
                peptide_name = re.sub(position_regex, '', ptm)
                if len(peptide_unimod) > 1:
                    temp_ptm = ptm
                    peptide_position = []
                    for unimod_pos in peptide_unimod:
                        multiposition = re.search(unimod_pos, temp_ptm)
                        peptide_position.append(multiposition.start() - 1)
                        temp_ptm = re.sub("\("+multiposition.group()+"\)", '', temp_ptm, count=1)
                else:
                    peptide_position = [re.search(peptide_unimod[0], ptm).start() - 1]
                if peptide_name in peptide_sequences.keys():
                    peptide_sequences[peptide_name].update({ptm:peptide_position})
                else:
                    peptide_sequences[peptide_name] = {ptm:peptide_position}
        return peptide_sequences
    except EOFError:
        print("Unable to parse Input file") 
        return False

# Function for parsing peptide names from TSV file into a Python dictionary
def sequencePairing(sequenceFile, inputFile):
    try:
        protein_sequence = parseFASTA(sequenceFile)
        peptide_sequence = parseTSV(inputFile)
        if protein_sequence == False or peptide_sequence == False:
            raise Exception
        ppp_sequence = {}
        for protein_name, sequence in protein_sequence.items():
            for peptide_name, position in peptide_sequence.items():
                matching = re.search(peptide_name, sequence)
                if matching is not None:
                    for peptide_vSequences, peptide_vPositions in position.items():
                        protein_positions = []
                        for peptide_vPositions in peptide_vPositions:
                            protein_position = peptide_vPositions + matching.start()
                            protein_positions.append(protein_position)
                        protein_positions = ', '.join(str(str_pos) for str_pos in protein_positions)
                        if protein_name in ppp_sequence.keys():
                            ppp_sequence[protein_name].update({peptide_vSequences:protein_positions})
                        else:
                            ppp_sequence[protein_name] = {peptide_vSequences:protein_positions}
        return ppp_sequence
    except EOFError:
        print("Unable to pair protein and peptide sequences") 
        return False


def outputTable():
    sequenceFile = "proteins.fasta"
    inputFile = "input.tsv"
    data = sequencePairing(sequenceFile, inputFile)
    df = pd.DataFrame.from_dict({(i, j): (i, j, k) for i in data.keys() for j, k in data[i].items() }, orient='index')
    df.reset_index().drop(columns=['index'])
    json_package = df.to_json(index=False, orient='split')
    return json_package

app = FastAPI()

@app.get("/")
async def main():
    content = """
        <body>
            <form action="/uploadfile/" method="post">
                <input name="files" type="file">
                <input type="submit">
            </form>
        </body>
    """
    return HTMLResponse(content=content)

@app.post("/uploadfile/")
async def create_upload_file(file: UploadFile = File(...)):
    return {"filename": file.filename}