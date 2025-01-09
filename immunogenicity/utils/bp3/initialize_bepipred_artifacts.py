#!/usr/bin/env python

import os
import requests
import tarfile


def initialize_bepipred30_artifacts(): 

    os.makedirs("BP3Models")
    id_ffnn_url = "https://www.dropbox.com/scl/fi/00dtksfzib6ll9z02isbo/BP3C50IDFFNN.tar.gz?rlkey=vspz0kr5urs9wwb76jmysb9oj&st=b36o52p5&dl=1"
    seqlen_ffnn_url = "https://www.dropbox.com/scl/fi/e3z0zszjq5nz0c5w0j2nz/BP3C50IDSeqLenFFNN.tar.gz?rlkey=ccqnw7ie6br17r6twkap8z4nn&st=h0vsygsr&dl=1"

    bp3_dir = os.path.join(os.path.dirname(__file__), 'BP3Models/')
    bp3_id_tar = os.path.join(bp3_dir, 'temp_id_ffnn_url.tar.gz')
    bp3_seqlen_tar = os.path.join(bp3_dir, 'temp_seqlen_ffnn_url.tar.gz')
    

    id_response = requests.get(id_ffnn_url)
    
    if id_response.status_code == 200:
        with open(bp3_id_tar, 'wb') as f:
            f.write(id_response.content)
        print(f"File downloaded to {bp3_id_tar}")

        # Unpack the tar file
        print(f"Unpacking {bp3_id_tar}...")
        with tarfile.open(bp3_id_tar, 'r') as tar:
            tar.extractall(path=bp3_dir)
        print(f"File unpacked to {bp3_dir}")
        # Optionally, delete the tar file after unpacking (if not needed)
        os.remove(bp3_id_tar)
        print(f"Tar file removed: {bp3_id_tar}")
    else:
        print(f"Failed to download file. HTTP status code: {id_response.status_code}")
    
    seqlen_response = requests.get(seqlen_ffnn_url)
    
    if seqlen_response.status_code == 200:
        with open(bp3_seqlen_tar, 'wb') as f:
            f.write(seqlen_response.content)
        print(f"File downloaded to {bp3_seqlen_tar}")

        # Unpack the tar file
        print(f"Unpacking {bp3_seqlen_tar}...")
        with tarfile.open(bp3_seqlen_tar, 'r') as tar:
            tar.extractall(path=bp3_dir)
        print(f"File unpacked to {bp3_dir}")
        # Optionally, delete the tar file after unpacking (if not needed)
        os.remove(bp3_seqlen_tar)
        print(f"Tar file removed: {bp3_seqlen_tar}")
    else:
        print(f"Failed to download file. HTTP status code: {seqlen_response.status_code}")



if __name__ == "__main__": 
    initialize_bepipred30_artifacts()
    pass