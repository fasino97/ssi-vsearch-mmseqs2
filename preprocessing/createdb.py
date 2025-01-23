import subprocess
import os

def create_mmseqs_db(input_fasta, db_dir, db_name, mmseqs_path="mmseqs"):
    try:
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

        db_path = os.path.join(db_dir, db_name)
        create_db_cmd = [
            mmseqs_path, "createdb", input_fasta, db_path
        ]
        
        print(f"Creating MMseqs2 database from {input_fasta}...")
        result = subprocess.run(create_db_cmd, capture_output=True, text=True, check=True)
        print("Database creation output:", result.stdout)
        print("Database creation error (if any):", result.stderr)
        
        print(f"MMseqs2 database created successfully at {db_path}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error creating MMseqs2 database: {e}")
        raise

if __name__ == "__main__":
    fasta_file = "C:/Users/afasi/OneDrive/Desktop/MMSeqs/Thesis/fastafiles/batch1.fasta"
    db_output_dir = "C:/Users/afasi/OneDrive/Desktop/MMSeqs/Thesis/databases"
    mmseqs_path = "C:/Users/afasi/OneDrive/Desktop/MMSeqs/Thesis/mmseqs/mmseqs.bat"
    db_name = "batch1db"  

    create_mmseqs_db(fasta_file, db_output_dir, db_name, mmseqs_path=mmseqs_path)
