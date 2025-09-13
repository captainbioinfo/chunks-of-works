# It contains QC (fastqc and multiqc) for both (paired-end and single-end) reads.
import os
import subprocess

# Function to run FastQC on all .fastq files in the input directory
def run_fastqc(raw, fastqc_out1):
    # Check if the output directory exists, create it if it doesn't
    if not os.path.exists(fastqc_out1):
        os.makedirs(fastqc_out1)

    # List all .fastq or .fq files in the raw directory
    file_to_process = [f for f in os.listdir(raw) if f.endswith(".fastq") or f.endswith(".fq")]

    # Check if there are any files to process
    if not file_to_process:
        print("No .fastq or .fq files found in the raw directory.")
        return

    # Loop through all the .fastq files in the raw directory
    for filename in file_to_process:
        file_path = os.path.join(raw, filename)
        print(f"Running FastQC on {filename}...")

        # Run the FastQC command on the file and save the output in fastqc_out1 directory
        command = f"fastqc {file_path} --outdir {fastqc_out1}"
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"FastQC analysis completed for {filename}")
        except subprocess.CalledProcessError as e:
            print(f"Error running FastQC on {filename}: {e}")

# Function to run MultiQC on the FastQC output directory
def run_multiqc(fastqc_out1, multiqc_out):
    # Check if the MultiQC output directory exists, create it if it doesn't
    if not os.path.exists(multiqc_out):
        os.makedirs(multiqc_out)

    print("Running MultiQC...")
    command = f"multiqc {fastqc_out1} -o {multiqc_out}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"MultiQC report generated at {multiqc_out}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MultiQC: {e}")

# Example usage
raw = '/home/farha/raw'  # Replace with your raw input directory
fastqc_out1 = '/home/farha/fastqc_out1'  # FastQC output directory
multiqc_out = '/home/farha/multiqc_report'  # MultiQC output directory

run_fastqc(raw, fastqc_out1)
run_multiqc(fastqc_out1, multiqc_out)
