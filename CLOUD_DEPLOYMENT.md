# Cloud Deployment Guide for RRRM-2 Pipeline

This guide covers running the RRRM-2 kidney transcriptome analysis pipeline on various cloud platforms.

## Quick Start (Any Cloud)

```bash
# 1. Build Docker image
docker build -t rrrm2-pipeline .

# 2. Run with your data mounted
docker run -v /path/to/data:/app/data rrrm2-pipeline --phases 2 3 5 6 7
```

---

## Option 1: Google Cloud Platform (GCP) - Recommended for Beginners

### 1A. Using Google Cloud Shell + Compute Engine

```bash
# SSH into Cloud Shell, then create a VM
gcloud compute instances create rrrm2-vm \
    --zone=us-central1-a \
    --machine-type=n2-standard-8 \
    --boot-disk-size=100GB \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud

# SSH into the VM
gcloud compute ssh rrrm2-vm --zone=us-central1-a
```

### 1B. Setup on the VM

```bash
# Install Docker
sudo apt-get update
sudo apt-get install -y docker.io git
sudo usermod -aG docker $USER
newgrp docker

# Clone repository
git clone https://github.com/ibrahimshahid1/RRRM2_Kidney_Transcriptome.git
cd RRRM2_Kidney_Transcriptome

# Upload your data (from local machine)
# gcloud compute scp --recurse ./data/* rrrm2-vm:~/RRRM2_Kidney_Transcriptome/data/ --zone=us-central1-a

# Build and run
docker build -t rrrm2-pipeline .
docker run -v $(pwd)/data:/app/data rrrm2-pipeline --phases 2 3 5 6 7
```

### 1C. Using Google Cloud Batch (for full production runs)

Create `batch_job.json`:
```json
{
  "taskGroups": [{
    "taskSpec": {
      "runnables": [{
        "container": {
          "imageUri": "gcr.io/YOUR_PROJECT/rrrm2-pipeline:latest",
          "commands": ["--phases", "0", "1", "2", "3", "5", "6", "7"]
        }
      }],
      "computeResource": {
        "cpuMilli": 8000,
        "memoryMib": 32768
      }
    },
    "taskCount": 1
  }],
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}
```

```bash
# Submit batch job
gcloud batch jobs submit rrrm2-full-run \
    --location=us-central1 \
    --config=batch_job.json
```

---

## Option 2: Amazon Web Services (AWS)

### 2A. Using EC2

```bash
# Launch EC2 instance (from AWS CLI)
aws ec2 run-instances \
    --image-id ami-0c55b159cbfafe1f0 \
    --instance-type r5.2xlarge \
    --key-name your-key \
    --security-group-ids sg-xxxxxxxx \
    --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":100}}]'

# SSH into instance
ssh -i your-key.pem ubuntu@<instance-ip>

# Install Docker
sudo apt-get update && sudo apt-get install -y docker.io
sudo usermod -aG docker ubuntu && newgrp docker

# Clone and run
git clone https://github.com/ibrahimshahid1/RRRM2_Kidney_Transcriptome.git
cd RRRM2_Kidney_Transcriptome
docker build -t rrrm2-pipeline .
docker run -v $(pwd)/data:/app/data rrrm2-pipeline
```

### 2B. Using AWS Batch

Create `batch-job-definition.json`:
```json
{
  "jobDefinitionName": "rrrm2-pipeline",
  "type": "container",
  "containerProperties": {
    "image": "YOUR_ECR_REPO/rrrm2-pipeline:latest",
    "vcpus": 8,
    "memory": 32768,
    "command": ["--phases", "0", "1", "2", "3", "5", "6", "7"],
    "mountPoints": [{
      "sourceVolume": "data",
      "containerPath": "/app/data"
    }]
  }
}
```

---

## Option 3: Terra.bio (Best for Genomics)

Terra.bio is built on GCP and is optimized for bioinformatics workflows.

### Setup on Terra

1. Go to https://terra.bio and create a workspace
2. Upload your data to the workspace bucket
3. Create a new notebook or workflow

### Terra Notebook Setup

```python
# In a Terra Jupyter notebook
!git clone https://github.com/ibrahimshahid1/RRRM2_Kidney_Transcriptome.git
%cd RRRM2_Kidney_Transcriptome

# Install dependencies
!pip install -r requirements.txt
!Rscript install_r_packages.R

# Run pipeline
!python src/run_all_phases.py --phases 2 3 5 6 7 --run-id terra_run_001
```

---

## Option 4: GitHub Codespaces (Quickest Setup)

For development and testing, GitHub Codespaces provides instant cloud development:

1. Go to your repo on GitHub
2. Click **Code** → **Codespaces** → **Create codespace on main**
3. Wait for environment to build
4. In terminal:

```bash
pip install -r requirements.txt
Rscript install_r_packages.R
python src/run_all_phases.py --dry-run
```

---

## Option 5: Paperspace Gradient (ML-Optimized)

```bash
# Create a Gradient notebook with Python 3.11 runtime
# In the terminal:
git clone https://github.com/ibrahimshahid1/RRRM2_Kidney_Transcriptome.git
cd RRRM2_Kidney_Transcriptome
pip install -r requirements.txt

# Note: R installation may require building from source on Gradient
```

---

## Data Transfer

### Uploading Data to Cloud

```bash
# GCP (Google Cloud Storage)
gsutil -m cp -r data/ gs://your-bucket/rrrm2-data/

# AWS (S3)
aws s3 sync data/ s3://your-bucket/rrrm2-data/

# From cloud VM, download:
# GCP
gsutil -m cp -r gs://your-bucket/rrrm2-data/ data/
# AWS
aws s3 sync s3://your-bucket/rrrm2-data/ data/
```

---

## Recommended Instance Sizes

| Phase | Minimum | Recommended | Notes |
|-------|---------|-------------|-------|
| Phase 0-1 (R preprocessing) | 4 vCPU, 16GB RAM | 8 vCPU, 32GB RAM | Memory-bound |
| Phase 2-3 (Network + Embedding) | 8 vCPU, 32GB RAM | 16 vCPU, 64GB RAM | CPU-bound |
| Phase 5-6 (Permutation) | 8 vCPU, 32GB RAM | 16 vCPU, 64GB RAM | Many iterations |
| **Full Pipeline** | 8 vCPU, 32GB RAM | **16 vCPU, 64GB RAM** | ~2-4 hours |

### Recommended Instance Types

| Cloud | Instance Type | vCPU | RAM | Cost/hr |
|-------|---------------|------|-----|---------|
| GCP | n2-standard-8 | 8 | 32GB | ~$0.39 |
| GCP | n2-standard-16 | 16 | 64GB | ~$0.78 |
| AWS | r5.2xlarge | 8 | 64GB | ~$0.50 |
| AWS | r5.4xlarge | 16 | 128GB | ~$1.01 |

---

## Cost Estimation

For a full pipeline run with K=2000 permutations/bootstraps:

| Platform | Instance | Time | Cost |
|----------|----------|------|------|
| GCP n2-standard-16 | 16 vCPU, 64GB | ~3 hours | ~$2.50 |
| AWS r5.4xlarge | 16 vCPU, 128GB | ~3 hours | ~$3.00 |
| Spot/Preemptible | Same specs | ~3 hours | ~$0.75-1.00 |

**Tip:** Use spot/preemptible instances for 70-80% cost savings. The pipeline saves checkpoints, so interruptions are recoverable.

---

## Running the Full Pipeline

```bash
# Full production run with all phases and full statistical power
python src/run_all_phases.py \
    --phases 0 1 2 3 5 6 7 \
    --max-genes 2500 \
    --topk 80 \
    --num-seeds 10 \
    --run-id production_v1

# Check the Phase 6 permutation/bootstrap settings in the code
# Default is K=200, but for publication you want K=2000
```

### Modifying Permutation/Bootstrap Counts

Edit `src/statistics/permutation_bootstrap.py` line 92-94:
```python
ap.add_argument("--K_perm", type=int, default=2000)  # Change from 200
ap.add_argument("--B_boot", type=int, default=2000)  # Change from 200
```

Or pass as arguments:
```bash
python -m src.statistics.permutation_bootstrap --K_perm 2000 --B_boot 2000
```

---

## Monitoring and Logging

```bash
# View live logs
docker logs -f <container_id>

# Save logs to file
docker logs <container_id> > pipeline_run.log 2>&1

# Monitor resource usage
docker stats <container_id>
```

---

## Downloading Results

```bash
# From cloud VM to local
scp -r user@cloud-vm:~/RRRM2_Kidney_Transcriptome/data/results ./results_cloud

# Or via cloud storage
gsutil -m cp -r gs://your-bucket/results/ ./results_cloud/  # GCP
aws s3 sync s3://your-bucket/results/ ./results_cloud/      # AWS
```

---

## Troubleshooting

### R Package Installation Fails
```bash
# Install system dependencies first
apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev
```

### Out of Memory
- Reduce `--max-genes` (e.g., 2000 instead of 2500)
- Use a larger instance
- Run phases separately

### Permission Denied on Docker
```bash
sudo chmod 666 /var/run/docker.sock
# Or use: sudo docker ...
```

---

## Quick Reference Commands

```bash
# Build image
docker build -t rrrm2-pipeline .

# Run full pipeline
docker run -v $(pwd)/data:/app/data rrrm2-pipeline --phases 0 1 2 3 5 6 7

# Run specific phases
docker run -v $(pwd)/data:/app/data rrrm2-pipeline --phases 2 3

# Dry run (check commands without executing)
docker run -v $(pwd)/data:/app/data rrrm2-pipeline --dry-run

# Interactive shell
docker run -it -v $(pwd)/data:/app/data rrrm2-pipeline /bin/bash
```
