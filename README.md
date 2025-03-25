# AB1 Sequence Annotator

## Overview
AB1 Sequence Annotator is a web application that helps researchers process AB1 sequencing files, convert them to FASTA, extract the best protein frames, and annotate protein sequences using different immunoglobulin numbering schemes.

## Features
- Upload AB1 sequencing files
- Convert AB1 to FASTA format
- Extract best protein frames
- Annotate protein sequences using Chothia or Kabat schemes
- Download processed files

## Prerequisites
- Python 3.7+
- pip

## Installation
1. Clone the repository
```bash
git clone https://github.com/HaroonMH/AB1-Sequence-Annotator.git
cd AB1-Sequence-Annotator
```

2. Create a virtual environment
```bash
python -m venv venv
source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
```

3. Install dependencies
```bash
pip install -r requirements.txt
```

## Running the Application
```bash
python app.py
```

Navigate to `http://localhost:5000` in your web browser.

## Usage
1. Select an AB1 file
2. Choose annotation scheme (Chothia or Kabat)
3. Click "Process File"
4. View results and download processed files
