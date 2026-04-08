#!/usr/bin/env python3
import argparse
import sys
import os

def check_file_exists(filepath, desc):
    if not os.path.exists(filepath):
        print(f"FAIL: {desc} '{filepath}' does not exist (Current Dir: {os.getcwd()})")
        return False
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--phenotype_file", required=True)
    parser.add_argument("--covariate_file", required=True)
    parser.add_argument("--include_studies", default="")
    parser.add_argument("--include_proteins", default="")
    parser.add_argument("--covariates", default="")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    report_path = os.path.join(args.outdir, "validation_report.txt")
    
    error_found = False
    with open(report_path, "w") as f:
        # 1. Samplesheet
        if not check_file_exists(args.samplesheet, "Samplesheet"):
            f.write(f"FAIL: Samplesheet missing\n")
            sys.exit(1)
            
        studies = []
        with open(args.samplesheet, "r") as s:
            header = s.readline().strip().split(",")
            for line in s:
                row = line.strip().split(",")
                if len(row) > 0:
                    studies.append(dict(zip(header, row)))
                    
        req_studies = [x.strip() for x in args.include_studies.split(",")] if args.include_studies else []
        
        # 2. Check study files
        for s in studies:
            sid = s['study_id']
            if req_studies and sid not in req_studies:
                continue
                
            f.write(f"Checking study {sid}...\n")
            
            # Genetic data detection
            bfile_prefix = s['plink_bfile']
            is_pfile = os.path.exists(bfile_prefix + ".pgen")
            is_bfile = os.path.exists(bfile_prefix + ".bed")
            
            if not (is_pfile or is_bfile):
                f.write(f"FAIL: No genetic data found for {sid} with prefix {bfile_prefix}\n")
                f.write(f"      (Looked for {bfile_prefix}.pgen or .bed)\n")
                error_found = True
            else:
                fmt = "PLINK2 (pfile)" if is_pfile else "PLINK1.9 (bfile)"
                f.write(f"Found {fmt} data for {sid}\n")

        # 3. Phenotype Header
        if not check_file_exists(args.phenotype_file, "Phenotype file"):
            f.write(f"FAIL: Missing phenotype file\n")
            error_found = True
        else:
            with open(args.phenotype_file, "r") as pf:
                line = pf.readline().strip()
                pheno_header = line.split("\t") if "\t" in line else line.split(",")
                
                req_proteins = []
                if args.include_proteins:
                    if os.path.exists(args.include_proteins):
                        with open(args.include_proteins, "r") as f_prots:
                            req_proteins = [x.strip() for x in f_prots if x.strip()]
                    else:
                        req_proteins = [x.strip() for x in args.include_proteins.split(",")]
                
                if req_proteins:
                    missing_proteins = [p for p in req_proteins if p not in pheno_header]
                    if missing_proteins:
                        # Only report first 10 missing to avoid flooding
                        show_missing = missing_proteins[:10]
                        f.write(f"FAIL: Missing requested proteins (showing {len(show_missing)}/{len(missing_proteins)}): {', '.join(show_missing)}\n")
                        error_found = True

        # 4. Covariate Header (Fallback to phenotype file if covariates.txt doesn't exist)
        cov_source = args.covariate_file if os.path.exists(args.covariate_file) else args.phenotype_file
        
        if os.path.exists(cov_source):
            with open(cov_source, "r") as cf:
                line = cf.readline().strip()
                cov_header = line.split("\t") if "\t" in line else line.split(",")
                
                if args.covariates:
                    req_covs = [x.strip() for x in args.covariates.split(",")]
                    missing_covs = [c for c in req_covs if c not in cov_header]
                    if missing_covs:
                        f.write(f"FAIL: Missing requested covariates: {', '.join(missing_covs)}\n")
                        error_found = True
        else:
            f.write(f"FAIL: Covariate source '{cov_source}' not found.\n")
            error_found = True
            cov_header = []

        for s in studies:
            if req_studies and s['study_id'] not in req_studies:
                continue

            # Only check group_column if it's actually specified in samplesheet
            gcol = s.get('group_column', '').strip()
            if not gcol:
                continue

            in_cov = gcol in cov_header
            in_sample = False
            sample_file = s.get('sample_file', '')
            if sample_file and os.path.exists(sample_file):
                if sample_file.endswith((".fam", ".sam")):
                    sample_header = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
                    if gcol in sample_header:
                        in_sample = True
                else:
                    with open(sample_file, "r") as sf:
                        sample_header = sf.readline().strip().split()
                        if gcol in sample_header:
                            in_sample = True

            if not in_cov and not in_sample:
                f.write(f"FAIL: Group column '{gcol}' missing for {s['study_id']} (covariate/sample files)\n")
                error_found = True
        
        if error_found:
            f.write("Validation finished with ERRORS.\n")
        else:
            f.write("Validation PASS.\n")

    if error_found:
        print("Validation finished with ERRORS. Check reports.")
        sys.exit(1)
    else:
        print("Validation PASS.")

if __name__ == "__main__":
    main()
