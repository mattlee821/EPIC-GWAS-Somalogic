#!/usr/bin/env python3
import argparse
import sys
import os

from config import read_sample_file, resolve_column_name

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
        def report(message, *, is_error=False):
            print(message, file=f)
            print(message, file=sys.stderr if is_error else sys.stdout)

        # 1. Samplesheet
        if not check_file_exists(args.samplesheet, "Samplesheet"):
            report("FAIL: Samplesheet missing", is_error=True)
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
                
            report(f"Checking study {sid}...")
            
            # Genetic data detection
            bfile_prefix = s['plink_bfile']
            is_pfile = os.path.exists(bfile_prefix + ".pgen")
            is_bfile = os.path.exists(bfile_prefix + ".bed")
            
            if not (is_pfile or is_bfile):
                report(f"FAIL: No genetic data found for {sid} with prefix {bfile_prefix}", is_error=True)
                report(f"      (Looked for {bfile_prefix}.pgen or .bed)", is_error=True)
                error_found = True
            else:
                fmt = "PLINK2 (pfile)" if is_pfile else "PLINK1.9 (bfile)"
                report(f"Found {fmt} data for {sid}")

        # 3. Phenotype Header
        if not check_file_exists(args.phenotype_file, "Phenotype file"):
            report("FAIL: Missing phenotype file", is_error=True)
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
                        report(
                            f"FAIL: Missing requested proteins (showing {len(show_missing)}/{len(missing_proteins)}): {', '.join(show_missing)}",
                            is_error=True,
                        )
                        error_found = True

        # 4. Covariate Header (Fallback to phenotype file if covariates.txt doesn't exist)
        cov_source = args.covariate_file if os.path.exists(args.covariate_file) else args.phenotype_file
        
        if os.path.exists(cov_source):
            with open(cov_source, "r") as cf:
                line = cf.readline().strip()
                cov_header = line.split("\t") if "\t" in line else line.split(",")
                
                if args.covariates:
                    req_covs = [x.strip() for x in args.covariates.split(",")]
                    missing_covs = [c for c in req_covs if resolve_column_name(cov_header, c) is None]
                    if missing_covs:
                        report(f"FAIL: Missing requested covariates: {', '.join(missing_covs)}", is_error=True)
                        error_found = True
        else:
            report(f"FAIL: Covariate source '{cov_source}' not found.", is_error=True)
            error_found = True
            cov_header = []

        for s in studies:
            if req_studies and s['study_id'] not in req_studies:
                continue

            # Only check group_column if it's actually specified in samplesheet
            gcol = s.get('group_column', '').strip()
            if not gcol:
                continue

            in_cov = resolve_column_name(cov_header, gcol) is not None
            in_sample = False
            sample_file = s.get('sample_file', '')
            if sample_file and os.path.exists(sample_file):
                try:
                    sample_df = read_sample_file(sample_file)
                except Exception as exc:
                    report(f"FAIL: Could not read sample file for {s['study_id']}: {exc}", is_error=True)
                    error_found = True
                else:
                    if resolve_column_name(sample_df.columns, gcol) is not None:
                        in_sample = True

            if not in_cov and not in_sample:
                report(
                    f"FAIL: Group column '{gcol}' missing for {s['study_id']} (covariate/sample files)",
                    is_error=True,
                )
                error_found = True
        
        if error_found:
            report("Validation finished with ERRORS.", is_error=True)
        else:
            report("Validation PASS.")

    if error_found:
        print("Validation finished with ERRORS. Check reports.")
        sys.exit(1)
    else:
        print("Validation PASS.")

if __name__ == "__main__":
    main()
