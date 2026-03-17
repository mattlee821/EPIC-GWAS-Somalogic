#!/bin/bash

# Environment Variables
export DIR_IN="/data/Epic/subprojects/Genetics/sources/Gwas"
export DIR_OUT="/data/Epic/subprojects/Genetics/work/data"

echo "Starting Manual Documented Migration..."

# -----------------------------------------------------------------------------
# 1. BREAST
# -----------------------------------------------------------------------------
echo "Processing: Breast"

# Brea_01_Erneg
mkdir -p "$DIR_OUT/Breast/Brea_01_Erneg/"
rsync -avP "$DIR_IN/Breast/Brea_01_Erneg/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Breast/Brea_01_Erneg/" 2>/dev/null

# Brea_02_Onco
mkdir -p "$DIR_OUT/Breast/Brea_02_Onco/"
rsync -avP "$DIR_IN/Breast/Brea_02_Onco/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Breast/Brea_02_Onco/" 2>/dev/null

# -----------------------------------------------------------------------------
# 2. COLONRECTUM - NOTE: DATA IS MISSING
# -----------------------------------------------------------------------------
echo "Processing: Colonrectum"

# Colo_01_GECCO
mkdir -p "$DIR_OUT/Colonrectum/Colo_01_GECCO/"
rsync -avP "$DIR_IN/Colonrectum/Colo_01_GECCO/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Colonrectum/Colo_01_GECCO/" 2>/dev/null

# -----------------------------------------------------------------------------
# 3. CORPUS_UTERI
# -----------------------------------------------------------------------------
echo "Processing: Corpus_Uteri"

# Corp_01_Endo
mkdir -p "$DIR_OUT/Corpus_Uteri/Corp_01_Endo/"
rsync -avP "$DIR_IN/Corpus_Uteri/GSA2022_922_025_V3/zCall_data/my.report.txt" "$DIR_OUT/Corpus_Uteri/Corp_01_Endo/" 2>/dev/null

# -----------------------------------------------------------------------------
# 4. EPIC_CVD
# -----------------------------------------------------------------------------
echo "Processing: Epic_Cvd"

# Ecvd_01
mkdir -p "$DIR_OUT/Epic_Cvd/Ecvd_01/"
rsync -avP "$DIR_IN/Epic_Cvd/Ecvd_01/Data_Received/"*.{bed,bim,fam,txt} "$DIR_OUT/Epic_Cvd/Ecvd_01/" 2>/dev/null

# Ecvd_02
mkdir -p "$DIR_OUT/Epic_Cvd/Ecvd_02/"
rsync -avP "$DIR_IN/Epic_Cvd/Ecvd_02/Data_Received/"*.{bed,bim,fam,txt} "$DIR_OUT/Epic_Cvd/Ecvd_02/" 2>/dev/null

# Ecvd_03
mkdir -p "$DIR_OUT/Epic_Cvd/Ecvd_03/"
rsync -avP "$DIR_IN/Epic_Cvd/Ecvd_03/Data_Received/"*.{bed,bim,fam,txt} "$DIR_OUT/Epic_Cvd/Ecvd_03/" 2>/dev/null

# -----------------------------------------------------------------------------
# 5. GALLBLADDER
# -----------------------------------------------------------------------------
echo "Processing: Gallbladder"

# Glbd_01
mkdir -p "$DIR_OUT/Gallbladder/Glbd_01/"
rsync -avP "$DIR_IN/Gallbladder/Glbd_01/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Gallbladder/Glbd_01/" 2>/dev/null

# -----------------------------------------------------------------------------
# 6. INTERACT
# -----------------------------------------------------------------------------
echo "Processing: Interact"

# Inte_01
mkdir -p "$DIR_OUT/Interact/Inte_01/"
rsync -avP "$DIR_IN/Interact/Inte_01/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Interact/Inte_01/" 2>/dev/null

# Inte_02
mkdir -p "$DIR_OUT/Interact/Inte_02/"
rsync -avP "$DIR_IN/Interact/Inte_02/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Interact/Inte_02/" 2>/dev/null

# Inte_03
mkdir -p "$DIR_OUT/Interact/Inte_03/"
rsync -avP "$DIR_IN/Interact/Inte_03/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Interact/Inte_03/" 2>/dev/null

# -----------------------------------------------------------------------------
# 7. KIDNEY
# -----------------------------------------------------------------------------
echo "Processing: Kidney"

# Kidn_01
mkdir -p "$DIR_OUT/Kidney/Kidn_01/"
rsync -avP "$DIR_IN/Kidney/Kidn_01/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Kidney/Kidn_01/" 2>/dev/null




head /data/Epic/subprojects/Genetics/sources/Gwas/Kidney/Kidn_01/Data_Received/epic_09dec2008.csv




# -----------------------------------------------------------------------------
# 8. LUNG
# -----------------------------------------------------------------------------
echo "Processing: Lung"

# Lung_01_Icare
mkdir -p "$DIR_OUT/Lung/Lung_01_Icare/"
rsync -avP "$DIR_IN/Lung/Lung_01_Icare/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Lung/Lung_01_Icare/" 2>/dev/null

# -----------------------------------------------------------------------------
# 9. LYMPHOMA
# -----------------------------------------------------------------------------
echo "Processing: Lymphoma"

# Lymp_01_Inter
mkdir -p "$DIR_OUT/Lymphoma/Lymp_01_Inter/"
rsync -avP "$DIR_IN/Lymphoma/Lymp_01_Inter/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Lymphoma/Lymp_01_Inter/" 2>/dev/null

# -----------------------------------------------------------------------------
# 10. NEURO
# -----------------------------------------------------------------------------
echo "Processing: Neuro"

# Neur_01_Park
mkdir -p "$DIR_OUT/Neuro/Neuro_01/"
rsync -avP "$DIR_IN/Neuro/Neuro_01/Data_Received/Genetic data before QC/"*.{map,ped} "$DIR_OUT/Neuro/Neuro_01/"

# -----------------------------------------------------------------------------
# 11. OVARY
# -----------------------------------------------------------------------------
echo "Processing: Ovary"

# Ovar_01_Ocac
mkdir -p "$DIR_OUT/Ovary/Ovar_01_Ocac/"
rsync -avP "$DIR_IN/Ovary/Ovar_01_Ocac/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Ovary/Ovar_01_Ocac/" 2>/dev/null

# -----------------------------------------------------------------------------
# 12. PANCREAS
# -----------------------------------------------------------------------------
echo "Processing: Pancreas"

# Panc_01_Panc
mkdir -p "$DIR_OUT/Pancreas/Panc_01_Panc/"
rsync -avP "$DIR_IN/Pancreas/Panc_01_Panc/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Pancreas/Panc_01_Panc/" 2>/dev/null

# -----------------------------------------------------------------------------
# 13. PROSTATE
# -----------------------------------------------------------------------------
echo "Processing: Prostate"

# Pros_01_Prac
mkdir -p "$DIR_OUT/Prostate/Pros_01_Prac/"
rsync -avP "$DIR_IN/Prostate/Pros_01_Prac/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Prostate/Pros_01_Prac/" 2>/dev/null

# -----------------------------------------------------------------------------
# 14. STOMACH
# -----------------------------------------------------------------------------
echo "Processing: Stomach"

# Stom_01_Sto
mkdir -p "$DIR_OUT/Stomach/Stom_01_Sto/"
rsync -avP "$DIR_IN/Stomach/Stom_01_Sto/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Stomach/Stom_01_Sto/" 2>/dev/null

# -----------------------------------------------------------------------------
# 15. UADT
# -----------------------------------------------------------------------------
echo "Processing: Uadt"

# Uadt_01_Arc
mkdir -p "$DIR_OUT/Uadt/Uadt_01_Arc/"
rsync -avP "$DIR_IN/Uadt/Uadt_01_Arc/Data_Received/"*.{bed,bim,fam} "$DIR_OUT/Uadt/Uadt_01_Arc/" 2>/dev/null

echo "Migration Complete."