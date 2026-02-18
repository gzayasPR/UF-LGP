#!/bin/bash
#==============================================================================
# DETECT INPUT GENOTYPE FORMAT
#==============================================================================
detect_genotype_format() {
    local geno_path=$1
    local plink_flag=""
    
    # STRATEGY 1: Check if input is a complete file with extension
    # Handle VCF files (complete files)
    if [[ -f "${geno_path}" ]]; then
        case "${geno_path}" in
            *.vcf.gz)
                plink_flag="--vcf ${geno_path}"
                echo "INFO: Detected VCF format (gzipped): ${geno_path}" >&2
                ;;
            *.vcf)
                plink_flag="--vcf ${geno_path}"
                echo "INFO: Detected VCF format: ${geno_path}" >&2
                ;;
            *.bcf)
                plink_flag="--bcf ${geno_path}"
                echo "INFO: Detected BCF format: ${geno_path}" >&2
                ;;
            *.bed)
                # For .bed file, strip extension to get prefix
                local prefix="${geno_path%.bed}"
                if [[ -f "${prefix}.bim" && -f "${prefix}.fam" ]]; then
                    plink_flag="--bfile ${prefix}"
                    echo "INFO: Detected binary PLINK format (.bed/.bim/.fam): ${prefix}" >&2
                else
                    echo "ERROR: Found .bed file but missing .bim or .fam: ${prefix}" >&2
                    exit 1
                fi
                ;;
            *.ped)
                # For .ped file, strip extension to get prefix
                local prefix="${geno_path%.ped}"
                if [[ -f "${prefix}.map" ]]; then
                    plink_flag="--file ${prefix}"
                    echo "INFO: Detected PLINK text format (.ped/.map): ${prefix}" >&2
                else
                    echo "ERROR: Found .ped file but missing .map: ${prefix}" >&2
                    exit 1
                fi
                ;;
            *.pgen)
                # For .pgen file, strip extension to get prefix
                local prefix="${geno_path%.pgen}"
                if [[ -f "${prefix}.pvar" && -f "${prefix}.psam" ]]; then
                    plink_flag="--pfile ${prefix}"
                    echo "INFO: Detected PLINK2 format (.pgen/.pvar/.psam): ${prefix}" >&2
                else
                    echo "ERROR: Found .pgen file but missing .pvar or .psam: ${prefix}" >&2
                    exit 1
                fi
                ;;
            *)
                echo "ERROR: Unrecognized file format: ${geno_path}" >&2
                echo "Supported formats: .vcf, .vcf.gz, .bcf, .bed, .ped, .pgen" >&2
                exit 1
                ;;
        esac
        echo "$plink_flag"
        return 0
    fi
    
    # STRATEGY 2: Treat input as a prefix and look for file sets
    local base_path="${geno_path}"
    
    # Check for binary PLINK format (.bed/.bim/.fam)
    if [[ -f "${base_path}.bed" && -f "${base_path}.bim" && -f "${base_path}.fam" ]]; then
        plink_flag="--bfile ${base_path}"
        echo "INFO: Detected binary PLINK format (.bed/.bim/.fam): ${base_path}" >&2
    
    # Check for PLINK 1 text format (.ped/.map)
    elif [[ -f "${base_path}.ped" && -f "${base_path}.map" ]]; then
        plink_flag="--file ${base_path}"
        echo "INFO: Detected PLINK text format (.ped/.map): ${base_path}" >&2
    
    # Check for PLINK 2 binary format (.pgen/.pvar/.psam)
    elif [[ -f "${base_path}.pgen" && -f "${base_path}.pvar" && -f "${base_path}.psam" ]]; then
        plink_flag="--pfile ${base_path}"
        echo "INFO: Detected PLINK2 format (.pgen/.pvar/.psam): ${base_path}" >&2
    
    # Check for VCF with prefix (less common but possible)
    elif [[ -f "${base_path}.vcf.gz" ]]; then
        plink_flag="--vcf ${base_path}.vcf.gz"
        echo "INFO: Detected VCF format (gzipped): ${base_path}.vcf.gz" >&2
    
    elif [[ -f "${base_path}.vcf" ]]; then
        plink_flag="--vcf ${base_path}.vcf"
        echo "INFO: Detected VCF format: ${base_path}.vcf" >&2
    
    # Check for BCF with prefix
    elif [[ -f "${base_path}.bcf" ]]; then
        plink_flag="--bcf ${base_path}.bcf"
        echo "INFO: Detected BCF format: ${base_path}.bcf" >&2
    
    else
        echo "ERROR: Could not detect valid genotype format for: ${geno_path}" >&2
        echo "Searched for:" >&2
        echo "  - Complete file: ${geno_path}" >&2
        echo "  - Prefix-based formats:" >&2
        echo "    * ${base_path}.bed/.bim/.fam (binary PLINK)" >&2
        echo "    * ${base_path}.ped/.map (text PLINK)" >&2
        echo "    * ${base_path}.pgen/.pvar/.psam (PLINK2)" >&2
        echo "    * ${base_path}.vcf or ${base_path}.vcf.gz (VCF)" >&2
        echo "    * ${base_path}.bcf (BCF)" >&2
        exit 1
    fi
    
    echo "$plink_flag"
}