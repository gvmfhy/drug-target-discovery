    # Load required packages
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    
    required_packages <- c("AnnotationDbi", "hgu133plus2.db", "GEOquery")
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE))
            BiocManager::install(pkg, ask = FALSE)
    }
    
    library(AnnotationDbi)
    library(hgu133plus2.db)
    library(GEOquery)
    
    # Read probe IDs
    args <- commandArgs(trailingOnly = TRUE)
    probe_file <- args[1]
    output_file <- args[2]
    
    probes <- readLines(probe_file)
    cat("Read", length(probes), "probe IDs\\n")
    
    # Try direct mapping first
    mapping_df <- NULL
    
    # Method 1: Using select with direct package
    tryCatch({
        cat("Trying direct select method...\\n")
        valid_keys <- keys(hgu133plus2.db, keytype="PROBEID")
        valid_probes <- probes[probes %in% valid_keys]
        cat("Valid probes:", length(valid_probes), "out of", length(probes), "\\n")
        
        if (length(valid_probes) > 0) {
            result <- select(hgu133plus2.db,
                             keys = valid_probes,
                             columns = "SYMBOL",
                             keytype = "PROBEID")
            if (nrow(result) > 0) {
                mapping_df <- result
                cat("Direct mapping found", nrow(result), "mappings\\n")
            }
        }
    }, error = function(e) {
        cat("Error in direct mapping:", conditionMessage(e), "\\n")
    })
    
    # Method 2: Using GEOquery platform annotation
    if (is.null(mapping_df) || nrow(mapping_df) == 0) {
        tryCatch({
            cat("Trying GEOquery platform annotation...\\n")
            gpl <- getGEO("GPL570")
            anno_table <- Table(gpl)
            
            # Find ID and Symbol columns
            id_col <- NULL
            for (col_name in c("ID", "SPOT_ID", "PROBEID", "Probe Set ID")) {
                if (col_name %in% colnames(anno_table)) {
                    id_col <- col_name
                    break
                }
            }
            
            symbol_col <- NULL
            for (col_name in c("Gene Symbol", "SYMBOL", "gene_symbol", "Gene_Symbol")) {
                if (col_name %in% colnames(anno_table)) {
                    symbol_col <- col_name
                    break
                }
            }
            
            if (!is.null(id_col) && !is.null(symbol_col)) {
                cat("Using columns:", id_col, "->", symbol_col, "\\n")
                
                # Create mapping dataframe
                mapping_df <- data.frame(
                    PROBEID = anno_table[[id_col]],
                    SYMBOL = anno_table[[symbol_col]],
                    stringsAsFactors = FALSE
                )
                
                # Filter to only include requested probes
                mapping_df <- mapping_df[mapping_df$PROBEID %in% probes,]
                cat("Platform annotation found", nrow(mapping_df), "mappings\\n")
            } else {
                cat("Could not find appropriate ID/Symbol columns in platform annotation\\n")
            }
        }, error = function(e) {
            cat("Error in GEOquery mapping:", conditionMessage(e), "\\n")
        })
    }
    
    # Method 3: Using environment mapping
    if (is.null(mapping_df) || nrow(mapping_df) == 0) {
        tryCatch({
            cat("Trying environment mapping...\\n")
            symbol_list <- as.list(hgu133plus2SYMBOL)
            probe_ids <- names(symbol_list)
            symbols <- as.character(symbol_list)
            
            mapping_df <- data.frame(
                PROBEID = probe_ids,
                SYMBOL = symbols,
                stringsAsFactors = FALSE
            )
            
            # Filter to only include requested probes
            mapping_df <- mapping_df[mapping_df$PROBEID %in% probes,]
            cat("Environment mapping found", nrow(mapping_df), "mappings\\n")
        }, error = function(e) {
            cat("Error in environment mapping:", conditionMessage(e), "\\n")
        })
    }
    
    # Write output file
    if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
        # Remove rows with empty symbols
        mapping_df <- mapping_df[mapping_df$SYMBOL != "" & !is.na(mapping_df$SYMBOL),]
        write.csv(mapping_df, output_file, row.names = FALSE)
        cat("Wrote", nrow(mapping_df), "mappings to output file\\n")
    } else {
        # Write empty file
        write.csv(data.frame(PROBEID=character(), SYMBOL=character()), output_file, row.names = FALSE)
        cat("No mappings found, wrote empty output file\\n")
    }