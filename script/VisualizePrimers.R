setwd(paste0(Sys.getenv("PROJECTS_ROOT"), "/Bioinformatics/workspace/other/primers_visualization_tool"))

# Load necessary libraries
library(tidyverse)
library(ggtranscript)    # For plotting transcripts
library(rtracklayer)
library(GenomicFeatures) # For TxDb and related functions
library(Rsamtools)       # For FaFile
library(Biostrings)      # For DNAString, vmatchPattern, etc.
library(GenomicRanges)   # For GRanges operations
library(readxl)

# Set working directory (optional, adjust as needed)
# setwd("/Users/tucos/Downloads/primer plot")

# --- Configuration ---
# Paths to your files:
gtf_file   <- "/Volumes/TucciSSD/Bioinformatics/resources/genomes/gencode/hsapiens/GRCh38/gencode.v38.primary_assembly.annotation.gtf"
fasta_file <- "/Volumes/TucciSSD/Bioinformatics/resources/genomes/gencode/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa"

# TxDb file (will be created if it doesn't exist to speed up subsequent runs)
txdb_path <- "~/.bioinfoCache/GRCh38_gencode.v38_txdb.sqlite" # Changed variable name for clarity
txdb_dir <- dirname(txdb_path)
if (!dir.exists(txdb_dir)) {
  message("Creating directory for TxDb cache: ", txdb_dir)
  dir.create(txdb_dir, recursive = TRUE, showWarnings = FALSE)
}


# --- Load Annotations and Fasta sequence ---
message("Loading GTF annotations (this might take a moment)...")
gtf_data <- import(gtf_file) # Loaded once globally

message("Preparing transcript database (TxDb)...")
# Build or load TxDb
if (file.exists(txdb_path)) {
  txdb <- loadDb(txdb_path)
  message("Loaded existing TxDb from ", txdb_path)
} else {
  message("Building TxDb from GTF (this may take a while for the first run)...")
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  saveDb(txdb, file = txdb_path)
  message("Saved TxDb to ", txdb_path)
}

# Load FASTA file
fasta <- FaFile(fasta_file) # Loaded once globally

# Define main function ---------
mapAndPlotPrimers <- function(gene_symbol_param, fw_primer, rev_primer, output_file_label = "", max_mismatch = 0) {
  
  # --- 1. Filter the Annotations for the Gene of Interest (GoI) ---
  message(paste0("Filtering annotations for gene: ", gene_symbol_param, "..."))
  GoI_annotation_full <- gtf_data %>%
    as_tibble() %>%
    filter(
      !is.na(gene_name),
      gene_name == gene_symbol_param, # Use function parameter
      type %in% c("exon", "transcript"), 
      tag %in% c("CCDS", "basic"),
      !is.na(havana_transcript),
      transcript_type == "protein_coding"
    )
  
  GoI_annotation_full <- GoI_annotation_full %>% filter(!is.na(transcript_id))
  selected_transcript_ids <- GoI_annotation_full %>%
    filter(type == "transcript") %>% 
    pull(transcript_id) %>%
    unique()
  
  if (length(selected_transcript_ids) == 0) {
    warning(paste("No transcripts found for gene", gene_symbol_param, "with the specified criteria. Plot will be empty."))
    return(NULL) # Exit function if no transcripts
  }
  
  GoI_annotation <- GoI_annotation_full %>%
    filter(transcript_id %in% selected_transcript_ids)
  
  GoI_exons_original_df <- GoI_annotation %>%
    filter(type == "exon") %>%
    mutate(
      seqnames = as.character(seqnames),
      strand = as.character(strand),
      start = as.integer(start),
      end = as.integer(end),
      transcript_id = as.character(transcript_id)
    ) %>%
    mutate(original_exon_identifier = paste(seqnames, start, end, strand, sep="_"))
  
  if (nrow(GoI_exons_original_df) == 0) {
    warning(paste("No exons found for the selected transcripts of gene", gene_symbol_param, ". Plot will be empty."))
    return(NULL) # Exit function if no exons
  }
  
  # --- 2. Prepare Transcript Sequences and Match Primers ---
  # convert primer sequences into DNA strings
  fw_primer <- DNAString(fw_primer)
  rev_primer <- DNAString(rev_primer)
  
  message("Extracting transcript sequences and matching primers...")
  exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
  exons_by_tx_selected <- exons_by_tx[names(exons_by_tx) %in% selected_transcript_ids]
  
  if (length(exons_by_tx_selected) == 0) {
    stop("No exons found in TxDb for the selected transcript IDs. Ensure transcript IDs match between GTF parsing and TxDb.")
  }
  
  tx_seqs <- extractTranscriptSeqs(fasta, exons_by_tx_selected)
  
  fwd_hits_on_tx <- vmatchPattern(fw_primer, tx_seqs, max.mismatch=max_mismatch, with.indels=FALSE) 
  rev_hits_on_tx <- vmatchPattern(reverseComplement(rev_primer), tx_seqs, max.mismatch=max_mismatch, with.indels=FALSE)
  
  # --- 2b. Calculate Amplicon Lengths ---
  message("Calculating amplicon lengths...")
  amplicon_data <- list()
  for (tx_id in names(tx_seqs)) {
    fwd_hits_current_tx <- fwd_hits_on_tx[[tx_id]]
    rev_hits_current_tx <- rev_hits_on_tx[[tx_id]]
    
    if (length(fwd_hits_current_tx) > 0 && length(rev_hits_current_tx) > 0) {
      min_fwd_start_on_tx <- min(start(fwd_hits_current_tx))
      max_rev_end_on_tx <- max(end(rev_hits_current_tx))
      
      if (min_fwd_start_on_tx < max_rev_end_on_tx) { # Ensure FWD is upstream of REV complement
        amplicon_length <- max_rev_end_on_tx - min_fwd_start_on_tx + 1
        amplicon_data[[tx_id]] <- list(transcript_id = tx_id, amplicon_length = amplicon_length)
      }
    }
  }
  amplicon_lengths_df <- bind_rows(amplicon_data)
  
  if (nrow(amplicon_lengths_df) > 0) {
    message("Calculated amplicon lengths (on spliced transcript):")
    # Prepare a nice table for printing
    amplicon_summary_print <- amplicon_lengths_df %>%
      mutate(transcript_id_short = str_remove(transcript_id, "\\..*$")) %>%
      dplyr::select(Transcript = transcript_id_short, `Amplicon Length (bp)` = amplicon_length) %>%
      arrange(Transcript)
    print(as.data.frame(amplicon_summary_print)) # Print as data.frame for cleaner console output
  } else {
    message("No valid amplicons formed by the provided primers for the selected transcripts.")
    # skip the execution of the rest of the function
    return(NULL)
    
  }
  
  # --- 3. Map Transcript-Relative Primer Hits to Genomic Coordinates ---
  message("Mapping primer hits to genomic coordinates...")
  all_primer_genomic_segments <- GRanges()
  
  for (tx_id in names(exons_by_tx_selected)) {
    exons_for_this_tx <- exons_by_tx_selected[[tx_id]] 
    tx_strand <- as.character(strand(exons_for_this_tx)[1]) 
    
    current_tx_seq_pos <- 1
    exon_map_list <- list()
    for (i in seq_along(exons_for_this_tx)) {
      exon_gr <- exons_for_this_tx[i]
      exon_len <- width(exon_gr)
      exon_map_list[[i]] <- list(
        exon_gr = exon_gr,
        tx_seq_start = current_tx_seq_pos,
        tx_seq_end = current_tx_seq_pos + exon_len - 1
      )
      current_tx_seq_pos <- current_tx_seq_pos + exon_len
    }
    
    map_hit_to_genomic <- function(tx_hit_range, primer_type_label) {
      hit_genomic_segments <- GRanges()
      if (length(tx_hit_range) == 0) return(hit_genomic_segments)
      
      for (j in 1:length(tx_hit_range)) { 
        hit_on_tx <- tx_hit_range[j]
        for (map_entry in exon_map_list) {
          exon_gr_from_map <- map_entry$exon_gr
          overlap_tx_start <- max(start(hit_on_tx), map_entry$tx_seq_start)
          overlap_tx_end <- min(end(hit_on_tx), map_entry$tx_seq_end)
          
          if (overlap_tx_start <= overlap_tx_end) { 
            hit_start_in_exon_seq_part <- overlap_tx_start - map_entry$tx_seq_start + 1
            hit_end_in_exon_seq_part <- overlap_tx_end - map_entry$tx_seq_start + 1
            
            genomic_seg_start <- NA
            genomic_seg_end <- NA
            if (tx_strand == "+") {
              genomic_seg_start <- start(exon_gr_from_map) + hit_start_in_exon_seq_part - 1
              genomic_seg_end <- start(exon_gr_from_map) + hit_end_in_exon_seq_part - 1
            } else if (tx_strand == "-") {
              genomic_seg_start <- end(exon_gr_from_map) - hit_end_in_exon_seq_part + 1
              genomic_seg_end <- end(exon_gr_from_map) - hit_start_in_exon_seq_part + 1
            }
            
            current_gen_segment <- GRanges(
              seqnames = seqnames(exon_gr_from_map),
              ranges = IRanges(genomic_seg_start, genomic_seg_end),
              strand = strand(exon_gr_from_map),
              transcript_id = tx_id,
              primer_type = primer_type_label,
              original_exon_identifier = paste(as.character(seqnames(exon_gr_from_map)), start(exon_gr_from_map), end(exon_gr_from_map), as.character(strand(exon_gr_from_map)), sep="_")
            )
            hit_genomic_segments <- c(hit_genomic_segments, current_gen_segment)
          }
        }
      }
      return(hit_genomic_segments)
    }
    
    fwd_gen_segs <- map_hit_to_genomic(fwd_hits_on_tx[[tx_id]], "Forward Primer")
    if (length(fwd_gen_segs) > 0) all_primer_genomic_segments <- c(all_primer_genomic_segments, fwd_gen_segs)
    
    rev_gen_segs <- map_hit_to_genomic(rev_hits_on_tx[[tx_id]], "Reverse Primer")
    if (length(rev_gen_segs) > 0) all_primer_genomic_segments <- c(all_primer_genomic_segments, rev_gen_segs)
  }
  
  # --- 4. Create "Virtual" Exons by Splitting Original Exons ---
  message("Creating virtual exons based on primer matches...")
  new_virtual_exon_rows <- list()
  
  for (i in 1:nrow(GoI_exons_original_df)) {
    current_original_exon_row <- GoI_exons_original_df[i, ]
    current_original_exon_gr <- GRanges(
      seqnames = current_original_exon_row$seqnames,
      ranges = IRanges(current_original_exon_row$start, current_original_exon_row$end),
      strand = current_original_exon_row$strand
    )
    
    primer_segments_on_this_exon <- GRanges()
    if (length(all_primer_genomic_segments) > 0) {
      primer_segments_on_this_exon <- subset(all_primer_genomic_segments, 
                                             original_exon_identifier == current_original_exon_row$original_exon_identifier &
                                               transcript_id == current_original_exon_row$transcript_id)
    }
    
    if (length(primer_segments_on_this_exon) == 0) {
      new_row <- current_original_exon_row
      new_row$feature_type <- "Other Exons"
      new_virtual_exon_rows <- append(new_virtual_exon_rows, list(new_row))
    } else {
      for (k in 1:length(primer_segments_on_this_exon)) {
        primer_seg_gr <- primer_segments_on_this_exon[k]
        new_row <- current_original_exon_row 
        new_row$start <- start(primer_seg_gr)
        new_row$end <- end(primer_seg_gr)
        new_row$feature_type <- primer_seg_gr$primer_type 
        new_virtual_exon_rows <- append(new_virtual_exon_rows, list(new_row))
      }
      
      if (length(primer_segments_on_this_exon) > 0) {
        covered_by_primers_on_exon <- GenomicRanges::reduce(primer_segments_on_this_exon)
        non_primer_parts_gr <- GenomicRanges::setdiff(current_original_exon_gr, covered_by_primers_on_exon)
        
        for (k in 1:length(non_primer_parts_gr)) {
          non_primer_part_gr_single <- non_primer_parts_gr[k]
          if (width(non_primer_part_gr_single) > 0) {
            new_row <- current_original_exon_row
            new_row$start <- start(non_primer_part_gr_single)
            new_row$end <- end(non_primer_part_gr_single)
            new_row$feature_type <- "Target Exon(s)"
            new_virtual_exon_rows <- append(new_virtual_exon_rows, list(new_row))
          }
        }
      } else { 
        new_row <- current_original_exon_row
        new_row$feature_type <- "Other Exons" 
        new_virtual_exon_rows <- append(new_virtual_exon_rows, list(new_row))
      }
    }
  }
  
  GoI_exons_split_df <- dplyr::bind_rows(new_virtual_exon_rows)
  
  if (nrow(GoI_exons_split_df) == 0) { # Should not happen if GoI_exons_original_df was not empty
    message("No exon segments (original or virtual) to plot. This might indicate an issue in the splitting logic or that all exons were consumed by primer matches without non-matching parts, which is unlikely.")
    # Add the original exons back if split_df is empty but original had content
    if(nrow(GoI_exons_original_df) > 0 && !"feature_type" %in% colnames(GoI_exons_original_df)){
      GoI_exons_split_df <- GoI_exons_original_df %>% mutate(feature_type = "Other Exons")
      message("Using original exons as fallback for plotting.")
    } else if (nrow(GoI_exons_original_df) > 0) {
      GoI_exons_split_df <- GoI_exons_original_df # if feature_type was already there
      message("Using original exons as fallback for plotting.")
    } else {
      warning("No exons available for plotting after splitting.")
      return(NULL)
    }
  }
  
  
  if ("exon_id" %in% colnames(GoI_exons_split_df)) {
    GoI_exons_split_df$exon_id_derived <- paste(GoI_exons_split_df$exon_id, GoI_exons_split_df$feature_type, GoI_exons_split_df$start, GoI_exons_split_df$end, sep = "_")
  } else {
    GoI_exons_split_df$exon_id_derived <- paste(GoI_exons_split_df$transcript_id, GoI_exons_split_df$feature_type, GoI_exons_split_df$start, GoI_exons_split_df$end, sep = "_")
  }
  GoI_exons_split_df$exon_id_derived <- make.unique(GoI_exons_split_df$exon_id_derived)
  GoI_exons_split_df$type <- "exon"
  
  # --- 5. Rescale Coordinates and Plot ---
  message("Rescaling coordinates for plotting...")
  
  GoI_introns_for_rescaling <- to_intron(GoI_exons_split_df, group_var = "transcript_id")
  
  GoI_rescaled <- shorten_gaps(
    exons = GoI_exons_split_df,
    introns = GoI_introns_for_rescaling,
    group_var = "transcript_id"
  )
  
  GoI_rescaled_exons <- GoI_rescaled %>% dplyr::filter(type == "exon")
  GoI_rescaled_introns <- GoI_rescaled %>% dplyr::filter(type == "intron")
  
  # Prepare y-axis labels with amplicon length
  y_labels_map_df <- GoI_rescaled_exons %>%
    distinct(transcript_id) %>%
    mutate(transcript_id_short = str_remove(transcript_id, "\\..*$")) %>%
    left_join(amplicon_lengths_df, by = "transcript_id") %>%
    mutate(
      y_display_label = ifelse(
        !is.na(amplicon_length),
        paste0(transcript_id_short, " (Amplicon: ", amplicon_length, " bp)"),
        transcript_id_short
      )
    ) %>%
    dplyr::select(transcript_id, y_display_label)
  
  # Join display labels to plotting data
  GoI_rescaled_exons_for_plot <- GoI_rescaled_exons %>%
    left_join(y_labels_map_df, by = "transcript_id")
  
  GoI_rescaled_introns_for_plot <- GoI_rescaled_introns %>%
    left_join(y_labels_map_df, by = "transcript_id")
  
  
  message("Generating plot...")
  feature_colors <- c(
    "Forward Primer" = "orangered",
    "Reverse Primer" = "dodgerblue",
    "Target Exon(s)" = "lightgoldenrodyellow",
    "Other Exons" = "lightgrey" 
  )
  # prepone "_" to the file name labelf it it was specified
  if(output_file_label != "") {
    output_file_label <- paste0("_", output_file_label)
  }
  
  # extract number of transcripts to plot
  n_of_transcripts = length(unique(GoI_rescaled_exons_for_plot$transcript_id))
  plot_length = abs(min(GoI_rescaled_exons_for_plot$start,GoI_rescaled_exons_for_plot$end) - max(GoI_rescaled_exons_for_plot$start,GoI_rescaled_exons_for_plot$end))
  # plot
  final_plot <- GoI_rescaled_exons_for_plot %>%
    ggplot(aes(
      xstart = start, 
      xend = end,     
      y = y_display_label, # Use the new label with amplicon info
      fill = feature_type 
    )) +
    geom_range(
      colour = "black", 
      linewidth = 0.3 
    ) +
    geom_intron(
      data = GoI_rescaled_introns_for_plot, # Use augmented intron data
      aes(strand = strand), # y is inherited
      arrow.min.intron.length = 200
    ) +
    scale_fill_manual(
      values = feature_colors,
      name = "", # Legend title removed for brevity, or customize
      breaks = names(feature_colors)[names(feature_colors) %in% unique(GoI_rescaled_exons_for_plot$feature_type)]
    ) +
    labs(
      title = paste("Gene:", gene_symbol_param, str_replace(output_file_label, "_", " - ")), # Use function parameter
      x = "Rescaled Genomic Position (bp)",
      y = "Transcript ID" # y-axis title
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10), # Adjust size if y_display_label is long
      legend.position = "bottom",
      legend.title = element_text(face="bold")
    ) +
    # Define plot size based on transcript length and number of different transcripts
    ggh4x::force_panelsizes(rows = unit(n_of_transcripts*1.5, "cm"), cols = unit(plot_length/300, "cm"))

  # export the plot
  ggsave(filename = paste0("results/primers_", gene_symbol_param, output_file_label, ".pdf"), 
         plot = final_plot, 
         device = "pdf", 
         height = (n_of_transcripts*1.5) + 5, 
         width = (plot_length/300) + 9,
         units = "cm")
}

# --- Example Run ---
# Your target gene symbol:
gene_symbol <- "KLK3"

# Define primers as DNAString objects
fw_primer <- "CGCAAGTTCACCCTCAGAAGGT"
rev_primer <- "GACGTGATACCTTGAAGCACACC"

# Call the main function
plot_object <- mapAndPlotPrimers(gene = gene_symbol, 
                                 fw = fw_primer, 
                                 rev = rev_primer, max_mismatch = 0, output_file_label = "Origene HP227909")

# Run over list of primers in xlsx fle --------
primers_list <- read_excel("data/Primers.xlsx") %>%
  mutate(gene_name = str_remove(Name, "^(h|m)") %>% str_extract("^[^_]+") %>% recode("CGA" = "CHGA"),
         primer = ifelse(str_detect(Name, "Fw"), "Forward", "Reverse")) %>%
  dplyr::select(-Name, -`Product Length`) %>%
  spread(key = "primer", value = "Sequence")

for(i in 1:nrow(primers_list)) {
  mapAndPlotPrimers(gene_symbol_param = primers_list$gene_name[i], 
                    fw_primer = primers_list$Forward[i], 
                    rev_primer = primers_list$Reverse[i],
                    output_file_label =primers_list$Source[i]
                    )
}

