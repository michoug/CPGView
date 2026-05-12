irDetect <- function(genome, seed.size = 1000) {
  tick = 0

  # Matching seeds to genome ------------------------------------------------
  genome_rc <- Biostrings::reverseComplement(genome)
  l <- Biostrings::nchar(genome_rc)
  seed_starts <- seq(1, (l - seed.size + 1))
  other_letter <- Biostrings::letterFrequency(genome, letters = "ATCG") != l

  m <- map_genome(
    genome,
    seed_starts = seed_starts,
    seed.size = seed.size,
    other_letter = other_letter
  )

  if (nrow(m) == 0) {
    ir_table <- data.frame(
      chr = "chr1",
      start = 0,
      end = l,
      center = round(l / 2),
      name = "Genome",
      text = paste("Genome:", l),
      stringsAsFactors = FALSE
    )
    ir_table <- gc_count_ir(genome, ir_table)
    return(list(ir_table = ir_table, indel_table = NULL))
  }

  # If IRA covers genome start point, shift genome backward
  while (m$group[1] == 1) {
    tick <- tick + seed.size
    if (tick >= l) {
      break
    } # safety guard
    genome <- c(genome[(l - tick + 1):l], genome[1:(l - tick)])
    m <- map_genome(
      genome,
      seed_starts = seed_starts,
      seed.size = seed.size,
      other_letter = other_letter
    )
  }

  # If IRB covers genome end point, shift genome forward
  while (m$group[nrow(m)] + seed.size - 1 == l) {
    tick <- tick - seed.size

    if (tick == 0) {
      # No shift needed — genome already correctly oriented
      m <- map_genome(
        genome,
        seed_starts = seed_starts,
        seed.size = seed.size,
        other_letter = other_letter
      )
      break
    } else if (tick < 0) {
      shift <- -tick # positive number of bases to move to front
      if (shift >= l) {
        break
      } # safety guard
      genome <- c(genome[(l - shift + 1):l], genome[1:(l - shift)])
    } else {
      # tick still positive
      if ((-tick + 1) < 1 || (-tick) < 1) {
        break
      } # safety guard
      genome <- c(genome[(-tick + 1):l], genome[1:-tick])
    }

    m <- map_genome(
      genome,
      seed_starts = seed_starts,
      seed.size = seed.size,
      other_letter = other_letter
    )
  }

  m <- m[!duplicated(m$group), ]

  df <- data.frame(
    group_before = m$group[-nrow(m)],
    group_after = m$group[-1],
    group_diff = m$group[-1] - m$group[-nrow(m)] - 1,
    start_diff = m$start[-1] - m$start[-nrow(m)] - 1,
    start_before = m$start[-nrow(m)],
    start_after = m$start[-1],
    end_diff = m$end[-1] - m$end[-nrow(m)] - 1,
    end_before = m$end[-nrow(m)],
    end_after = m$end[-1],
    stringsAsFactors = FALSE
  )

  pos <- df[c(1, nrow(df)), ] %>%
    rbind.data.frame(dplyr::filter(df, group_diff != 0)) %>%
    dplyr::arrange(group_before)

  if (sum(pos$start_diff < 0) > 0) {
    tmp <- pos[pos$start_diff < 0, ]
    pos <- rbind.data.frame(
      pos,
      df[which(df$group_diff %in% tmp$group_diff) - 1, ],
      df[which(df$group_diff %in% tmp$group_diff) + 1, ]
    ) %>%
      dplyr::filter(start_diff >= 0) %>%
      dplyr::arrange(group_before)

    mismatch_group <- which(pos$group_diff > 0)
    ira_s <- pos$group_before[mismatch_group[1] - 1]
    ira_e <- pos$group_before[which.max(pos$group_diff)] + seed.size - 1
    irb_s <- l - pos$start_before[which.max(pos$group_diff)] - seed.size + 2
    irb_e <- pos$group_after[mismatch_group[length(mismatch_group)] + 1] +
      seed.size -
      1
  } else {
    ira_s <- pos$group_before[1]
    ira_e <- pos$group_before[(nrow(pos) + 1) / 2] + seed.size - 1
    irb_s <- pos$group_after[(nrow(pos) + 1) / 2]
    irb_e <- pos$group_after[nrow(pos)] + seed.size - 1
  }

  ira_len <- ira_e - ira_s + 1
  irb_len <- irb_e - irb_s + 1
  lsc_len <- ira_s - 1 + l - irb_e
  ssc_len <- irb_s - ira_e - 1

  # Detecting indels and replacements in IR ----------------------------------
  if (nrow(pos) > 3) {
    ira_seq <- genome[ira_s:ira_e]
    irb_seq <- genome[irb_s:irb_e]
    other_letter <- (Biostrings::letterFrequency(ira_seq, letters = "ATCG") !=
      Biostrings::nchar(ira_seq)) |
      (Biostrings::letterFrequency(irb_seq, letters = "ATCG") !=
        Biostrings::nchar(irb_seq))
    indel_table <- detect_mismatch(
      ira_seq = ira_seq,
      irb_seq = irb_seq,
      ira_s = ira_s,
      irb_s = irb_s,
      other_letter = other_letter
    )
  } else {
    indel_table <- NULL
  }

  if (ira_len != irb_len) {
    if (
      is.null(indel_table) |
        sum(indel_table$mismatch_type %in% c("insert", "delete")) == 0 |
        abs(ira_len - irb_len) > 100
    ) {
      stop("The IR regions are not in similar length and no indel was detected")
    }
  }

  # Recover original coordinates (0-base) ------------------------------------
  if (tick == 0) {
    ir_table <- data.frame(
      chr = rep("chr1", 5),
      start = c(0, ira_s - 1, ira_e, irb_s - 1, irb_e),
      end = c(ira_s - 1, ira_e, irb_s - 1, irb_e, l),
      name = c("LSC", "IRA", "SSC", "IRB", "LSC"),
      text = c(
        paste("LSC:", lsc_len),
        paste("IRA:", ira_len),
        paste("SSC:", ssc_len),
        paste("IRB:", irb_len),
        ""
      ),
      stringsAsFactors = FALSE
    )
  } else if (tick > 0) {
    if ((tick - ira_s) == 0) {
      ir_table <- data.frame(
        chr = rep("chr1", 4),
        start = c(0, ira_e - tick, irb_s - tick, irb_e - tick),
        end = c(ira_e - tick, irb_s - tick, irb_e - tick, l),
        name = c("IRA", "SSC", "IRB", "LSC"),
        text = c(
          paste("IRA:", ira_len),
          paste("SSC:", ssc_len),
          paste("IRB:", irb_len),
          paste("LSC:", lsc_len)
        ),
        stringsAsFactors = FALSE
      )
    } else {
      ir_table <- data.frame(
        chr = rep("chr1", 5),
        start = c(
          0,
          ira_e - tick,
          irb_s - tick - 1,
          irb_e - tick,
          l - (tick - ira_s) - 1
        ),
        end = c(
          ira_e - tick,
          irb_s - tick - 1,
          irb_e - tick,
          l - (tick - ira_s) - 1,
          l
        ),
        name = c("IRA", "SSC", "IRB", "LSC", "IRA"),
        text = c(
          paste("IRA:", ira_len),
          paste("SSC:", ssc_len),
          paste("IRB:", irb_len),
          paste("LSC:", lsc_len),
          ""
        ),
        stringsAsFactors = FALSE
      )
    }
  } else if (tick < 0) {
    if ((irb_e - tick - l) == 0) {
      ir_table <- data.frame(
        chr = rep("chr1", 4),
        start = c(0, ira_s - tick - 1, ira_e - tick, irb_s - tick - 1),
        end = c(ira_s - tick - 1, ira_e - tick, irb_s - tick - 1, l),
        name = c("LSC", "IRA", "SSC", "IRB"),
        text = c(
          paste("LSC:", lsc_len),
          paste("IRA:", irb_len),
          paste("SSC:", ssc_len),
          paste("IRB:", irb_len)
        ),
        stringsAsFactors = FALSE
      )
    } else {
      ir_table <- data.frame(
        chr = rep("chr1", 5),
        start = c(
          0,
          (irb_e - tick - l),
          ira_s - tick - 1,
          ira_e - tick,
          irb_s - tick - 1
        ),
        end = c(
          (irb_e - tick - l),
          ira_s - tick - 1,
          ira_e - tick,
          irb_s - tick - 1,
          l
        ),
        name = c("IRB", "LSC", "IRA", "SSC", "IRB"),
        text = c(
          paste("IRB:", irb_len),
          paste("LSC:", lsc_len),
          paste("IRA:", irb_len),
          paste("SSC:", ssc_len),
          ""
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  # Add text coordinates
  ir_table$center <- round((ir_table$start + ir_table$end) / 2, 0)
  if (nrow(ir_table) == 5) {
    ir_table$center[1] <- ir_table$center[1] +
      round(l / 2) +
      round(ir_table$start[5] / 2)
  }
  if (ir_table$center[1] > l) {
    ir_table$center[1] <- ir_table$center[1] - l
  }

  # Add GC count
  ir_table <- gc_count_ir(genome, ir_table)

  if (lsc_len < ssc_len) {
    ir_table$name <- sub("LSC", "tmp", ir_table$name, fixed = TRUE)
    ir_table$name <- sub("SSC", "LSC", ir_table$name, fixed = TRUE)
    ir_table$name <- sub("tmp", "SSC", ir_table$name, fixed = TRUE)
    ir_table$text <- sub("LSC", "tmp", ir_table$text, fixed = TRUE)
    ir_table$text <- sub("SSC", "LSC", ir_table$text, fixed = TRUE)
    ir_table$text <- sub("tmp", "SSC", ir_table$text, fixed = TRUE)
  }

  return(list(ir_table = ir_table, indel_table = indel_table))
}
