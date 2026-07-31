library(xml2)

p <- "E:/20241004_MCI/20260709/Manuscript_FIO_CNS_story_logic_revised_20260720_104738.html"
d <- read_html(p)
imgs <- xml_attr(xml_find_all(d, "//img"), "src")
missing <- imgs[!file.exists(file.path(dirname(p), imgs))]
captions <- trimws(xml_text(xml_find_all(d, "//figcaption")))
pmid_tags <- xml_find_all(d, "//*[contains(., 'PMID:')]")
old_markers <- c("Figure 8.", "Figure5_CSF_immune_landscape", "Figure6_CellChat_group_level", "figures_story_20260716_213246/Figure", "direct GSE63060-TSMR overlap markers, <em>RPS27L</em>")
html <- paste(readLines(p, warn = FALSE), collapse = "\n")
present_old_markers <- old_markers[vapply(old_markers, function(x) grepl(x, html, fixed = TRUE), logical(1))]

cat("HTML_PARSE=PASS\n")
cat("FIGURE_COUNT=", length(imgs), "\n", sep = "")
cat("MISSING_IMAGE_COUNT=", length(missing), "\n", sep = "")
if (length(missing)) cat(paste(missing, collapse = "\n"), "\n")
cat("PMID_TAG_COUNT=", length(pmid_tags), "\n", sep = "")
cat("OLD_MARKER_COUNT=", length(present_old_markers), "\n", sep = "")
if (length(present_old_markers)) cat(paste(present_old_markers, collapse = "\n"), "\n")
cat("CAPTIONS\n", paste(captions, collapse = "\n"), "\n", sep = "")
cat("FIG4_RESULT\n", trimws(xml_text(xml_find_first(d, "//section[h2='Results']//h3[contains(., 'donor-aware')]/following-sibling::p[1]"))), "\n", sep = "")

stopifnot(length(missing) == 0L, length(present_old_markers) == 0L, length(imgs) == 7L)
