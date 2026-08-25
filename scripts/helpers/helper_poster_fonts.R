# scripts/helpers/helper_poster_fonts.R
#
# Ensures the poster font (Source Sans Pro) is available to ggplot2/ggsave,
# via showtext + sysfonts - this downloads and embeds the Google Font
# (cached locally by sysfonts after the first run), so poster figures
# render with a consistent font regardless of whether it happens to be
# installed as a system font on the machine running the pipeline.
#
# Report figures are unaffected: only theme() calls that explicitly set
# `family = POSTER_FONT_FAMILY` (the poster branches in
# 04_summarize_wsdb_effort.R and 06_summarize_wsdb_seasonality.R) use it.

POSTER_FONT_FAMILY <- "Source Sans Pro"

# Google Fonts now catalogs this typeface as "Source Sans 3" (Adobe's 2021
# rename/relaunch of Source Sans Pro as a variable font); font_add_google()
# below fetches it under that name but registers it locally as
# POSTER_FONT_FAMILY, so the rest of the pipeline can keep referring to it
# as "Source Sans Pro".
POSTER_GOOGLE_FONT_NAME <- "Source Sans 3"

#' Make POSTER_FONT_FAMILY available for ggplot2/ggsave via showtext
#'
#' Call once, before ggsave()-ing any poster figure. `dpi` must match the
#' `dpi` argument used in that ggsave() call - showtext sizes glyphs for a
#' target resolution, so a mismatch makes text render too small/large
#' relative to the image. Safe to call more than once (e.g. once from each
#' poster-producing script); repeat calls are cheap no-ops.
#'
#' @param dpi Resolution (dots per inch) the poster figure will be saved at.
#' @return Invisibly, TRUE if the font is ready to use, FALSE if it fell
#'   back to the default ggplot2 font (a warning explains why).
use_poster_font <- function(dpi = 300) {

  font_pkgs    <- c("sysfonts", "showtext")
  missing_pkgs <- font_pkgs[!vapply(font_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    message("  Installing poster font packages: ", paste(missing_pkgs, collapse = ", "))
    tryCatch(install.packages(missing_pkgs), error = function(e) NULL)
  }

  if (!requireNamespace("sysfonts", quietly = TRUE) ||
      !requireNamespace("showtext", quietly = TRUE)) {
    warning(
      "'sysfonts'/'showtext' unavailable - poster figures will fall back to ",
      "the default ggplot2 font instead of ", POSTER_FONT_FAMILY, ". ",
      "Install with install.packages(c('sysfonts','showtext'))."
    )
    return(invisible(FALSE))
  }

  if (!POSTER_FONT_FAMILY %in% sysfonts::font_families()) {
    tryCatch(
      sysfonts::font_add_google(POSTER_GOOGLE_FONT_NAME, POSTER_FONT_FAMILY),
      error = function(e) {
        warning(
          "Could not fetch '", POSTER_GOOGLE_FONT_NAME, "' from Google Fonts ",
          "(no internet access?). Poster figures will fall back to the ",
          "default ggplot2 font. Original error: ", conditionMessage(e)
        )
      }
    )
  }

  showtext::showtext_auto()
  showtext::showtext_opts(dpi = dpi)

  invisible(POSTER_FONT_FAMILY %in% sysfonts::font_families())
}
