# Dispatcher for the data loaders. Accepts a local path or a Salmobase URL.
# URL downloads are cached under cache/ so repeated knits do not re-fetch.

read_table_anywhere <- function(path_or_url, ...) {
  if (!is_url(path_or_url)) {
    return(readr::read_tsv(path_or_url, show_col_types = FALSE, ...))
  }
  readr::read_tsv(fetch_to_cache(path_or_url),
                  show_col_types = FALSE, ...)
}

# Download a URL once and return its on-disk path. The file is stored under
# cache/{url-hash}.tsv and re-used on subsequent calls. Local paths pass
# through unchanged. process_tissue_data() uses this to feed data.table::fread
# without losing its speed advantage.
fetch_to_cache <- function(path_or_url) {
  if (!is_url(path_or_url)) return(path_or_url)
  cache_dir <- "cache"
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  ext <- tools::file_ext(path_or_url)
  if (!nzchar(ext)) ext <- "tsv"
  cache_file <- file.path(cache_dir,
                          paste0(url_key(path_or_url), ".", ext))
  if (!file.exists(cache_file)) {
    utils::download.file(path_or_url, cache_file, mode = "wb", quiet = FALSE)
  }
  cache_file
}

is_url  <- function(x) grepl("^https?://", x)
url_key <- function(u) paste0("url_", substr(rlang::hash(u), 1, 16))
