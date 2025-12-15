check_and_fix_rcpp_build <- function(
    verbose = FALSE,
    rtools_root = NULL,
    rtools_dir = NULL, # backward compat
    
    prefer = c("ucrt64", "mingw64", "x86_64-w64-mingw32", "mingw32", "i686-w64-mingw32"),
    try_pkgbuild = TRUE,
    run_rcpp_test = TRUE,
    linux_make_candidates = c("make", "gmake"),
    linux_compiler_candidates = c("g++", "c++", "clang++")
) {
  prefer <- match.arg(prefer)
  
  say <- function(...) if (isTRUE(verbose)) message(...)
  cat_kv <- function(k, v) if (isTRUE(verbose)) message(sprintf("%-18s %s", paste0(k, ":"), v))
  
  # ----- helpers -----
  which_first <- function(cmds) {
    for (cmd in cmds) {
      p <- Sys.which(cmd)
      if (nzchar(p)) return(p)
    }
    ""
  }
  
  is_windows <- identical(.Platform$OS.type, "windows")
  sep <- if (is_windows) ";" else ":"
  
  find_rtools_dirs <- function() {
    drives <- paste0(LETTERS, ":/")
    drives <- drives[dir.exists(drives)]
    hits <- character(0)
    for (d in drives) {
      sub <- tryCatch(list.dirs(d, recursive = FALSE, full.names = TRUE), error = function(e) character(0))
      if (length(sub)) hits <- c(hits, sub[grepl("^rtools", basename(sub), ignore.case = TRUE)])
      hits <- c(hits,
                file.path(d, "rtools"),
                file.path(d, "Rtools"),
                file.path(d, "Rtools44"),
                file.path(d, "Rtools45"),
                file.path(d, "rtools44"),
                file.path(d, "rtools45"))
    }
    hits <- unique(hits)
    hits[dir.exists(hits)]
  }
  
  pick_candidate <- function(rt) {
    usrbin <- file.path(rt, "usr", "bin")
    make  <- file.path(usrbin, "make.exe")
    
    toolbin_pref <- file.path(rt, prefer, "bin")
    ok_pref <- file.exists(make) && dir.exists(toolbin_pref)
    if (ok_pref) {
      return(list(root = rt, usrbin = usrbin, toolbin = toolbin_pref, toolchain = prefer, ok = TRUE))
    }
    
    toolchains <- c("ucrt64", "mingw64", "x86_64-w64-mingw32", "mingw32", "i686-w64-mingw32")
    toolchains <- toolchains[dir.exists(file.path(rt, toolchains, "bin"))]
    
    if (file.exists(make) && length(toolchains) > 0) {
      tc <- toolchains[[1]]
      return(list(root = rt, usrbin = usrbin, toolbin = file.path(rt, tc, "bin"), toolchain = tc, ok = TRUE))
    }
    
    list(root = rt, usrbin = usrbin, toolbin = NA_character_, toolchain = NA_character_, ok = FALSE)
  }
  
  prepend_path <- function(paths) {
    paths <- paths[!is.na(paths) & nzchar(paths) & dir.exists(paths)]
    old <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(paste(paths, collapse = sep), old, sep = sep))
    invisible(old)
  }
  
  rcpp_smoke_test <- function() {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      return(list(ok = FALSE, error = "Rcpp is not installed."))
    }
    
    code <- '
      #include <Rcpp.h>
      using namespace Rcpp;

      // [[Rcpp::export]]
      int add_ints(int a, int b) {
        return a + b;
      }
    '
    
    ok <- TRUE
    err <- NULL
    out <- NULL
    
    tryCatch({
      Rcpp::sourceCpp(code = code, verbose = isTRUE(verbose))
      out <- add_ints(2L, 3L)
      if (!identical(out, 5L)) stop("Unexpected result from compiled code: ", out)
    }, error = function(e) {
      ok <<- FALSE
      err <<- conditionMessage(e)
    })
    
    list(ok = ok, error = err, value = out)
  }
  
  # ----- begin -----
  result <- list(
    ok = FALSE,
    os = .Platform$OS.type,
    sysname = Sys.info()[["sysname"]],
    r_version = R.version.string,
    r_home = R.home(),
    rstudio = nzchar(Sys.getenv("RSTUDIO")) || nzchar(Sys.getenv("RSTUDIO_SESSION_PID")),
    initial = list(
      make = which_first(linux_make_candidates),
      gcc  = Sys.which(if (is_windows) "gcc" else "gcc"),
      gpp  = which_first(linux_compiler_candidates)
    ),
    selected = NULL,
    pkgbuild = NULL,
    rcpp_test = NULL,
    warnings = character(0)
  )
  
  if (isTRUE(verbose)) {
    cat_kv("OS", paste(result$sysname, "(", result$os, ")", sep = ""))
    cat_kv("R", result$r_version)
    cat_kv("R.home", result$r_home)
    cat_kv("RStudio", if (isTRUE(result$rstudio)) "yes" else "no")
    cat_kv("make (pre)", ifelse(nzchar(result$initial$make), result$initial$make, "<not found>"))
    cat_kv("gcc  (pre)", ifelse(nzchar(result$initial$gcc),  result$initial$gcc,  "<not found>"))
    cat_kv("c++  (pre)", ifelse(nzchar(result$initial$gpp),  result$initial$gpp,  "<not found>"))
  }
  
  if (!is_windows) {
    # Linux / macOS: do not patch PATH aggressively; just validate tool presence and run smoke test.
    # If tools missing, provide actionable hint in warnings.
    if (!nzchar(result$initial$make)) {
      result$warnings <- c(result$warnings, "No make found on PATH. Install build-essential (Debian/Ubuntu) or equivalent.")
      say("No make found on PATH")
    }
    if (!nzchar(result$initial$gpp)) {
      result$warnings <- c(result$warnings, "No C++ compiler found on PATH. Install g++ (Debian/Ubuntu: build-essential) or clang++.")
      say("No C++ compiler found on PATH")
    }
    
    if (isTRUE(try_pkgbuild) && requireNamespace("pkgbuild", quietly = TRUE)) {
      ok <- FALSE
      err <- NULL
      tryCatch({
        ok <- pkgbuild::has_build_tools(debug = isTRUE(verbose))
      }, error = function(e) {
        ok <<- FALSE
        err <<- conditionMessage(e)
      })
      result$pkgbuild <- list(ok = ok, error = err)
      if (isTRUE(verbose)) {
        cat_kv("pkgbuild", if (isTRUE(ok)) "OK" else "FAILED")
        if (!isTRUE(ok) && !is.null(err)) say("pkgbuild error: ", err)
      }
    } else {
      result$pkgbuild <- list(ok = NA, error = if (isTRUE(try_pkgbuild)) "pkgbuild not installed" else "skipped")
      if (isTRUE(try_pkgbuild)) say("pkgbuild not installed; skipped pkgbuild validation")
    }
    
    if (isTRUE(run_rcpp_test)) {
      result$rcpp_test <- rcpp_smoke_test()
      if (isTRUE(verbose)) {
        cat_kv("Rcpp test", if (isTRUE(result$rcpp_test$ok)) "OK" else "FAILED")
        if (!isTRUE(result$rcpp_test$ok)) say("Rcpp error: ", result$rcpp_test$error)
      }
    }
    
    # Determine overall success for Linux/macOS
    if (isTRUE(run_rcpp_test)) {
      result$ok <- isTRUE(result$rcpp_test$ok)
    } else {
      result$ok <- nzchar(result$initial$make) && nzchar(result$initial$gpp)
    }
    
    if (!isTRUE(verbose) && !isTRUE(result$ok)) {
      msg <- if (length(result$warnings)) paste(result$warnings, collapse = " | ") else "Build check failed"
      message(msg)
      if (!is.null(result$rcpp_test$error) && nzchar(result$rcpp_test$error)) message("Rcpp error: ", result$rcpp_test$error)
    }
    
    return(invisible(result))
  }
  
  # Windows: locate and link Rtools if needed
  manual_root <- if (!is.null(rtools_root)) rtools_root else rtools_dir
  
  candidate_dirs <- if (!is.null(manual_root)) {
    if (length(manual_root) == 1L) {
      if (!dir.exists(manual_root)) stop("Provided rtools_root does not exist: ", manual_root)
      manual_root
    } else {
      manual_root[dir.exists(manual_root)]
    }
  } else {
    find_rtools_dirs()
  }
  
  if (length(candidate_dirs) == 0) {
    result$warnings <- c(result$warnings, "No Rtools directories found.")
    say("No Rtools directories found")
    if (isTRUE(run_rcpp_test)) {
      result$rcpp_test <- rcpp_smoke_test()
      if (isTRUE(verbose)) {
        cat_kv("Rcpp test", if (isTRUE(result$rcpp_test$ok)) "OK" else "FAILED")
        if (!isTRUE(result$rcpp_test$ok)) say("Rcpp error: ", result$rcpp_test$error)
      }
    }
    return(invisible(result))
  }
  
  picks <- lapply(candidate_dirs, pick_candidate)
  ok_picks <- Filter(function(x) isTRUE(x$ok), picks)
  
  if (length(ok_picks) == 0) {
    result$warnings <- c(result$warnings, "Rtools directories found, but no usable make.exe + toolchain bin pair detected.")
    say("Rtools directories found, but no usable make.exe + toolchain bin pair detected")
    if (isTRUE(run_rcpp_test)) {
      result$rcpp_test <- rcpp_smoke_test()
      if (isTRUE(verbose)) {
        cat_kv("Rcpp test", if (isTRUE(result$rcpp_test$ok)) "OK" else "FAILED")
        if (!isTRUE(result$rcpp_test$ok)) say("Rcpp error: ", result$rcpp_test$error)
      }
    }
    return(invisible(result))
  }
  
  chosen <- ok_picks[[1]]
  result$selected <- chosen
  
  if (isTRUE(verbose)) {
    say("Linking Rtools into PATH for this session")
    cat_kv("Rtools root", chosen$root)
    cat_kv("Toolchain", chosen$toolchain)
  }
  
  prepend_path(c(chosen$usrbin, chosen$toolbin))
  
  # re-check tools
  result$after <- list(
    make = Sys.which("make"),
    gcc  = Sys.which("gcc"),
    gpp  = Sys.which("g++")
  )
  
  if (isTRUE(verbose)) {
    cat_kv("make (post)", ifelse(nzchar(result$after$make), result$after$make, "<not found>"))
    cat_kv("gcc  (post)", ifelse(nzchar(result$after$gcc),  result$after$gcc,  "<not found>"))
    cat_kv("g++  (post)", ifelse(nzchar(result$after$gpp),  result$after$gpp,  "<not found>"))
  }
  
  if (isTRUE(try_pkgbuild) && requireNamespace("pkgbuild", quietly = TRUE)) {
    ok <- FALSE
    err <- NULL
    tryCatch({
      ok <- pkgbuild::has_build_tools(debug = isTRUE(verbose))
    }, error = function(e) {
      ok <<- FALSE
      err <<- conditionMessage(e)
    })
    result$pkgbuild <- list(ok = ok, error = err)
    cat_kv("pkgbuild", if (isTRUE(ok)) "OK" else "FAILED")
    if (!isTRUE(ok) && !is.null(err)) say("pkgbuild error: ", err)
  } else {
    result$pkgbuild <- list(ok = NA, error = if (isTRUE(try_pkgbuild)) "pkgbuild not installed" else "skipped")
    if (isTRUE(try_pkgbuild)) say("pkgbuild not installed; skipped pkgbuild validation")
  }
  
  if (isTRUE(run_rcpp_test)) {
    result$rcpp_test <- rcpp_smoke_test()
    cat_kv("Rcpp test", if (isTRUE(result$rcpp_test$ok)) "OK" else "FAILED")
    if (!isTRUE(result$rcpp_test$ok)) say("Rcpp error: ", result$rcpp_test$error)
  }
  
  # Determine overall success
  if (isTRUE(run_rcpp_test)) {
    result$ok <- isTRUE(result$rcpp_test$ok)
  } else if (is_windows) {
    result$ok <- !is.null(result$selected) && nzchar(Sys.which("make"))
  } else {
    result$ok <- nzchar(result$initial$make) && nzchar(result$initial$gpp)
  }
  
  if (!isTRUE(verbose) && !isTRUE(result$ok)) {
    msg <- if (length(result$warnings)) paste(result$warnings, collapse = " | ") else "Build check failed"
    message(msg)
    if (!is.null(result$rcpp_test$error) && nzchar(result$rcpp_test$error)) message("Rcpp error: ", result$rcpp_test$error)
  }
  
  invisible(result)
}

# Example usage:
# res <- check_and_fix_rcpp_build(verbose = TRUE, prefer = "ucrt64")
# str(res)
