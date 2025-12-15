png_gallery_to_pdf_paginated <- function(
    folder,
    output_pdf = "gallery_A4.pdf",
    images_per_row = 2,
    rows_per_page = 3,
    title_fontsize = 10
) {
  if (!requireNamespace("png", quietly = TRUE)) stop("Please install the 'png' package.")
  
  files <- list.files(folder, pattern = "\\.png$", full.names = TRUE)
  if (length(files) == 0) stop("No PNG files found in the folder.")
  
  total_per_page <- images_per_row * rows_per_page
  num_pages <- ceiling(length(files) / total_per_page)
  
  # A4 portrait in inches: 8.27 x 11.69
  grDevices::pdf(output_pdf, width = 8.27, height = 11.69)
  
  for (page in seq_len(num_pages)) {
    grid::grid.newpage()
    start <- (page - 1) * total_per_page + 1
    end <- min(page * total_per_page, length(files))
    page_files <- files[start:end]
    
    layout <- grid::grid.layout(rows_per_page, images_per_row)
    grid::pushViewport(grid::viewport(layout = layout))
    
    for (i in seq_along(page_files)) {
      row <- ceiling(i / images_per_row)
      col <- i - (row - 1) * images_per_row
      
      img <- png::readPNG(page_files[i])
      grob_img <- grid::rasterGrob(img, interpolate = TRUE)
      
      grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
      
      # Add image with space above for title
      grid::pushViewport(grid::viewport(y = unit(0.45, "npc"), height = unit(0.85, "npc")))
      grid::grid.draw(grob_img)
      grid::popViewport()
      
      # Title
      grid::grid.text(basename(page_files[i]),
                      y = unit(0.98, "npc"),
                      just = "top",
                      gp = grid::gpar(fontsize = title_fontsize))
      grid::popViewport()
    }
    
    grid::popViewport()
  }
  
  dev.off()
  
  message("Gallery saved as PDF: ", normalizePath(output_pdf))
}


