#' create hex sticker
#'
#' Supply path to image and this package will resize for use 
#' in pkgdown website. 
#' 
#' @param imgurl Full path to image file
#' @param exportPath Path to location of hex image. Final hex should be saved in man/figures folde
#' 
#' @return png file. r
#'

create_hex <- function(imgurl = NULL, exportPath = here::here("man/figures")) {
  
  if(is.null(imgurl)){
    stop("Please provide a path to an image file")
  }
  #imgurl <- here::here("data-raw","testlogo.jpg")
  hexSticker::sticker(imgurl,
                      package = NULL,
                      p_size=24,
                      s_x=1,
                      s_y=1,
                      s_width=1.28,
                      p_x = 1,
                      p_y = 1.2,
                      h_size = 1,
                      h_fill="#FFFFFF",
                      h_color="#bc4700",
                      #angle = 30,
                      spotlight=T,
                      l_x = 1.5,
                      l_y = 1.1,
                      l_width = 3,
                      l_height = 3,
                      l_alpha = 0.5,
                      u_x = .99,
                      u_y=.05,
                      url = "noaa-edab.github.io/Rpath",
                      u_size = 5.5,
                      u_color = "#FFFFFF",
                      white_around_sticker = T,
                      #filename = here::here("man/figures", "logo.png"),
                      #filename = here::here("data-raw", "logo.png")
                      filename = paste0(exportPath, "/logo.png")
                      )
  
}
