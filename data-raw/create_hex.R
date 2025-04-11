#' create hex sticker
#'
#' Supply path to image and this package will resize for use 
#' in pkgdown website. 
#' 
#' @param imgurl Full path to image file
#' @param exportPath Path to location of hex image. Final hex should be saved in man/figures folde
#' 
#' @return png file. r
#' this function is highly costume for the current hex image. Please do not change it. 

create_hex <- function(imgurl = NULL, exportPath = here::here("man/figures"), hex_color= "#062825", hex_border="#083F22") {
  
  if(is.null(imgurl)){
    stop("Please provide a path to an image file")
  }
  #imgurl <- here::here("data-raw","testlogo.jpg")
  hexSticker::sticker(imgurl,
                      package = NULL,
                      p_size=24,
                      s_x=0.99,
                      s_y=1,
                      s_width = 1,
                      s_height = 1,
                      p_x = 1,
                      p_y = 1.2,
                      h_size = 1.5,
                      h_fill= hex_color,
                      h_color= hex_border,     #"#bc4700", "#083F22"
                      #angle = 30,
                      spotlight=T,
                      l_x = 1,
                      l_y = 1.75,
                      l_width = 1.5,
                      l_height = 1.5,
                      l_alpha = 0.15,
                      u_x = 1.78,
                      u_y=.7,
                      url = "noaa-edab.github.io/Rpath",
                      u_size = 5,
                      u_color = "#9DF2EF",
                      u_angle = 90,
                      white_around_sticker = FALSE,
                      #filename = here::here("man/figures", "logo.png"),
                      #filename = here::here("data-raw", "logo.png")
                      filename = paste0(exportPath, "/logo.png")
                      )
  
}
