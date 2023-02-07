plot_flows <- function(
    y, # value of the flow
    index_o, # label of the origin
    index_d,  # label of the destination
    geo_site,  # the geographic coordinates of the sites: polygon or points of class sf
    column_name, # the column name of the regions in geo_site sf object and/or in xy_sf
    add = F,   # if add = T, the arcs and only the arcs will printed on the existing device
    type = "barchart", # if T plot barplot for outloflows/inflows and intra flows
    xy_sf = NULL, # possibility to give different coordinates of the starting and ending flows (by default, it corresponds to the centroid of geo_site)
    select_o = NULL, # the origin sites to print
    select_d = NULL, # the destination site to print
    select_bar = NULL, # the sites to print the barplot
    select_flows = NULL, # a vector of logical of size N indicating the flows to be printed
    col_site = NULL, # a vector of colors with attribute names equal to the name of the sites
    transparency = list(alpha = c(1, 1), 
                        range = c(min(y), max(y))), # two scalars to define the values of the minimum and maximum transparency
    width_arc = 1, # width of the arc that is also equals to the height of the bar
    size_head_arc = 1, # the size of the ending and starting of the arrows
    reduce_arc = 0, # the percentage of reduction of the flow
    percent_arc = 1, # the percentage of flows to represent
    width_bar = 1, # the width of the barplot
    lwd_bar = 0.2,   # the lwd of the bar
    percent_bar = 1, # the percentage of barplot to print
    col_direction = "outflows", # a character among "outflows", "inflows", "none"; by default, the outflows will be plotted in color
    col_geometry = "lightgrey",  # color inside the polygons
    border_geometry = "white", # color of the contours 
    lwd_geometry = 0.2, # width of the contours,
    print_names = F, # print the names of the sites on the map
    size_names = 1, # size of the names if printed
    print_legend = F, # legend for the width of the flows
    pos_legend = c("topright", "topleft", "bottomright", "bottomleft"), # position of the legend
    values_legend = NULL, # values of the flows to be printed on the legend
    size_legend = 1, # size of the characters in the legend
    horiz_legend = F, # by default vertical
    title_legend = NULL, # the name of the legend
    digit_legend # number of digits to keep in legend
    ) {
  
  ############### Initialisation   
  # number of flows 
  N <- length(y)
  # index of the origin
  O <- levels(unique(factor(index_o)))
  n_o <- length(O)
  D <- levels(unique(factor(index_d)))
  n_d <- length(D)
  # number of unique site in S
  S <- union(O, D)
  n <- length(S)
  # the sites selected to be printed
  if (is.null(select_o)) 
    select_o <- O
  if (is.null(select_d)) 
    select_d <- D 
  if (is.null(select_bar)) 
    select_bar <- S 
  # the flows to be printed
  if (is.null(select_flows)) {
    select_flows <- !logical(N)
  }
  # Transform the polygons in points 
  sf_use_s2(FALSE)
  if (is.null(xy_sf)) {
    xy_sf <- st_centroid(geo_site)
  }  
  # graphical parameters the bbox aims to find 
  # * the width of the arc
  # * the size of the arc
  my_bbox <- st_bbox(geo_site)
  my_side_x <- my_bbox$xmax - my_bbox$xmin
  my_side_y <- my_bbox$ymax - my_bbox$ymin
  my_side <- max(my_side_x, my_side_y)
  # the width of the arc is relative to the sum of the outflows
  my_sum_flows <- sum(y)
  width_arc <- y / my_sum_flows * my_side / 5 * width_arc 
  # The width of the bar is proportional to the bbox of the x-axis
  width_bar <- my_side_x / 100 * width_bar
  # The size of the head of the arcs
  size_head_arc <- my_side / 100 * size_head_arc
  # Get the coordinates of the centroid or the given spatial object
  coords_xy <- st_coordinates(xy_sf)

  # the coordinates will depend on the presence or not of the barplot
  my_coords <- st_coordinates(xy_sf)
  
  # give the row names
  row.names(my_coords) <- unlist(st_drop_geometry(xy_sf)[, column_name], use.names= F)

  # separate origin and destination
  my_coords_O <- my_coords[O, ]
  my_coords_D <- my_coords[D, ]
  
  # the region are ordered s.t. the bars in the barplot are ordered from the North to the South
  region_ordered_o <- row.names(my_coords_O[order(my_coords_O[, 2]), ])
  region_ordered_d <- row.names(my_coords_D[order(my_coords_D[, 2]), ])
  
  # colors : we create a vector of colors
  if((col_direction == "outflows" & !all(O %in% names(col_site))) |
     (col_direction == "inflows" & !all(D %in% names(col_site)))) {
    
    q4 <- colorspace::qualitative_hcl(ifelse(col_direction == "outflows", n_o, n_d),
                                      palette = "Dark 3")
    if (col_direction == "outflows") {
      names(q4) <- region_ordered_o 
      vec_col <- q4[index_o]
    } else {
      names(q4) <- region_ordered_d 
      vec_col <- q4[index_d]
    }
  } else { 
    if (col_direction == "none") {
      if (is.null(col_site)) {
        q4 <- rep("grey", n)
        names(q4) <- S
      } else {
        q4 <- rep(col_site[1], n)
        names(q4) <- S
      }
      vec_col <- rep(q4[1], N)
    } else {
      q4 <- col_site
      if (col_direction == "outflows") {
        vec_col <- q4[index_o]
      } else {
        vec_col <- q4[index_d]
      }
    }
  }
  # use the transparency in case all flows are not selected
  vec_col_bar <- ifelse(select_flows, vec_col, adjustcolor(vec_col, 0.25))
  vec_col_bar_border <- ifelse(select_flows, vec_col, adjustcolor(vec_col, 0.01))
  
  # the location of the intra flow
  ind_intra <- index_o == index_d
  
  #############################
  # selection of the sites
  # we select the biggest in terms of volume (intra + outflows + inflows)
  total_flows <- rep(0, length(S))
  names(total_flows) <- S
  sum_outflows <- sapply(split(x = y[!ind_intra], f = index_o[!ind_intra]), sum)
  sum_inflows <- sapply(split(x = y[!ind_intra], f = index_d[!ind_intra]), sum)
  total_flows[names(sum_outflows)] <- total_flows[names(sum_outflows)] + sum_outflows
  total_flows[names(sum_inflows)] <- total_flows[names(sum_inflows)] + sum_inflows
  if (length(which(ind_intra)) > 0) {
    sum_intra <-  sapply(split(x = y[ind_intra], f = index_d[ind_intra]), sum)
    total_flows[names(sum_intra)] <- total_flows[names(sum_intra)] + sum_intra
  }
  selected_bar <- intersect(select_bar, S[total_flows >= quantile(total_flows, 1 - percent_bar)])

  if (type == "barchart") {
    #############################
    # initialisation
    x_out <- y_out <- x_in <- y_in <- numeric(N)
    # stock the cumulative
    my_cum_sum_out <- rep(0, n_o)
    names(my_cum_sum_out) <- region_ordered_o
    my_cum_sum_in <- rep(0, n_d)
    names(my_cum_sum_in) <- region_ordered_d
    # stock the barplot of the outflows or inflows
    rect_out <- vector("list", n_o)
    names(rect_out) <- region_ordered_o
    rect_in <- vector("list", n_d)
    names(rect_in) <- region_ordered_d
    
    # the location of the intra flow
    if (length(which(ind_intra)) != 0) {
      x_out[ind_intra] <- x_in[ind_intra] <- my_coords_O[index_o[ind_intra], 1]
      y_out[ind_intra] <- y_in[ind_intra] <- my_coords_O[index_o[ind_intra], 2]
      my_cum_sum_out[index_o[ind_intra]] <- my_cum_sum_in[index_d[ind_intra]] <- 
        width_arc[ind_intra] / 2
    }
    
    # fill from the south to the north
    for(in_country in region_ordered_d) {
      ind <- which(index_d == in_country & index_o != in_country)
      if (length(ind) != 0) {
        # we order
        extract_o <- factor(index_o[ind], levels = region_ordered_o)
        extract_width <- (width_arc[ind])[order(extract_o)]
        ind <- ind[order(extract_o)]
        extract_o <- as.character(sort(extract_o))
        k <- 1
        for (out_country in extract_o) {
          # we fill the outflows
          my_cum_sum_out[out_country] <- my_cum_sum_out[out_country] + width_arc[ind[k]] / 2
          x_out[ind[k]] <- my_coords_O[out_country, 1] + ifelse(my_coords_O[out_country, 1] < my_coords_D[in_country, 1], 0, -1) * width_bar 
          y_out[ind[k]] <- my_coords_O[out_country, 2] + my_cum_sum_out[out_country]
          my_cum_sum_out[out_country] <- my_cum_sum_out[out_country] + width_arc[ind[k]] / 2
          # and we prepare the barplot
          #if (col_direction == "inflows") {
            rect_out[[out_country]] <- rbind(rect_out[[out_country]],
                                           data.frame(
                                             xleft = my_coords_O[out_country, 1] - width_bar,
                                             ybottom = y_out[ind[k]] - width_arc[ind[k]] / 2,
                                             xright = my_coords_O[out_country, 1],
                                             ytop = y_out[ind[k]] + width_arc[ind[k]] / 2,
                                             col = vec_col_bar[ind[k]],
                                             border = vec_col_bar_border[ind[k]])) # q4[in_country]
          #}
          # we fill the inflows
          my_cum_sum_in[in_country] <- my_cum_sum_in[in_country] + width_arc[ind[k]] / 2
          x_in[ind[k]] <- my_coords[in_country, 1] + ifelse(my_coords[out_country, 1] < my_coords[in_country, 1], 0, 1) * width_bar
          y_in[ind[k]] <- my_coords[in_country, 2] + my_cum_sum_in[in_country]
          my_cum_sum_in[in_country] <- my_cum_sum_in[in_country] + width_arc[ind[k]] / 2
          # and we prepare the barplot
          #if (col_direction == "outflows") {
            rect_in[[in_country]] <- rbind(rect_in[[in_country]],
                                           data.frame(
                                             xleft = my_coords[in_country, 1],
                                             ybottom = y_in[ind[k]] - width_arc[ind[k]] / 2,
                                             xright = my_coords[in_country, 1] + width_bar,
                                             ytop = y_in[ind[k]] + width_arc[ind[k]] / 2,
                                             col = vec_col_bar[ind[k]],
                                             border = vec_col_bar_border[ind[k]])) # q4[out_country]))
          #}
          k <- k + 1
        }
      }
    }
    my_coords_flows <- cbind(x_out, y_out, x_in, y_in)
  } else {
    my_coords_flows <- cbind(my_coords[index_o, ], 
                             my_coords[index_d, ], 
                             my_coords[index_o, ], 
                             my_coords[index_d, ])
  }
  
  # plot the arc
  # we order to print them from the smallest to the highest
  my_ord <- rev(order(y))
  ord_index_o <- index_o[my_ord]
  ord_index_d <- index_d[my_ord]
  ord_y <- y[my_ord]
  my_coords_flows <- my_coords_flows[my_ord, ]
  width_arc <- width_arc[my_ord]
  vec_col <- vec_col[my_ord]
  vec_col_bar <- vec_col_bar[my_ord]
  # transparency and colors
  my_transpa <- function(x, transpa = transparency) {
    my_range <- transpa$range
    my_alpha <- transpa$alpha
    K <- length(my_range)
    if (K != length(my_alpha)) 
      stop("In transparency object, range and alpha must have the same length")
    res <- numeric(length(x))
    for (k in 1:length(x)) {
      if (x[k] <= my_range[1]) {
        res[k] <- my_alpha[1]
      } else {
        if (x[k] >= my_range[K]) {
          res[k] <- my_alpha[K]
        } else {
        ind_int <- which(x[k] >= my_range[1:(K-1)] & x[k] < my_range[2:K]) 
        res[k] <- my_alpha[ind_int] + (x[k] - my_range[ind_int]) * 
          (my_alpha[ind_int + 1] - my_alpha[ind_int]) / (my_range[ind_int + 1] - my_range[ind_int])
      }
     }
    }
    return(res)
  }
  alpha_col <- my_transpa(ord_y)
  
  # initialisation of the plot
  if(!add)
    plot(st_geometry(geo_site), 
         col = col_geometry,
         border = border_geometry, 
         lwd = lwd_geometry)
  # plot the arcs
  ind_quantile <- round(quantile(1:N, percent_arc))
  ind_select <- which(ord_index_o %in% select_o & 
                        ord_index_d %in% select_d & 
                        (1:N) <= ind_quantile & 
                        select_flows[my_ord])
  
  for(j in ind_select) {
    if (ord_y[j] != 0)
      polygon(my_arc(my_coords_flows[j, c(1, 2)], 
                     my_coords_flows[j, c(3, 4)], 
                     width_arc = width_arc[j], 
                     width_intra = width_bar,
                     size_head_arc = size_head_arc,
                     reduce_arc = reduce_arc),
              border = NA,
              col = adjustcolor(vec_col[j], alpha_col[j]),
              lwd = 1)
  }
  # add the barplot
  if (type == "barchart") {
    # index of the intra
    index_intra <- which(ord_index_o == ord_index_d & ord_index_o %in% selected_bar)
    # start to plot the bar of the intra flow
    if(length(index_intra) != 0) {
      for(ind in index_intra) {
        polygon(my_arc(my_coords_flows[ind, c(1, 2)], 
                       my_coords_flows[ind, c(3, 4)], 
                       width_arc = width_arc[ind], 
                       width_intra = width_bar,
                       size_head_arc = size_head_arc,
                       reduce_arc = reduce_arc),
                border = "black",
                col = vec_col_bar[ind],
                lwd = lwd_bar)
      }
    }
    # select the sites to be printed
    rect_out <- rect_out[selected_bar]
    rect_in <- rect_in[selected_bar]
    
    # plot the bar of the outflows
     temp <- lapply(rect_out, function(x)
        rect(xleft = x$xleft, 
         ybottom = x$ybottom, 
         xright = x$xright, 
         ytop = x$ytop, 
         col = x$col, 
         border = x$border, #ifelse(rep(col_direction, length(x$col)) == "outflows", x$col, NA), 
         lwd = lwd_bar)
      )
     temp <- lapply(rect_in, function(x)
       rect(xleft = x$xleft, 
            ybottom = x$ybottom, 
            xright = x$xright, 
            ytop = x$ytop, 
            col = x$col, 
            border = x$border, #ifelse(rep(col_direction, length(x$col)) == "inflows", x$col, NA), 
            lwd = lwd_bar)
     ) 
     # add the black contours
     temp <- lapply(rect_in, function(x)
       rect(xleft = x$xleft[1], 
            ybottom = x$ybottom[1], 
            xright = x$xright[1], 
            ytop = x$ytop[length(x$ytop)], 
            col = NA, 
            border = "black", 
            lwd = lwd_bar)
     ) 
     
     temp <- lapply(rect_out, function(x)
         rect(xleft = x$xleft[1], 
              ybottom = x$ybottom[1], 
              xright = x$xright[1], 
              ytop = x$ytop[length(x$ytop)], 
              col = NA, 
              border = "black", 
              lwd = lwd_bar)
       ) 
     
  }
  
  if (print_names) {
    # select the sites to be printed
    text(my_coords[selected_bar, 1], my_coords[selected_bar, 2], 
         selected_bar, pos = 1, cex = size_names)
  }
  
  if (print_legend) {
    diff_x <- (par()$usr[2] - par()$usr[1]) / 20
    diff_y <- (par()$usr[4] - par()$usr[3]) / 20
    x_min <- par()$usr[1] + diff_x
    x_max <- par()$usr[2] - diff_x
    y_min <- par()$usr[3] + diff_y
    y_max <- par()$usr[4] - diff_y 
    if (is.null(values_legend)) {
      # the values of the plot and the true values
      values_legend <- round(c(max(y), max(y) / 2, max(y) / 4), digit_legend)
      width_arc_legend <- c(max(width_arc), max(width_arc) / 2, max(width_arc) / 4)
    }
    # position of the 
    if (pos_legend == "topleft") {
     for (k in 1:length(values_legend)) {
       coord_A <- c(x_min, y_max - k * (y_max - y_min) / 20)
       polygon(my_arc(coord_A, 
                      coord_A, 
                     width_arc = width_arc_legend[k], 
                     width_intra = 2 * width_bar,
                     size_head_arc = size_head_arc,
                     reduce_arc = reduce_arc,
                     side = "right"),
              border = rgb(0.3, 0.3, 0.3),
              col = rgb(0.3, 0.3, 0.3),
              lwd = lwd_bar)
     }
      
      text(x_min + 2 * width_bar, y_max - (1:3) * (y_max - y_min) / 20,
           values_legend, adj = c(0, 0), cex = size_legend)
      
      if (!is.null(title_legend)) {
        text(x_min + width_bar, y_max, title_legend, cex = size_legend)
      }
    }
  }
  
}


# Function that plots an arc
my_arc <- function(A, B, width_arc = 0.1, size_head_arc = 0.1, 
                   width_intra = 0.1, reduce_arc = 0, side = "both") {
  
  # width_arc : width of the arc
  # size_head_arc : size of the arrow head and bottom
  # width_intra : the width of intra flow (the height of intra flow is equal to width_arc)
  # reduce_arc : the percentage of reduccion of the flow
  
  # if A = B, we plot an intra flow
  if (all(A == B)) {
    arc_x <- seq(A[1] - width_intra/2, A[1] + width_intra/2, length.out = 100)
    arc_y <- rep(A[2], 100)
    
    diff_x <- diff(arc_x[c(1, 10)])
    
    arc <- cbind(c(arc_x, 
                   arc_x[100] + diff_x, 
                   rev(arc_x), 
                   arc_x[1] + ifelse(side == "both", -1, 1) * diff_x, 
                   arc_x[1]),
                 c(arc_y + width_arc/2,
                   A[2],
                   arc_y - width_arc/2,
                   A[2],
                   A[2] + width_arc/2))

    return(arc)
  }
  
  
  # find the coordinates of O
  mid_AB <- c(A[1] + B[1], A[2] + B[2]) / 2
  R <- sqrt((B[2] - A[2]) ^ 2 + (B[1] - A[1])^2)
  if ((B[2] != A[2]) & (B[1] != A[1])) { 
    slope <- (B[2] - A[2])/(B[1] - A[1])
    s <- -1 / slope
    k <- sqrt(3) / 2 * R
    # arrow above or below
    if (B[1] > A[1]) {
      O <- mid_AB - c(sqrt(k^2 / (s^2 + 1)), k * s / sqrt(s^2 + 1))
    } else {
      O <- mid_AB + c(sqrt(k^2 / (s^2 + 1)), k * s / sqrt(s^2 + 1))   
    }
  } else {
    # case A and B are on the same horizontal line
    if (B[2] == A[2]) {
      O <- c(mid_AB[1], A[2] + ifelse(B[1] > A[1], 1, -1) * sqrt(3) * R / 2)
    } else {
      O <- c(A[1] + ifelse(B[2] > A[2], 1, -1) * sqrt(3) * R / 2, mid_AB[2])
    }
  }
  # compute the radius between A and B 
  # for one flow 
  # center to compute trigonometry formula
  A_center <- A - O 
  B_center <- B - O 
  angle_A <- ifelse(A_center[2] > 0, acos(A_center[1] / R), 
                    2*pi - acos(A_center[1] / R))
  angle_B <- ifelse(B_center[2] > 0, acos(B_center[1] / R), 
                    2*pi - acos(B_center[1] / R))
  
  # we define the direction on the circle for going from A to B
  if (abs(angle_B - angle_A) >= pi) {
    angle_B <- angle_B + ifelse(angle_A > angle_B, 1, -1) * 2 * pi 
  }

  # reduccion of the angle
  theta <- width_intra / R
  diff_angle <- angle_B - angle_A
  if (2 * theta >= abs(diff_angle)) {
    theta <- abs(diff_angle) / 2
  }
  angle_A <- angle_A + ifelse(diff_angle > 0, 1, -1) * theta * reduce_arc 
  angle_B <- angle_B + ifelse(diff_angle > 0, -1, 1) * theta * reduce_arc 
  
  # we add the start and the end of the arrow
  ending_x <- R * cos(angle_B)
  ending_y <- R * sin(angle_B)

  starting_x <- R * cos(angle_A + size_head_arc / R * ifelse(angle_B > angle_A, 1, -1))
  starting_y <- R * sin(angle_A + size_head_arc / R * ifelse(angle_B > angle_A, 1, -1))

  # we define the sequence of the arrow
  arc_x <- cos(seq(angle_A, angle_B - size_head_arc / R * ifelse(angle_B > angle_A, 1, -1), length.out = 100))
  arc_y <- sin(seq(angle_A, angle_B - size_head_arc / R * ifelse(angle_B > angle_A, 1, -1), length.out = 100)) 
  
  arc <- cbind(c((R - width_arc/2) * arc_x, ending_x, 
                 rev((R + width_arc/2)  * arc_x),
                 starting_x,
                 (R - width_arc/2) * arc_x[1]) + O[1],
               c((R - width_arc/2) * arc_y, ending_y,  
                 rev((R + width_arc/2) * arc_y),
                 starting_y,
                 (R - width_arc/2) * arc_y[1]) + O[2])
  
  return(arc)
}

