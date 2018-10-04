# library for geo-spatial related R code

# note that R geosphere and sp packages contain many
# related functions.  E.g., geosphere contains distance
# matrix function 'distm'

# there is a conflict between the representation here of
# a bounding box and the representation in package sp.
# In sp the rows are 'x','y' and the columns are 'min','max'.
# Here, to make bounding boxes a special case of a location
# matrix, cols are 'x','y' and rows are 'min','max'

library(plotrix)
library(ggmap)
library(grid)
library(geosphere)

# convert kilometers to degrees (only approximate)
km_to_deg = function(x) { return(x/111.0)}
deg_to_km = function(x) { return(x * 111.0)}

km_to_mi  = function(x) { return(x*0.6214)}

#
# some utility functions
#

# remove "consecutive duplicates" from a vector
remove_conseq_dups = function(v) {
  dups = (c(NA, v[1:(length(v)-1)]) == v) # compare v to an "offset" version of itself
  dups[is.na(dups)] = F
  return(v[!dups])
}

# remove "consecutive duplicate" BSIDs from a vector of cells
remove_conseq_dup_bsids = function(v) {
  v1 = sapply(v, function(x) substr(x, 1, nchar(x)-1))
  dups = (c(NA, v1[1:(length(v1)-1)]) == v1) # compare v to an "offset" version of itself
  dups[is.na(dups)] = F
  return(v[!dups])
}

# return sequence of long/lat values that is like the input
# except with duplicates removed
#
remove_dup_latlongs = function(lseq) {
  dseq = diff(lseq)
  i = intersect(which(dseq[,1] == 0), which(dseq[,2] == 0))
  if (length(i) > 0) {
    result = lseq[-i,] 
  } else {
    result = lseq
  }
  return(result)
}

#
# basic geometry
#

# generate points for a circle
circleFun <- function(center = c(0,0), radius = 1, npoints = 20){
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + radius * cos(tt)
  yy <- center[2] + radius * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# distance (in km) on surface of earth between (long1,lat1) and (long2,lat2)
# (note: R package "lmap" could be used for this)
geodetic.distance = function(point1, point2) {
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))  
  d <- acos(pmin(pmax(d,-1.0),1.0))
  R*d
}

#
# functions for bounding boxes and coordinates
#

# some function below use concept of "bounding box" and
# "coords".  These are specified as matrices and vectors,
# respectively, with the following format:
#
# bounding box
#    left    bottom
#    top     right
#
# coords
#    left  bottom  right  top
#
# (i.e.  min_long  min_lat  max_long  max_lat)
#

bbox_to_coords = function(bb) {
  v = c(bb[1,],bb[2,])
  names(v) = c("left","bottom","right","top")
  return(v)
}

coords_to_bbox = function(coords) { 
  bb = matrix(coords, byrow=TRUE, ncol=2)
  rownames(bb) = c("x","y")
  colnames(bb) = c("min", "max")
  return(bb)
}

# given long/lat matrix m, return smallest bounding
# box that encloses the points in m
# (note: this is similar to function 'bbox' in package sp)
matrix_to_bbox = function(m) {
  bb = matrix(c(min(m[,1]), min(m[,2]), max(m[,1]), max(m[,2])), byrow=TRUE, ncol=2)
  colnames(bb) = c("x","y")
  rownames(bb) = c("min", "max")
  return(bb)
}
  
# scale given bounding box in each direction by scaling
# factor s (s=0.15 means 15% added to each side)
scale_bbox = function(bb, s) {
  dx = bb[2,1] - bb[1,1]
  dy = bb[2,2] - bb[1,2]
  bbs = bb + matrix(c(-s*dx,s*dx,-s*dy,s*dy),ncol=2)
  return(bbs)
}

# scale given bounding box in each direction by scaling
# factor s (s=1.15 means 15% added to each side)
scale_bbox_mult = function(bb, s) {
  center = c(sum(bb[,1])/2, sum(bb[,2])/2)
  dx = (bb[2,1] - bb[1,1])/2
  dy = (bb[2,2] - bb[1,2])/2
  bbs = matrix(c(center[1]-s*dx, center[1]+s*dx, center[2]-s*dy, center[2]+s*dy), ncol=2)
  return(bbs)
}

# extend a bbox by given number of kms on each side
extend_bbox = function(bb, margin) {
  dx_deg = bb[2,1] - bb[1,1]
  dy_deg = bb[2,2] - bb[1,2]
  dx_km = geodetic.distance(bb[1,], c(bb[2,1], bb[1,2]))
  dy_km = geodetic.distance(bb[1,], c(bb[1,1], bb[2,2]))
  x_deg_per_km = dx_deg/dx_km
  y_deg_per_km = dy_deg/dy_km
  x_margin = x_deg_per_km * margin
  y_margin = y_deg_per_km * margin
  bb_new = bb + matrix(c(-x_margin, x_margin, -y_margin, y_margin), ncol=2)
  return(bb_new)
}

# get bound box for a matrix and scale it
bound_box = function(m, adjust=0) scale_bbox(matrix_to_bbox(m), adjust)

# return center point of a bounding box as a (x,y) vector
bbox_center = function(bb) c(mean(bb[,1]), mean(bb[,2]))

# return bounding box defined by center long/lat and radius (in meters)
point_to_bbox = function(long, lat, radius) {
  p = c(long,lat)
  top = destPoint(p, 0, radius)
  right = destPoint(p, 90, radius)
  bottom = destPoint(p, 180, radius)
  left = destPoint(p, 270, radius)
  bb = matrix(c(left[1],right[1],bottom[2],top[2]), ncol=2)
  colnames(bb) = c("x","y")
  rownames(bb) = c("min", "max")
  return(bb)
}

# deprecated
local_area = point_to_bbox

# return the bounds (min_lon, min_lat, max_lon, max_lat)
# of a square enclosing the given points, with an adjustment
# relative to the length of a side of the square that just 
# fits the given points
# (m is a long/lat matrix of at least two rows)
bound_square = function(m, adjust=1.0) {
  min_lat  = min(m[,2])
  min_long = min(m[,1])
  max_lat  = max(m[,2])
  max_long = max(m[,1])
  d = geodetic.distance(c(min_long, min_lat), c(max_long, max_lat)) * 1000
  d = d * adjust
  center = c(min_long + max_long, min_lat + max_lat)/2
  la = point_to_bbox(center[1], center[2], d/2)
  return(la)
}

# return TRUE iff bbox1 contains bbox2
box_contains = function(bbox1, bbox2) {
  # min of bbox2 must be greater than min of bbox1
  # max of bbox2 must be less than max of bbox1
  return(all(bbox2[1,] >= bbox1[1,]) && all(bbox2[2,] <= bbox1[2,]))
}

# return a matrix containing the rows of long/lat matrix m
# that lie within the given bounding box
filter_points = function(m, bbox) {
  m1 = m[m[,1] >= bbox[1,1] & m[,1] <= bbox[2,1] &
           m[,2] >= bbox[2,1] & m[,2] <= bbox[2,2],]
  return(m1)
}

#
# functions for getting maps
#
# to obtain a google map via the google api, you provide a
# location and a zoom level.  We want to specify maps by a
# bounding box, not a location and zoom level.  ggmap's get_map
# function tries to estimate zoom level, but does a poor job,
# often giving maps that don't include the bounding box given.
# By making lots of get_map calls with various parameters, I was
# able to reverse engineering how the zoom value works and
# do a better job than get_map (see file map-zoom.R)

# get zoom as a function of x (longitude) range in degrees
zoom_x = function(xr) 9.815 - log2(xr)

# get zoom as a function of y (latitude) range in degrees,
# and y (latitude) in degrees
zoom_y = function(yr, y) {
  a = 5.584E-2
  b = -9.424E-5
  c = -6.255E-6
  r = log2((a + b*y + c*y^2)/yr) + 14.01
}

# get zoom value needed to get map that encloses the area
# in the given bounding box
get_zoom = function(bb) {
  xrange = bb[2,1] - bb[1,1]
  yrange = bb[2,2] - bb[1,2]
  ycenter = (bb[2,2] + bb[1,2])/2
  xzoom = ceiling(zoom_x(xrange))
  yzoom = ceiling(zoom_y(yrange, ycenter))
  zoom = min(xzoom, yzoom)
  return(zoom)
}

# get a map that encloses the given bounding box
box_map = function(bb, retry=TRUE) {
  # it would be good to really get the zoom calculation right
  z = get_zoom(scale_bbox(bb, 0.2))
  map = get_map(location=bbox_center(bb), zoom=z, maptype="road", color="color")
  # zoom out and try again if necessary
  if (retry && !box_contains(map_box(map), bb)) {
    map = get_map(location=bbox_center(bb), zoom=z-1, maptype="road", color="color")
  }
  return(map)
}

# get plot of map, given long/lat matrix (which may be just a bounding box)
get_map_plot = function(m, adjust=0.0) {
  map = box_map(bound_box(m, adjust))
  p = ggmap(map)
  return(p)
}

# get map from lat/long vectors
latlongs_map = function(lats, longs) {
  return(box_map(bound_box(cbind(longs, lats))))
}

# given a map, get the bounding box
map_box = function(map) {
  v = attr(map, "bb")
  bb = matrix(c(v$ll.lon, v$ll.lat, v$ur.lon, v$ur.lat), byrow=TRUE, ncol=2)
  return(bb)
}
               
#
# functions for plotting on maps
#

# plot points on a map
# points is a data frame including columns "lat" and "long"
# coords is a vector c(min_long, min_lat, max_long, max_lat)
plot_map_points = function(points, map, coords, col="red") {
  # filter points outside coordinates
  fpoints = points[points$long >= coords[1] & points$long <= coords[3] &
                     points$lat >= coords[2] & points$lat <= coords[4],]
  
  # layer 1 is the map
  p = ggmap(map)   
  
  # layer 2 are the road points
  p = p + geom_point(data=fpoints, aes(x=long,y=lat), color=col, size=2)
  
  return(p)
}

# plot a graph on a map
# nodes is a data frame including columns "lat" and "long"
# edges is a data frame including columns "from.lat", "from.long", "to.lat", and "to.long",
#  and optionally "freq" (edge weighting value)
# coords is a vector c(min_long, min_lat, max_long, max_lat)
plot_map_graph = function(nodes, edges, map, coords, label_nodes = F) {
  # filter nodes outside coordinates
  fnodes = nodes[nodes$long >= coords[1] & nodes$long <= coords[3] &
                   nodes$lat >= coords[2] & nodes$lat <= coords[4],]
  
  # layer 1 is the map
  p = ggmap(map)   
  
  # layer 2 is the edges
  if ("freq" %in% names(edges)) {
    p = p + geom_segment(data=edges, aes(x=from.long, xend=to.long, y=from.lat, yend=to.lat, alpha=freq))
  } else {
    p = p + geom_segment(data=edges, aes(x=from.long, xend=to.long, y=from.lat, yend=to.lat))
  }
  
  # layer 3 is the nodes
  p = p + geom_point(data=fnodes, aes(x=long,y=lat), color="red", size=2)
  
  # layer 4 is the node labels
  if (label_nodes) {
    p = p + geom_text(data=fnodes, aes(x=long,y=lat,label=bsid), vjust=-2)
  }
  
  return(p)
}

# deprecated
extend_map_with_cells = function(p, cells, latlong, col="red") {
  p = add_cells_to_map(p, cells, latlong)
  return(p)
}

# given a vector of cells, return a "nodes" data frame
# containing fields "bsid", "lat", "long"
cells_to_nodes = function(path, latlong) {
  nodes = merge(data.frame(bsid=unique(path)), latlong, by="bsid")
  return(nodes)
}

# given a vector of cells, return an "edges" data frame
# containing fields "from", "to", "from.lat", "from.long", "to.lat", "to.long"
# (duplicate consecutive bsids in input path are removed)
path_to_edges = function(path, latlong) {
  path1 = remove_conseq_dup_bsids(path)
  n = length(path1)
  edges = data.frame(from=path1[1:(n-1)], to=path1[2:n])
  edges$seq_num = 1:(n-1)
  
  edges = merge(edges, latlong, by.x="from", by.y="bsid")
  names(edges)[names(edges) == "lat"] = "from.lat"
  names(edges)[names(edges) == "long"] = "from.long"
  
  edges = merge(edges, latlong, by.x="to", by.y="bsid")
  names(edges)[names(edges) == "lat"] = "to.lat"
  names(edges)[names(edges) == "long"] = "to.long"
  
  edges = edges[order(edges$seq_num),]
  edges = edges[,c("from", "to", "from.lat", "from.long", "to.lat", "to.long")]
  
  return(edges)
}

# given a "nodes" data frame, return an "edges" data frame
nodes_to_edges = function(nodes) {
  n = nrow(nodes)
  edges = data.frame(from.lat=nodes[1:(n-1),]$lat, 
                     from.long=nodes[1:(n-1),]$long,
                     to.lat=nodes[2:n,]$lat,
                     to.long=nodes[2:n,]$long)
  return(edges)
}

# plot a vector of bsids as a graph
plot_cell_path = function(path, latlong, label_nodes=F) {
  nodes = cells_to_nodes(path, latlong)
  edges = path_to_edges(path, latlong)
  coords = bound_box(nodes$lat, nodes$long, adjust=2)
  map = coords_map(coords)
  plot_map_graph(nodes, edges, map, coords, label_nodes)
}

# plot a vector of bsids as points
plot_cells = function(cells, latlong) {
  nodes = path_to_nodes(cells, latlong)
  coords = bound_box(nodes$lat, nodes$long, adjust=2)
  map = coords_map(coords)
  p = ggmap(map)
  p = add_nodes_to_map(nodes)
  return(p)
}

# add a cell sequence to a plot
add_cell_path = function(p, path, latlong) {
  nodes = cells_to_nodes(path, latlong)
  edges = path_to_edges(path, latlong)
  p = add_nodes_to_map(p, nodes)
  p = add_edges_to_map(p, edges)
  return(p)
}

add_nodes_to_map = function(p, nodes, col="red", sz=2) {
  p = p + geom_point(data=nodes, aes(x=long,y=lat), color=col, size=sz)
  return(p)
}

add_edges_to_map = function(p, edges, a=1.0, col="black", sz=1, arrows=FALSE) {
  if (arrows) {
    p = p + geom_segment(data=edges, color=col, size=sz, alpha=a, arrow=arrow(angle=35, ends="first", length=unit(0.25, "cm")), 
                         aes(x=to.long, xend=from.long, y=to.lat, yend=from.lat))
  } else {
    p = p + geom_segment(data=edges, aes(x=from.long, xend=to.long, y=from.lat, yend=to.lat))
  }
  return(p)
}

#
# the following functions work purely with locations
#

# add long/lat path to given plot p
add_path = function(p, m, col="black", sz=1.0, a=1.0, lty=1, arrows=TRUE) {
  lseq = as.data.frame(remove_dup_latlongs(m))
  n = nrow(lseq)
  edges = cbind(lseq[1:(n-1),], lseq[2:n,])
  names(edges) = c("from.long", "from.lat", "to.long", "to.lat")
  
  if (arrows) {
    # from/to are written backwards, with ends="first", to get arrow heads to work
    # (see http://stackoverflow.com/questions/18747706/adding-a-bunch-of-arrows-to-a-ggmap)
    p = p + geom_segment(data=edges, color=col, size=sz, linetype=lty, alpha=a, arrow=arrow(angle=35, ends="first", length=unit(0.25, "cm")), 
                         aes(x=to.long, xend=from.long, y=to.lat, yend=from.lat))
  } else {
    p = p + geom_segment(data=edges, color=col, size=sz, alpha=a, aes(x=from.long, xend=to.long, y=from.lat, yend=to.lat))
  }
  return(p)
}

# add long/lat path to given plot p, no arrowheads, dotted line
add_dotted_path = function(p, m, col="black") {
  lseq = remove_dup_latlongs(m)
  n = nrow(lseq)
  edges = data.frame(cbind(lseq[1:(n-1),], lseq[2:n,]))
  names(edges) = c("from.long", "from.lat", "to.long", "to.lat")
  
  # from/to are written backwards, with ends="first", to get arrow heads to work
  # (see http://stackoverflow.com/questions/18747706/adding-a-bunch-of-arrows-to-a-ggmap)
  p = p + geom_segment(data=edges, color=col, size=1, linetype=3, aes(x=to.long, xend=from.long, y=to.lat, yend=from.lat))
  return(p)
}

# add long/lat points to given plot p
# (ggplot params: col=color, shp=shape, sz=size, a=alpha)
# shape is usual R shape codes
add_points = function(p, m, col="black", shp=16, sz=1, a=1.0, labels=c()) {
  points = data.frame(long=m[,1], lat=m[,2])
  p = p + geom_point(data=points, aes(x=long,y=lat), color=col, shape=shp, size=sz, alpha=a)
  if (length(labels) == nrow(points)) {
    points$lab = labels
    # p = p + geom_text(data=points, color=col, aes(x=long,y=lat,label=lab), vjust=-2)
    p = p + geom_text(data=points, aes(x=long,y=lat,label=lab), size=6, vjust=-1.5)
  }
  return(p)
}

# add orientation line at each given point, per given angle
# in matrix m, column 1 is long, 2 is lat, 3 is angle
add_orientation = function(p, m, col="gray30", len_scale=0.05) {
  # scale length of line to area of plot
  ranges = ggplot_build(p)$panel$ranges[[1]]
  yrange = ranges$y.range
  xrange = ranges$x.range
  minloc = c(xrange[1],yrange[1])
  maxloc = c(xrange[2],yrange[2])
  dist = distCosine(maxloc, minloc)
  len = dist*len_scale
  
  df = data.frame(long=m[,1], lat=m[,2], azimuth=m[,3])
  dests = apply(m, 1, function(x) destPoint(c(x[1],x[2]), x[3], len))
  df$long_end = dests[1,]
  df$lat_end = dests[2,]
  p = p + geom_segment(data=df, aes(x=long,y=lat,xend=long_end,yend=lat_end), linetype=1, color=col)
  return(p)
}

# return a plot showing the given sequence of locations
# lseq is a matrix with longitude in first column, latitude in second
# (deprecated)
plot_loc_seq = function(lseq, col="black") {
  coords = bound_box(lseq[,2], lseq[,1], adjust=1)
  map = coords_map(coords)
  p = ggmap(map)
  p = add_loc_seq_to_plot(p, lseq, col)
  return(p)
}

# add a sequence of locations to a plot
# (deprecated)
add_loc_seq_to_plot = function(p, lseq, col="red") {
  lseq = remove_dup_latlongs(lseq)
  n = nrow(lseq)
  edges = data.frame(cbind(lseq[1:(n-1),], lseq[2:n,]))
  names(edges) = c("from.long", "from.lat", "to.long", "to.lat")

  # from/to are written backwards, with ends="first", to get arrow heads to work
  # (see http://stackoverflow.com/questions/18747706/adding-a-bunch-of-arrows-to-a-ggmap)
  p = p + geom_segment(data=edges, color=col, size=1.5, arrow=arrow(ends="first", length=unit(0.3, "cm")), aes(x=to.long, xend=from.long, y=to.lat, yend=from.lat))
  return(p)
}

# given a plot and the map used to create the plot (via ggmap),
# add a scale legend to the plot
add_scale = function(p, map) {
  bb = map_box(map)
  # left to right edge distance in miles, for map
  xdist = km_to_mi(geodetic.distance(bb[1,], bb[2,]))
  xincr = xdist/10
  scales = c(0.1, 0.25, 0.5, 1, 5, 10, 50, 100)
  xscale = ifelse(xincr >= max(scales), max(scales), min(scales[xincr < scales]))
  x_scale_left  = bb[1,1] + 0.1*(bb[2,1] - bb[1,1])
  x_scale_right = x_scale_left + (xscale/xdist)*(bb[2,1] - bb[1,1])
  y_scale_loc = bb[1,2] + 0.1*(bb[2,2] - bb[1,2])
  sline = data.frame(from.lon=x_scale_left, to.lon=x_scale_right, from.lat=y_scale_loc, to.lat=y_scale_loc)
  p = p + geom_segment(data=sline, size=1.5, aes(x=from.lon, xend=to.lon, y=from.lat, yend=to.lat))
  p = p + annotate("text", label=paste(xscale,"miles"), x=(x_scale_left + x_scale_right)/2, y=y_scale_loc-.02*(bb[2,2]-bb[1,2])) 
  return(p)
}





