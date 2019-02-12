library (osmdata)
library (sf) # for trim_osmdata
library (dodgr)

if (!dir.exists ("data"))
    dir.create ("data")

if (!file.exists ("data/osm-hamburg.Rds"))
{
    bdry <- opq ("hamburg de") %>%
        add_osm_feature (key = "boundary", value = "political") %>%
        add_osm_feature (key = "is_in:state_code", value = "HH") %>%
        osmdata_sf (quiet = FALSE)
    bdry <- bdry$osm_multipolygons
    library (mapdeck)
    set_token (Sys.getenv ("MAPBOX_TOKEN"))
    plot_poly <- function (poly) {
        mapdeck() %>%
            add_polygon (poly,
                         stroke_colour = "#000000",
                         stroke_width = 100,
                         fill_colour = "#22cccc22",
                         fill_opacity = 55)
    }
    plot_poly (bdry)
    # convert multipolygons to polygons and find the enclosing polygon. Note
    # that st_cast does not work here because one of the multipolygons is
    # ill-formed, so have to manually extract the enclosing polygon of each
    # multipolygon object
    bdry <- lapply (bdry$geometry, function (i)
                    st_polygon (list (i [[1]] [[1]]))) %>%
        st_sfc (crs = 4326) %>%
        st_union ()
    plot_poly (st_sf (bdry))

    bp_xy <- st_coordinates (bdry) [, 1:2]
    bbox <- as.vector (t (apply (bp_xy, 2, range)))

    hamburg <- opq (bbox) %>%
        add_osm_feature (key = "highway") %>%
        osmdata_sf (quiet = FALSE) %>%
        trim_osmdata (bb_poly = bp_xy) %>%
        osm_poly2line ()
    # This file is 1/4 GB!
    saveRDS (hamburg$osm_lines, file = "data/osm-hamburg.Rds")
}

if (!file.exists ("data/gtfs-hvv.zip"))
{
    # from http://suche.transparenz.hamburg.de/
    u <- paste0 ("http://daten.transparenz.hamburg.de/",
                 "Dataport.HmbTG.ZS.Webservice.GetRessource100/",
                 "GetRessource100.svc/32aedec8-c69f-4053-a5e8-9b1267eb25de/",
                 "Upload__HVV_Rohdaten_GTFS_Fpl_20190207.zip")
    download.file (u, destfile = "data/gtfs-hvv.zip")
}

if (!file.exists ("data/dodgr-flows-hh.Rds"))
{
    # Read the HVV GTFS data containing the timetables for all HVV Verkehr
    filename <- "data/gtfs-hvv.zip"
    flist <- utils::unzip (filename, list = TRUE)
    for (f in flist$Name)
    {
        fout <- data.table::fread (cmd = paste0 ("unzip -p \"", filename,
                                                 "\" \"", f, "\""),
                                   integer64 = "character",
                                   showProgress = FALSE)
        assign (gsub (".txt", "", basename (f)), fout, envir = .GlobalEnv)
    }
    library (dplyr)
    # Assume importance of stops is reflected in the total number of services
    # arriving at each stop in the entire timetable scheme:
    stop_counts <- stop_times %>%
        group_by (stop_id) %>%
        summarise (total = length (stop_id)) %>%
        right_join (stops) %>%
        filter (!is.na (total)) %>%
        select (stop_id, stop_lat, stop_lon, total) %>%
        as.data.frame ()
    # take the 1,000 most important stops
    cutoff <- sort (stop_counts$total, decreasing = TRUE) [1000]
    stop_counts <- stop_counts [stop_counts$total >= cutoff, ]

    hamburg <- readRDS ("data/osm-hamburg.Rds")
    net <- weight_streetnet (hamburg, wt_profile = "foot")
    net <- net [net$d_weighted < max (net$d_weighted), ]
    v <- dodgr_vertices (net)
    pts <- match_pts_to_graph (v, stop_counts [, 2:3], connected = TRUE)

    f <- dodgr_flows_disperse (net, from = v$id [pts],
                               dens = stop_counts$total, k = 1)
    saveRDS (f, file = "data/dodgr-flows-hh.Rds")
}

if (!file.exists ("data/dodgr-exposure-hh"))
{
    hamburg <- readRDS ("data/osm-hamburg.Rds")

    # functions from flowlayers:
    dist_to_lonlat_range <- function (verts, d = 20)
    {
        xy0 <- c (mean (verts$x), mean (verts$y))
        names (xy0) <- c ("x", "y")
        minf <- function (a, xy0) { abs (geodist::geodist (xy0, xy0 + a) - d) }
        optimise (minf, c (0, 0.1), xy0)$minimum
    }

    disperse_flows <- function (graph, d = 20)
    {
        # This actually works by mapping flows to vertices, and then applying a
        # Gaussian kernel to those points, before re-mapping them back to edges
        v <- dodgr::dodgr_vertices (graph)
        sig <- dist_to_lonlat_range (v, d)

        indxf <- match (v$id, graph$from_id)
        indxt <- match (v$id, graph$to_id)
        indxf [is.na (indxf)] <- indxt [is.na (indxf)]
        indxt [is.na (indxt)] <- indxf [is.na (indxt)]
        v$flow <- apply (cbind (graph$flow [indxf], graph$flow [indxt]), 1, max)

        # warnings are issued here if any (x, y) are duplicated
        xy <- suppressWarnings (spatstat::ppp (v$x, v$y, range (v$x), range (v$y)))
        d <- spatstat::density.ppp (xy, weights = v$flow, sigma = sig, at = "points")
        vsum <- sum (v$flow)
        d <- d * vsum / sum (d)
        indx <- which (d > v$flow)
        v$flow [indx] <- d [indx]
        v$flow <- v$flow * vsum / sum (v$flow)

        # Then map those vertex values back onto the graph
        indxf <- match (graph$from_id, v$id)
        indxt <- match (graph$to_id, v$id)
        fmax <- apply (cbind (v$flow [indxf], v$flow [indxt]), 1, max)
        fsum <- sum (graph$flow)
        fmax <- fmax * fsum / sum (graph$flow)
        indx <- which (fmax > graph$flow)
        # That again imbalances the flow of the graph, so needs to be standardised
        # back to original sum again
        fsum <- sum (graph$flow)
        graph$flow [indx] <- fmax [indx]
        graph$flow <- graph$flow * fsum / sum (graph$flow)

        return (graph)
    }

    net <- weight_streetnet (hamburg, wt_profile = "motorcar")
    net <- net [net$d_weighted < max (net$d_weighted), ]
    net_c <- dodgr_contract_graph (net)
    v <- dodgr_vertices (net_c$graph)
    pts <- sample (v$id, size = 1000)
    message ("aggregating vehicular flows ... ", appendLF = FALSE)
    fv <- dodgr_flows_aggregate (net, from = pts, to = pts,
                                 flow = matrix (1, nrow = 1000, ncol = 1000))
    message ("done\ndispersing ... ", appendLF = FALSE)
    fvd <- disperse_flows (fv, 20) %>%
        merge_directed_flows ()
    fv <- merge_directed_flows (fv)
    message ("done")

    f <- readRDS ("data/dodgr-flows-hh.Rds") %>%
        merge_directed_flows ()
    f <- f [which (f$edge_id %in% fvd$edge_id), ]
    indx1 <- match (f$edge_id, fvd$edge_id)
    indx2 <- match (fvd$edge_id [indx1], f$edge_id)
    f$flow <- f$flow / max (f$flow)
    fvd$flow <- fvd$flow / max (fvd$flow)
    f$flow <- f$flow * fvd$flow [indx1] [indx2]
    saveRDS (f, file = "data/dodgr-exposure-hh.Rds")
}
