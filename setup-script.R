# Script to open mapdeck tabs for actual presentation
library (mapdeck)
library (colourvalues)
set_token (Sys.getenv ("MAPBOX_TOKEN"))

plot1map <- function (f, scale = 10)
{
    if (is.character (f)) # name of file
    {
        f <- readRDS (f)
        f <- dodgr::merge_directed_flows (f)
    }
    library (mapdeck)
    f$flow <- scale * f$flow / max (f$flow)
    mapdeck (style = 'mapbox://styles/mapbox/dark-v9') %>%
        add_line (data = f,
                  origin = c("from_lon", "from_lat"),
                  destination = c("to_lon", "to_lat"),
                  stroke_colour = "flow",
                  stroke_width = "flow",
                  palette = colour_palettes("colorRamp")[7]) # green2red
}

plot1map ("data/dodgr-flows-aggregate-hh.Rds", scale = 10)
plot1map ("data/dodgr-flows-disperse-hh.Rds", scale = 5)
