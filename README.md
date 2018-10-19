# bioanalyzeR

Simple R functions for importing and graphing electrophoresis data from an Agilent 2100 Bioanalyzer or 2200 TapeStation.

## `read.bioanalyzer`
Reads data from one or more Bioanalyzer XML files (see below for how to export those) into a tall data frame containing all useful information. Each row is a single measurement so there are thousands of rows. The columns are:

* `batch` - the batch (instrument run) of the sample, from the file name
* `index` - the loading order of the sample
* `name` - the name of the sample (as a factor, whose levels are in the order observed in the input data)
* `time` - the time this data point was measured
* `fluorescence` - the fluorescence reading at this time
* `length` - the estimated length, in nucleotides, of the molecules being read at this time
* `molarity` - the estimate concentration, in picomolar (pM), of the molecules being read at this time

`length` and `molarity` will be NA if they are outside the range of interpolation from the standard curve. By default the standard curve is fit from the ladder according to linear interpolation between the bands, which is very simplistic but seems to match Agilent's approach. Optionally, set `fit = "spline"` to fit the curve more smoothly according to a parametric model.

## `read.tapestation`
Reads data from one or more TapeStation runs, requiring a list of metadata files and a second list of corresponding gel image files in the same order (see below). Reports a tall data frame similar to that from `read.bioanalyzer`, with two differences:
* Instead of `time` there is `distance` and `relative.distance` (normalized between marker bands).
* Instead of `index` there is `well.number`, which is also parsed into `well.col` and `well.row`.

## `plot.electropherogram`
Generates a nice `ggplot` object from a data frame of the type generated by `read.bioanalyzer` and `read.tapestation`.

Each sample is shown in a separate facet, except the ladder is not shown by default (`include.ladder`). `x` and `y` can be set to graph different variables, and `geom` can be any applicable ggplot2 geom function. By default, the vertical axis is scaled independently for each facet; the `scales` argument is passed to `facet_wrap`. You can override these settings by adding another `facet_wrap` etc.

## How to export Bioanalyzer data
In the 2100 Expert software, open your data file (`.xad`) in the "Data" context. Select "File->Export..." from the top menu. Check the "Export to XML" box and no others. Click "Export" and then save the file wherever you like. This will create an XML file that can be read by the `read.bioanalyzer` function.

## How to export TapeStation data
In the TapeStation Analysis Software, open your data file (`.D1000`, `.HSD1000`, `.RNA`., `.gDNA`, etc.). You need to export both the metadata (in XML format) and gel image (PNG).

### Metadata XML
Select "File->Export Data->Export to XML". You do not need to export the gel image or individual EPG images at this point. Select your destination and then click "Export" to save the file.

###Gel image PNG
It is important to follow these unusual directions carefully!

1. In the "Home" tab (top), verify that there is a ladder lane ("Electronic Ladder" is okay) and that the markers are correctly identified in every sample (except failed lanes, which are okay).
1. Select the "Gel" context (top left button).
1. Select "Show All Lanes" if the button is not grayed out.
1. Unselect "Aligned", "Scale to Sample", and "Scale to MW Range" if it is not grayed out.
1. Leave the contrast slider in the middle.
1. Maximize the window and drag the lower end of the gel image area to make it as tall as possible. At this point you should see all lanes from the run, with the marker bands present but unaligned.
1. Right-click on a lane near the left end of the gel (it doesn't matter which) and you will see a context menu.
1. Move your cursor over "Snapshot" but **also over a lane to the right of the one you right-clicked on**.
1. Left-click "Snapshot". This will copy the gel image to your system clipboard, but you should see also see the newly selected lane become highlighted in light blue.
1. Open any image editor (e.g. Paint) and paste the image from your clipboard. You should see one lane highlighted in light blue. If not, try taking the snapshot again.
1. Save the gel image as a file in PNG format, preferably with the same name as the XML file (but the `.png` extension).

## Why does `read.tapestation` need a gel image?
The TapeStation software cannot export the raw fluorescence data in any other format; the XML only contains metadata (sample ID, peak list, etc.). The data would be easily accessible inside the `.D1000`, `.HSD1000`, etc. file formats except they are encrypted. I have asked Agilent about this and they stated that they do not want the data to be accessible to users.

