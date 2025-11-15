# Dash GUI Geographic Features

## Overview

The PyPopART Dash GUI now includes full support for geographic visualization of haplotype networks. Users can upload metadata with geographic coordinates and overlay their networks on a coordinate space with proper map projections.

## New Features

### 1. Metadata File Upload

**Location**: Upload Data card (Section 1)

The upload card now has two file upload controls:

1. **Sequence File** (required)
   - FASTA, NEXUS, or PHYLIP format
   - Contains the DNA sequences for analysis

2. **Metadata File** (optional)
   - CSV format with `id`, `latitude`, and `longitude` columns
   - Additional metadata columns are also supported
   - Coordinates are automatically extracted and validated

**Example metadata CSV:**

```csv
id,population,location,latitude,longitude
Hap1,PopA,New York,40.7128,-74.0060
Hap2,PopB,London,51.5074,-0.1278
Hap3,PopC,Tokyo,35.6762,139.6503
```

**Upload Status:**

- Success message shows number of sequences with metadata
- Shows number of sequences with valid coordinates
- Error alerts for invalid file formats or parsing errors

### 2. Geographic Layout Option

**Location**: Layout Options card (Section 3)

The layout dropdown now includes:

- Spring (Force-Directed)
- Circular
- Radial
- Hierarchical
- Kamada-Kawai
- **Geographic (requires metadata)** ← NEW

When "Geographic" is selected, additional controls appear:

#### Map Projection

Choose how geographic coordinates are projected:

- **Mercator**: Standard web mapping, preserves angles (default)
- **PlateCarree**: Simple lat/lon mapping, good for global views
- **Orthographic**: 3D globe perspective

#### Zoom Level

Slider from 1-10:

- Lower values (1-3): Global view
- Medium values (4-7): Continental/regional view
- Higher values (8-10): Local/detailed view

### 3. Geographic Mode Visualization

When using geographic layout, the network visualization shows:

**Visual Indicators:**

- "Geographic Mode" badge in top-left corner (blue)
- X-axis labeled "Longitude"
- Y-axis labeled "Latitude"
- Network nodes positioned by their geographic coordinates
- Edges connect nodes showing genetic relationships

**Behavior:**

- If no metadata with coordinates is uploaded, falls back to spring layout
- Network overlay maintains all interactive features (zoom, pan, hover)
- Node sizes still reflect haplotype frequency
- Colors still indicate populations or other attributes

### 4. Workflow

**Step-by-step guide:**

1. **Upload Sequence Data**
   - Click "Select File" in Upload Data section
   - Choose FASTA, NEXUS, or PHYLIP file
   - Wait for success message

2. **Upload Metadata (Optional)**
   - Click "Select Metadata" in Upload Data section
   - Choose CSV file with latitude/longitude columns
   - Verify success message shows coordinate count

3. **Compute Network**
   - Select algorithm (MSN, TCS, MJN, etc.)
   - Configure parameters if needed
   - Click "Compute Network"
   - Wait for computation to complete

4. **Apply Geographic Layout**
   - Go to Layout Options section
   - Select "Geographic" from dropdown
   - Choose map projection (Mercator, PlateCarree, or Orthographic)
   - Adjust zoom level with slider
   - Click "Apply Layout"

5. **View Results**
   - Network tab shows nodes at geographic positions
   - Statistics tab shows network metrics
   - Alignment tab shows sequences

### 5. Technical Details

**Coordinate Requirements:**

- Latitude: -90 to 90 (degrees)
- Longitude: -180 to 180 (degrees)
- Must be in decimal degrees format
- Column names must be "latitude" and "longitude" (case-insensitive)

**Fallback Behavior:**

- If geographic layout selected without coordinates → spring layout
- If some nodes lack coordinates → placed using spring layout for missing nodes
- If coordinates are invalid → error message in metadata upload status

**Data Storage:**

- Metadata stored in `metadata-store` (Dash Store component)
- Coordinates extracted and stored separately
- Geographic mode flag stored in `geographic-mode` store
- Layout positions stored in `layout-store`

### 6. Integration with Existing Features

Geographic visualization works seamlessly with:

- All network construction algorithms (MST, MSN, TCS, MJN)
- Statistics calculations
- Network export functionality
- All distance metrics
- Population coloring (if metadata includes population column)

### 7. Limitations and Notes

**Current Limitations:**

- GUI uses plotly for visualization (not cartopy/folium)
- No actual map background tiles in Dash GUI
- Map projections affect node positioning only
- For full map backgrounds, use CLI `geo-visualize` command

**Why This Design:**

- Dash GUI focuses on interactive network analysis
- Plotly provides excellent interactivity for networks
- CLI tool provides full cartographic features
- Both approaches complement each other

**When to Use Each:**

- **Dash GUI**: Interactive exploration, trying different layouts, network analysis
- **CLI geo-visualize**: Publication-quality maps, presentations, map tiles

### 8. Future Enhancements

Potential future improvements:

- Map tile background in Dash (using plotly's mapbox)
- Multiple coordinate systems per node (for haplotypes in multiple locations)
- Color coding by region
- Great circle paths between related haplotypes
- Animated transitions between layouts

## Example Use Case

**Scenario**: Analyzing mitochondrial DNA from global human populations

1. Upload sequences from 50 individuals across 10 locations worldwide
2. Upload metadata CSV with sampling locations (latitude/longitude)
3. Compute MSN network to see genetic relationships
4. Apply geographic layout with Mercator projection
5. Observe how genetic similarity correlates with geographic distance
6. Export results for publication

**Benefits:**

- Visual pattern recognition
- Geographic vs. genetic distance comparison
- Identifying dispersal routes
- Planning future sampling locations

## Troubleshooting

**Issue**: "Geographic layout shows spring layout instead"

- **Solution**: Ensure metadata CSV has "latitude" and "longitude" columns
- Check metadata upload status for coordinate count
- Verify coordinates are in valid range

**Issue**: "Metadata upload fails"

- **Solution**: Ensure file is CSV format
- Check that first column is "id" matching sequence IDs
- Verify no special characters in coordinate values

**Issue**: "Geographic mode indicator doesn't show"

- **Solution**: Ensure you clicked "Apply Layout" after selecting geographic
- Refresh the page and try again
- Check browser console for JavaScript errors

## API Reference

### Dash Callbacks

**Metadata Upload**

```python
Input: 'upload-metadata' contents
Output: 'metadata-status' (status message), 'metadata-store' (data)
```

**Geographic Options Toggle**

```python
Input: 'layout-select' value
Output: 'geographic-options' style, 'geographic-mode' data
```

**Layout Computation with Geographic**

```python
Input: 'apply-layout-button' clicks, 'network-store' data
State: 'layout-select', 'metadata-store', 'map-projection'
Output: 'layout-store' data
```

**Visualization with Geographic Mode**

```python
Input: 'layout-store', 'network-store', 'geographic-mode'
State: 'metadata-store'
Output: 'network-graph' figure
```

## Conclusion

The Dash GUI now provides a complete interface for geographic haplotype network analysis, making it easy to visualize genetic data in a spatial context. Combined with the CLI tools for publication-quality maps, PyPopART offers a comprehensive solution for phylogeographic analysis.
