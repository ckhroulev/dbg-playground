import shapefile
import Image, ImageDraw

def rasterize(filename, width, height, bg, fg):
  # Read in a shapefile
  data = shapefile.Reader(filename)

  img = Image.new("L", (width, height), bg)
  draw = ImageDraw.Draw(img)

  # domain width and height
  x_max = data.bbox[2]
  x_min = data.bbox[0]
  xdist = x_max - x_min

  y_max = data.bbox[3]
  y_min = data.bbox[1]
  ydist = y_max - y_min

  xratio = width / xdist
  yratio = height / ydist

  def rasterize_poly(points, color):
    pixels = []
    for x,y in points:
      px = int(width - ((x_max - x) * xratio))
      py = int((y_max - y) * yratio)
      pixels.append((px,py))

    draw.polygon(pixels, fill=color)

  def rasterize_shape(shape, bg, fg):
    parts = shape.parts

    if len(parts) == 1:
      rasterize_poly(shape.points, fg)
    else:
      for j in range(0, len(parts)-1):
        if j != 0:
          fg = bg

        start = parts[j]
        end   = parts[j+1]
        rasterize_poly(shape.points[start:end], fg)

  for shape in data.shapes():
    rasterize_shape(shape, bg, fg)

  return img

img = rasterize("outlines/outlines_zurich_sp.shp", 4853, 3046, "white", "black")

img.save("chugach.png")
