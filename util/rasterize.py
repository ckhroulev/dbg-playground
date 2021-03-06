import shapefile
import Image, ImageDraw

def rasterize(filename, width, height, bg, fg, x_range = None, y_range = None):
    # Read in a shapefile
    data = shapefile.Reader(filename)

    img = Image.new("L", (width, height), bg)
    draw = ImageDraw.Draw(img)

    # domain width and height
    if x_range is None or y_range is None:
        x_min, x_max = data.bbox[0], data.bbox[2]
        y_min, y_max = data.bbox[1], data.bbox[3]
    else:
        x_min, x_max = x_range
        y_min, y_max = y_range

    xdist = x_max - x_min
    ydist = y_max - y_min

    xratio = width / xdist
    yratio = height / ydist

    def rasterize_poly(points, color):
        pixels = []
        for x,y in points:
            px = int((x - x_min) * xratio)
            py = int((y - y_min) * yratio)
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

if __name__ == "__main__":
    img = rasterize("outlines/outlines_zurich_sp.shp", 4853, 3046, "white", "black")

    img.save("chugach.png")
