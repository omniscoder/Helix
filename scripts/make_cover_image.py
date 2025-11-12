from PIL import Image, ImageDraw, ImageFont
import math, pathlib

def draw_helix(draw, width, height, n_turns=3):
    cx, cy = width/2, height/2
    amp = height/3
    points = []
    for x in range(width):
        y = cy + amp * math.sin((x/width)*n_turns*2*math.pi)
        points.append((x, y))
    draw.line(points, fill="#0d9488", width=4)
    # complementary strand
    points2 = [(x, y+10) for x,y in points]
    draw.line(points2, fill="#14b8a6", width=4)

def make_cover(path="docs/assets/cover.png"):
    W, H = 1200, 630  # social-share aspect ratio
    img = Image.new("RGB", (W, H), "#1e293b")
    draw = ImageDraw.Draw(img)
    draw_helix(draw, W, H)
    font = ImageFont.truetype("DejaVuSans-Bold.ttf", 72)
    sub  = ImageFont.truetype("DejaVuSans.ttf", 36)
    draw.text((60, 60), "Veri-Helix", fill="white", font=font)
    draw.text((60, 150), "Proves what it plots.", fill="#94f7ef", font=sub)
    img.save(path)
    print("Cover written to", pathlib.Path(path).resolve())

if __name__ == "__main__":
    make_cover()
