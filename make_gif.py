from PIL import Image
import os
import re

def create_gif(png_folder, output_gif, duration):
    # List all files in the directory
    png_files = [f for f in os.listdir(png_folder) if f.endswith('.png')]
    
    # Sort files based on the numerical part of the filename
    png_files.sort(key=lambda f: int(re.search(r'\d+', f).group()))
    
    # Create a list to store the image objects
    images = []

    for file in png_files:
        st = "verysmalleps"
        if file[:len(st)] != st:
            continue
        print(file)
        file_path = os.path.join(png_folder, file)
        images.append(Image.open(file_path))

    # Save the images as a GIF
    if images:
        images[0].save(output_gif, save_all=True, append_images=images[1:], duration=duration, loop=0)

output_gif = 'output_gifs/verysmalleps.gif'
duration = 100  # Duration between frames in milliseconds

create_gif("gif_pngs", output_gif, duration)